## Developing Custom Filters and Algorithms in VTK

Creating custom filters and algorithms opens up a world of possibilities for tailored data processing and visualization. By extending VTK’s capabilities, it becomes possible to carry out specialized techniques that meet the unique needs of scientific research, engineering, medical imaging, and data analysis workflows.

VTK comes with a broad range of built-in filters and classes that cover many common visualization tasks, but there are situations where more specific functionality is required. For instance, it may be necessary to process data from specialized scientific instruments, create a custom metric for point analysis, or experiment with novel geometry-manipulation algorithms. In those cases, developing a custom filter allows the following:

* Seamless integration of a custom algorithm into the native VTK pipeline.
* Reuse of existing VTK infrastructure for rendering, interaction, and I/O.
* Leveraging VTK’s performance-minded memory and execution model.
* Maintaining a consistent workflow without switching between external libraries.

What makes this especially valuable is that a custom filter, once written, behaves like any other VTK component: it can be connected to sources, chained with other filters, inspected, reused, and shared. Instead of being “one-off code,” it becomes a reusable building block that fits naturally into larger pipelines.

Below is a walkthrough of how VTK processes data, the steps for building a custom filter, and a complete example that computes Euclidean distance from each point in polygonal data to a specified target point.

### Understanding the VTK Pipeline

The VTK pipeline is the backbone of VTK’s data processing and visualization workflow. It is built around a demand-driven architecture, meaning computation typically happens only when something downstream, such as a mapper or renderer, requests output. This design keeps workflows modular and efficient while allowing complex visualization systems to remain understandable as a sequence of small, composable steps.

Here is a simple ASCII diagram illustrating the VTK pipeline:

```
[Source] -> [Filter 1] -> [Filter 2] -> ... -> [Mapper] -> [Actor] -> [Renderer]
```

* **Source**: Generates initial data, such as procedural geometry (e.g., a sphere or cube) or data read from a file (e.g., a volume dataset).
* **Filters**: Operate on the data to transform it or compute new information (e.g., smoothing, contour extraction, feature detection).
* **Mapper**: Converts processed data into a graphical representation the rendering engine can draw.
* **Actor**: Represents an object (geometry + visual properties) in the 3D scene.
* **Renderer**: Manages drawing the actor(s) onto the screen or into an off-screen buffer.

Creating a custom filter effectively inserts a new link into this chain. The filter accepts upstream data, performs specialized computations, and forwards either modified geometry or new derived attributes downstream. By following VTK conventions, especially around input/output ports and pipeline modification, custom filters remain interoperable with the rest of the ecosystem.

### Steps to Create a Custom Filter

Building a custom filter in VTK is a structured process that ensures the new component behaves consistently with the existing framework. The following outline describes a common development path.

#### I. Identify Data and Goals

Determine the data type to process and the computation or transformation to apply. For example, polygonal meshes (`vtkPolyData`) behave differently from image volumes (`vtkImageData`) or irregular finite-element meshes (`vtkUnstructuredGrid`). Choosing correctly at this stage simplifies implementation and improves performance.

#### II. Choose an Appropriate Base Class

VTK provides several base classes for filters, each tailored for different types of datasets. The most frequently used include:

| Base Class              | Purpose & Description                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------------- |
| `vtkAlgorithm`          | General-purpose algorithm class supporting multiple inputs/outputs and varying data types.     |
| `vtkPolyDataAlgorithm`  | Specialized for polygonal data such as meshes and surfaces; ideal for many 3D model workflows. |
| `vtkImageDataAlgorithm` | Specialized for image-based data, common in medical imaging and grid-based operations.         |

* For polygonal meshes, `vtkPolyDataAlgorithm` is typically the best fit.
* For pixel/voxel grids and image volumes, `vtkImageDataAlgorithm` provides conveniences specific to image extents and spacing.
* `vtkAlgorithm` is the most flexible option when supporting multiple data types or uncommon configurations.

#### III. Subclass the Chosen Base Class

In C++, this is typically done through header (`.h`) and implementation (`.cxx`) files. In Python, subclassing can be done directly. The subclass defines member variables (filter parameters) and overrides one or more key pipeline methods.

#### IV. Implement Core Methods

Most custom filters implement or override the following:

* **`RequestData()`**: The main execution method, where input is read, computation is performed, and output is written.
* **`RequestInformation()`** (optional): Used for meta-information such as extents, spacing, or requested output type (more common in image processing).
* **`RequestDataObject()`** (optional): Used when output type differs from default assumptions.

If a filter needs non-default port counts, it should configure them using `SetNumberOfInputPorts()` and `SetNumberOfOutputPorts()`.

#### V. Expose Parameters and Methods

Filters become practical when they are configurable. Parameters such as thresholds, coordinate inputs, and scaling factors should be implemented as setter/getter methods that follow VTK conventions (for example, `SetX()` and `GetX()`). When parameters change, `Modified()` should be called so downstream consumers know output must be recomputed.

#### VI. Compile and Use

* In C++, the filter is added to a project build system (often via CMake) and compiled as part of an executable or library.
* In Python, the filter can be imported and used directly once defined.

### Example: Creating a Custom Filter to Compute Point Distances

Consider polygonal data representing a 3D surface (a sphere, an imported CAD model, or an extracted anatomical mesh). A common analysis task is computing the distance from each vertex to a fixed target point. This can support use cases such as:

* **Computational geometry**: measuring deformation gradients or proximity fields.
* **Medical imaging**: examining how tissue points distribute around a reference.
* **Engineering**: checking tolerances relative to a design reference feature.

Once computed, distances are typically stored in a scalar array. This makes the results immediately usable for:

* Color mapping (visual inspection of near/far regions).
* Thresholding, contouring, and further derived-field computations.

#### Implementing the Custom Filter

The filter below subclasses `vtkPolyDataAlgorithm` and overrides `RequestData()`. The output is a shallow copy of the input `vtkPolyData` with an added point-data scalar array named `"DistanceToTarget"`.

```python
import math
import vtk


class DistanceToPointFilter(vtk.vtkPolyDataAlgorithm):
    """
    A custom filter that computes the Euclidean distance of each point in a vtkPolyData
    to a specified target point in 3D space. The distances are stored as a scalar array.

    Usage:
        distance_filter = DistanceToPointFilter()
        distance_filter.SetInputConnection(somePolyDataSource.GetOutputPort())
        distance_filter.SetTargetPoint(1.0, 2.0, 3.0)
        distance_filter.Update()
        output_polydata = distance_filter.GetOutput()
    """

    def __init__(self):
        super().__init__()
        self.TargetPoint = [0.0, 0.0, 0.0]

        # Explicit port configuration (often optional in Python, but clear and safe)
        self.SetNumberOfInputPorts(1)
        self.SetNumberOfOutputPorts(1)

    def SetTargetPoint(self, x, y, z):
        """Sets the 3D coordinates of the target point from which distances are calculated."""
        self.TargetPoint = [float(x), float(y), float(z)]
        self.Modified()

    def GetTargetPoint(self):
        """Returns the current target point as a tuple (x, y, z)."""
        return tuple(self.TargetPoint)

    def RequestData(self, request, inInfo, outInfo):
        """
        Main execution method: reads input vtkPolyData and produces output vtkPolyData
        with a scalar distance array named 'DistanceToTarget'.
        """
        input_data = vtk.vtkPolyData.GetData(inInfo[0], 0)
        output_data = vtk.vtkPolyData.GetData(outInfo, 0)

        if input_data is None or output_data is None:
            return 0

        # Preserve geometry/connectivity/attributes, then add a new scalar array
        output_data.ShallowCopy(input_data)

        num_points = input_data.GetNumberOfPoints()

        distances = vtk.vtkFloatArray()
        distances.SetName("DistanceToTarget")
        distances.SetNumberOfComponents(1)
        distances.SetNumberOfTuples(num_points)

        tx, ty, tz = self.TargetPoint
        for i in range(num_points):
            px, py, pz = input_data.GetPoint(i)
            dx = px - tx
            dy = py - ty
            dz = pz - tz
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            distances.SetValue(i, distance)

        output_data.GetPointData().AddArray(distances)
        output_data.GetPointData().SetScalars(distances)

        return 1
```

#### What Happens Inside `RequestData()`

* The input and output objects are retrieved from the VTK pipeline (`GetData`).
* `ShallowCopy` transfers geometry, connectivity, and existing attributes from input to output.
* A new `vtkFloatArray` is created to store distances.
* Each point’s Euclidean distance to the target point is computed and written into the array.
* The array is attached to the output’s point data and set as active scalars for coloring and downstream operations.

### Using the Custom Distance Filter

The following snippet demonstrates the filter in a simple pipeline using a sphere source. The computed distance scalars are then used for color mapping with a lookup table.

#### Creating and Configuring the Distance Filter

```python
import vtk

# Create a sphere source for demonstration
sphere = vtk.vtkSphereSource()
sphere.SetRadius(1.0)
sphere.SetThetaResolution(30)
sphere.SetPhiResolution(30)

# Create and configure the distance filter
distance_filter = DistanceToPointFilter()
distance_filter.SetInputConnection(sphere.GetOutputPort())
distance_filter.SetTargetPoint(2.0, 0.0, 0.0)  # Target point outside the sphere
distance_filter.Update()

output_data = distance_filter.GetOutput()
```

#### Setting Up the Visualization Pipeline with Color Mapping

```python
# Create mapper and configure scalar visualization
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(output_data)

# Set up color mapping
scalar_range = output_data.GetPointData().GetScalars().GetRange()
lut = vtk.vtkLookupTable()
lut.SetHueRange(0.667, 0.0)  # Blue to red color range
lut.SetTableRange(scalar_range)
lut.Build()

mapper.SetLookupTable(lut)
mapper.SetScalarRange(scalar_range)

# Create actor and set up visualization
actor = vtk.vtkActor()
actor.SetMapper(mapper)
```

With this setup, the output appears as a sphere colored by distance from the target point:

![Side by side comparison](https://github.com/user-attachments/assets/a7e0dd03-2a27-457f-bd9d-b500b61deac1)

* Left: Original sphere
* Right: Sphere colored by distance from target point (2.0, 0.0, 0.0)
* Blue: Points closer to the target
* Red: Points farther from the target
* Smooth gradient between based on actual distances
