## Developing Custom Filters and Algorithms in VTK

Creating custom filters and algorithms in the Visualization Toolkit (VTK) opens up a world of possibilities for tailored data processing and visualization. By extending VTK's capabilities, you can carry out specialized techniques that meet the unique needs of your projects—whether it's for scientific research, engineering, medical imaging, or data analysis. 

VTK comes with a broad range of built-in filters and classes that cover many common visualization tasks, but there may be occasions when you need more specific functionality. For instance, you might need to process data from specialized scientific instruments, create a custom metric for point analysis, or experiment with novel geometry-manipulation algorithms. In these cases, developing a custom filter allows you to:

- Easily integrate your algorithm into the native VTK pipeline.
- Reuse existing VTK infrastructure for rendering, interaction, and I/O.
- Use VTK’s optimized performance and memory management.
- Keep your workflow consistent without needing to switch between external libraries.

Below, we’ll walk through the fundamentals of how VTK processes data, the steps for building a custom filter, and a detailed example that demonstrates how you might carry out a distance-to-point calculation for polygonal data.

### Understanding the VTK Pipeline

The VTK pipeline is the backbone of VTK’s data processing and visualization workflow. It’s designed around a demand-driven architecture, meaning that data flows through the pipeline when something downstream (like a renderer) requests it. The pipeline consists of sources, filters, and mappers, connected in a sequence where each filter takes data from its predecessor, processes it, and passes it on to the next stage.

Here’s a simple ASCII diagram illustrating the VTK pipeline:

```
[Source] -> [Filter 1] -> [Filter 2] -> ... -> [Mapper] -> [Actor] -> [Renderer]
```

- Generates initial data, such as procedural geometry (e.g., a sphere or cube) or data read from a file (e.g., a volume dataset).
- Operate on the data to transform or compute new information (e.g., smoothing, contour extraction, feature detection).
- Converts the processed data into a graphical representation that the rendering engine can understand.
- Represents an object (geometry + properties) in the 3D scene.
- Manages rendering the actor(s) onto the screen or into an off-screen buffer.

When you create a custom filter, you add a new link in this pipeline. Your filter can accept data from upstream VTK filters or sources, perform specialized calculations, and pass new or modified data downstream. By following VTK’s design patterns, your custom filter remains interoperable with any other VTK component, preserving the modular, pipeline-oriented architecture that makes VTK powerful and flexible.

### Steps to Create a Custom Filter

Building a custom filter in VTK is a multi-step process that ensures your new filter integrates seamlessly with the existing framework. Below is an outline of the typical steps you’ll follow:

I. Identify Your Data and Goals  

Determine the type of data you want to process and the computations or transformations you intend to apply. For example, are you dealing with polygonal data (`vtkPolyData`), image data (`vtkImageData`), or unstructured grids (`vtkUnstructuredGrid`)? Identifying this will help you choose the appropriate base class.

II. Choose an Appropriate Base Class  

VTK provides several base classes for filters, each tailored for different data types. The most commonly used base classes are:

| Base Class             | Purpose & Description                                                                                  |
|------------------------|---------------------------------------------------------------------------------------------------------|
| `vtkAlgorithm`         | A general-purpose class capable of handling multiple inputs and outputs of various data types.          |
| `vtkPolyDataAlgorithm` | Specialized for polygonal data, such as meshes or surfaces, making it ideal for 3D models.              |
| `vtkImageDataAlgorithm`| Specialized for image-based data, commonly used in medical imaging or 2D/3D grid-based data operations. |

- If your algorithm needs to process polygonal meshes (e.g., STL files, surfaces), `vtkPolyDataAlgorithm` is usually the go-to choice.  
- For image data (like DICOM or 2D images), `vtkImageDataAlgorithm` provides convenient methods for working with pixel/voxel grids.  
- `vtkAlgorithm` is the most general and flexible if you have an unusual data structure or want to support multiple data types.

III. Subclass the Chosen Base Class  

In languages like C++, you’ll create a header (`.h`) and implementation (`.cxx`) file, then subclass the base class. In Python, you can just subclass directly. Your subclass should define any member variables you need (like parameters for your filter) and override the relevant methods.

IV. Carry out Core Methods  

Every custom filter has a few important methods you need to carry out (or override), with the most important typically being:

- Where the main processing occurs. This is where you read from the input data object, execute your custom algorithm, and populate the output data object.
- (Optional) Used to provide meta-information about the data, like extent or data type. More relevant for image-based filters.
- (Optional) Used to specify the exact data type of your output if it differs from the default.

Properly managing input and output ports is necessary. You need to specify how many inputs your filter expects (`SetNumberOfInputPorts`) and how many outputs it will provide (`SetNumberOfOutputPorts`) if you deviate from the defaults.

V. Expose Parameters and Methods  

If your filter has parameters (e.g., a threshold value, a coordinate, or a scaling factor), create setter and getter methods so users can configure your filter. Keep them in sync with the VTK naming conventions (e.g., `SetX`, `GetX`, etc.) if possible.

VI. Compile and Use  

- You’ll typically add your custom filter files to your project, update your CMakeLists, and compile them into a library or as part of your executable.
- You can use it directly once you’ve defined the subclass in your script or module. Just import your custom class and insert it into the pipeline.

### Example: Creating a Custom Filter to Compute Point Distances

Let’s say you have a 3D model (e.g., a sphere, a mesh from a CT scan, or a CAD object), and you need to calculate how far each vertex on the model is from a specific point in 3D space. This is a common operation in fields like:

- Computational Geometry – where you might want to find gradient fields or measure deformation.  
- Medical Imaging – to see how tissue or tumor points are distributed around a reference point.  
- Engineering – to measure distance from a design feature to a reference location, or to measure design tolerances.

When the distance for each point is computed, it’s often stored in a scalar array so it can be used for:

- Color mapping (to visually inspect where distances are small or large).  
- Further calculations, such as thresholding or contouring the distance field.

#### Implementing the Custom Filter

Below is an example demonstrating how to carry out this custom distance filter for polygonal data. We’ll subclass `vtkPolyDataAlgorithm` and override the `RequestData()` method, where we do all our distance computations.

```python
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
    outputPolyData = distance_filter.GetOutput()
"""

def __init__(self):
    super().__init__()

    # By default, VTK doesn't have automatic new-style class initialization in Python.
    # We'll make sure we set up the input and output ports properly.
    # Note: VTK may handle input/output ports automatically; in some cases,
    # you might need to call SetNumberOfInputPorts(1) / SetNumberOfOutputPorts(1).
    # Initialize the target point. Users can modify this via SetTargetPoint().

    self.TargetPoint = [0.0, 0.0, 0.0]

def SetTargetPoint(self, x, y, z):
    """
    Sets the 3D coordinates of the target point from which distances will be calculated.
    """
    self.TargetPoint = [float(x), float(y), float(z)]

    # Mark the pipeline as modified so changes propagate

    self.Modified()

def GetTargetPoint(self):
    """
    Returns the current target point as a tuple (x, y, z).
    """
    return tuple(self.TargetPoint)

def RequestData(self, request, inInfo, outInfo):
    """
    The main execution method where the filter processes input vtkPolyData
    and produces an output with a scalar distance array.
    """

    # 1. Retrieve the input and output data objects from the pipeline.
    input_data = vtk.vtkPolyData.GetData(inInfo[0])
    output_data = vtk.vtkPolyData.GetData(outInfo)

    # 2. Copy the input to the output. 
    output_data.ShallowCopy(input_data)

    # 3. Get the number of points in the polygonal data.
    num_points = input_data.GetNumberOfPoints()

    # 4. Create a new vtkFloatArray to store the distance values for each point.
    distances = vtk.vtkFloatArray()
    distances.SetName("DistanceToTarget")
    distances.SetNumberOfComponents(1)
    distances.SetNumberOfTuples(num_points)

    # 5. Compute distances for each point and store them in the array.
    tx, ty, tz = self.TargetPoint
    for i in range(num_points):
        px, py, pz = input_data.GetPoint(i)
        dx = px - tx
        dy = py - ty
        dz = pz - tz
        distance = (dx2 + dy2 + dz2)  0.5
        distances.SetValue(i, distance)

    # 6. Attach the array to the output's point data. Also set it as the active scalar.
    output_data.GetPointData().AddArray(distances)
    output_data.GetPointData().SetScalars(distances)

    return 1

# Optionally, you could override RequestInformation() if you need
# to specify extents or other meta-information.

```

Here’s a breakdown of what’s happening in `RequestData()`:

- In the `GetData` step, the `input_data` and `output_data` are fetched from the pipeline.  
- The `ShallowCopy` operation preserves geometry, connectivity, and attributes by copying the input to the output.  
- A `ShallowCopy` shares underlying data and is suitable unless a fully independent copy is required.  
- A `DeepCopy` creates a fully independent copy but is not needed if only new arrays are added.  
- A float array is created during distance calculation to store Euclidean distances.  
- Distances are computed by iterating over each point and calculating its distance to the target point.  
- The calculated distance array is attached to the output as part of the point data.  
- The distance array is set as the active scalar to enable use in coloring or further pipeline processing.

#### Using the Custom Distance Filter

This custom VTK filter calculates the Euclidean distance from each point in a vtkPolyData to a specified target point. Below is an example demonstrating how to use it with a sphere visualization.

**Creating and configuring the distance filter:**

```python
# Create a sphere source for demonstration
sphere = vtk.vtkSphereSource()
sphere.SetRadius(1.0)
sphere.SetThetaResolution(30)
sphere.SetPhiResolution(30)
sphere.Update()

# Create and configure the distance filter
filter = DistanceToPointFilter()
filter.SetInputData(sphere.GetOutput())
filter.SetTargetPoint(2.0, 0.0, 0.0)  # Target point outside the sphere
output_data = filter.ProcessDataObject(sphere.GetOutput())
```

**Setting up the visualization pipeline with color mapping:**

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

The output will look similar to the following:

![Side by side comparison](https://github.com/user-attachments/assets/a7e0dd03-2a27-457f-bd9d-b500b61deac1)

- Left: Original sphere  
- Right: Sphere colored by distance from target point (2.0, 0.0, 0.0)
- Blue: Points closer to the target
- Red: Points farther from the target
- Smooth gradient between based on actual distances
- Computes Euclidean distance from each point to a target point
- Stores distances as scalar values in the output's point data
- Supports visualization with customizable color mapping
- Works with any vtkPolyData input
- Adjust the target point location to highlight different distance patterns
- Modify the lookup table's hue range for different color schemes
- Use the orientation widget to better understand spatial relationships
