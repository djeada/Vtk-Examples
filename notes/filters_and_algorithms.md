## Filters and Algorithms

### Overview
* VTK filters and algorithms are used to process and manipulate data
* They can generate new data, extract features, or transform existing data
* Common types include:
  - Sources
  - Geometric filters
  - Topological filters
  - Scalars and attribute filters
  - Temporal filters

### vtkAlgorithm
* Base class for all VTK algorithms
* Subclasses include:
  - Source
    * Procedural source
    * Reader source
  - Filter

### Data Flow
* Source -> Data object
* Data object -> Filter -> Data object

### Sources
* Generate data objects or read data from files
* Examples: vtkSphereSource, vtkConeSource, vtkSTLReader, vtkXMLPolyDataReader

### Geometric Filters
* Modify the geometry of data objects
* Examples:
  - vtkShrinkFilter: shrinks the geometry of a dataset
  - vtkSmoothPolyDataFilter: smooths the surface of a polydata object
  - vtkDecimatePro: reduces the number of triangles in a mesh

### Topological Filters
* Modify the topology of data objects
* Examples:
  - vtkTriangleFilter: converts polygons to triangles
  - vtkDelaunay2D: constructs a 2D Delaunay triangulation
  - vtkContourFilter: generates contours (isosurfaces) from scalar values

### Scalars and Attribute Filters
* Modify or generate data attributes, such as scalars, vectors, or tensors
* Examples:
  - vtkGradientFilter: computes the gradient of a scalar field
  - vtkVectorNorm: computes the magnitude of vector data
  - vtkCurvatures: computes the Gaussian and mean curvatures of a surface

### Temporal Filters
* Process time-varying data or generate animations
* Examples:
  - vtkTemporalInterpolator: interpolates data between time steps
  - vtkTemporalShiftScale: shifts and scales time values
  - vtkTemporalStatistics: computes statistical information over time

## Example: Creating a Sphere Source and Applying a Shrink Filter

```python
import vtk

# Create a sphere source
sphere_source = vtk.vtkSphereSource()
sphere_source.SetRadius(1.0)

# Create a shrink filter
shrink_filter = vtk.vtkShrinkFilter()
shrink_filter.SetInputConnection(sphere_source.GetOutputPort())
shrink_filter.SetShrinkFactor(0.8)

# Update the filter to generate the output
shrink_filter.Update()
```
