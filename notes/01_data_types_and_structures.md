## Data Types and Structures

### vtkDataObject
* Base class for all data objects in VTK
* Stores data and metadata including:
  - Geometric information: This pertains to the shape, position, and size of the data elements.
  - Topological information: This refers to the connectivity of data elements, indicating how they are related or connected to each other. For instance, in a mesh, it can describe how points form lines or polygons.
* Common subclasses include:
  - vtkImageData
  - vtkRectilinearGrid
  - vtkStructuredGrid
  - vtkPolydata
  - vtkUnstructuredGrid

```
|
|--- vtkDataObject
     |
     |--- vtkDataSet
          |
          |--- vtkPointSet
          |    |
          |    |--- vtkPolyData
          |    |
          |    |--- vtkStructuredPoints (vtkImageData)
          |    |
          |    |--- vtkStructuredGrid
          |    |
          |    |--- vtkUnstructuredGrid
          |
          |--- vtkRectilinearGrid
          |
          |--- vtkCompositeDataSet
               |
               |--- vtkHierarchicalBoxDataSet
               |
               |--- vtkMultiBlockDataSet
               |
               |--- vtkHierarchicalDataSet
               |
               |--- vtkOverlappingAMR
               |
               |--- vtkNonOverlappingAMR
```

### Overview of VTK Data Structures

#### vtkImageData
- **Description**: Represents a regular grid with fixed topology and uniform spacing between points.
- **Applications**: Ideal for volumetric data like images or 3D scalar fields.
- **Examples**: CT and MRI scans, where vtkImageData stores pixel values at each grid point for 3D visualization.

#### vtkRectilinearGrid
- **Description**: A regular grid with fixed topology but non-uniform spacing between points.
- **Applications**: Suitable for data with varying resolution, like in climate or terrain elevation data.
- **Examples**: Climate models with varying altitude resolution, or terrain data with resolution changing with slope.

#### vtkPolyData
- **Description**: Represents a dataset comprising points, vertices, lines, polygons, and triangle strips.
- **Applications**: Ideal for surface meshes and 3D models.
- **Examples**: Models of 3D objects (e.g., vehicles, buildings), terrain surfaces (mountains, valleys).

#### vtkStructuredGrid
- **Description**: Curvilinear grid maintaining fixed topology.
- **Applications**: Best for data on curvilinear coordinate systems.
- **Examples**: Fluid flow simulations around objects, finite volume method simulations with a fixed cell division.

#### vtkUnstructuredGrid
- **Description**: An irregular grid with flexible topology, capable of containing various geometric primitives.
- **Applications**: Suitable for complex geometries and adaptive meshes.
- **Examples**: Finite element method simulations with complex, dynamic domains, models of intricate 3D geometries (human brain, aircraft wings).

#### Structured vs Unstructured Grids
- **Structured Grids**: Regular grids with fixed topology, easier to work with, more memory-efficient.
  - **Examples**: vtkImageData, vtkRectilinearGrid, vtkStructuredGrid.
- **Unstructured Grids**: Irregular grids with flexible topology, more versatile but less efficient.
  - **Example**: vtkUnstructuredGrid.

#### Multiblock Dataset
- **Description**: A collection of multiple datasets, organized hierarchically.
- **Applications**: Handles complex, hierarchical data structures effectively.
- **Examples**: Multi-domain simulations (each domain as a separate dataset), multi-resolution data (datasets at different detail levels).

A `vtkMultiBlockDataSet` organizes datasets (or blocks) hierarchically. Each block can be a vtkDataSet subclass or another composite dataset, allowing versatile dataset organization. For instance, a vtkMultiBlockDataSet might contain blocks like vtkPolyData, vtkStructuredGrid, and even another vtkMultiBlockDataSet. This structure is invaluable in large-scale simulations where data is segmented into blocks representing different simulation areas or system components, enabling VTK to manage complex, multi-part datasets coherently.

### Choosing the Right Data Structure
* Consider the geometry, topology, and resolution of your data
* The choice of data structure affects memory usage, processing time, and ease of use
* Some data structures can be converted to others, but not always without loss of information

| VTK Data Object       | Purpose and Use Cases                                                                 | Complexity                                                                                   | Flexibility                                                                                           | Efficiency                                                                                                  |
|-----------------------|---------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| PolyData              | Ideal for representing 3D surfaces and complex shapes like sculptures or terrain models. | Moderate: Comprised of polygons and polylines, allowing for manipulation and detailing.      | High: Excellent for complex, arbitrary shapes including intricate surface details.                   | High: Efficiently stores data specific to the shape, omitting unnecessary information.                      |
| Structured Grids      | Suitable for 3D grid-based data like volumetric images or regular spatial data.         | Low: Simple, regular structure with implicit connectivity, making it straightforward to handle. | Low: Confined to regular grid shapes, limiting its application to uniformly structured data.          | High for regular data: Efficiently handles uniformly spaced data but less suitable for irregular datasets. |
| Unstructured Grids    | Used for complex, irregular shapes like biological structures or geological formations. | High: Requires detailed definition of point connectivity, accommodating intricate geometries.   | High: Extremely versatile, capable of representing any shape, be it regular or irregular.             | Low to Moderate: Efficiency varies with the shape's complexity; generally less efficient than structured grids. |
| Structured Points (ImageData) | Perfect for regularly spaced data such as medical imaging (CT, MRI scans) or 3D textures. | Low: Simple and regular, with implicit connectivity for easy navigation.                       | Low: Restricted to regular grid shapes, ideal for uniformly spaced data.                             | Very High: Optimal for regular data thanks to implicit connectivity and uniform structure.                   |
| Rectilinear Grids     | Appropriate for data with variable resolution, like atmospheric or oceanic datasets.      | Moderate: Allows variable spacing while maintaining a regular grid structure.                  | Moderate: More adaptable than regular grids but less so than unstructured grids.                      | High: More efficient than unstructured grids for data that aligns with its variable but regular structure.  |
