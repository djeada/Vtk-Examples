## Input and Output

VTK offers a comprehensive suite of tools for reading and writing a variety of data formats. This includes the native VTK file formats (legacy and XML-based), as well as numerous third-party formats.

### Common File Formats

VTK supports an extensive range of data formats, including:

I. Legacy VTK File Format

- **Nature**: ASCII or binary.
- **Features**: Supports various data structures and attributes.
- **Structure**: Composed of five main sections - file version and identifier, header, dataset type, dataset structure, and data attributes.

Example:

```
# vtk DataFile Version 3.0
VTK Example Data
ASCII
DATASET POLYDATA
POINTS 8 float
0.0 0.0 0.0
1.0 0.0 0.0
1.0 1.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
1.0 0.0 1.0
1.0 1.0 1.0
0.0 1.0 1.0
POLYGONS 6 30
4 0 1 2 3
...
```

II. XML-Based VTK File Format

- **Nature**: ASCII or binary, offering enhanced flexibility and extensibility.
- **Features**: Supports a diverse range of data structures and attributes.
- **File Extensions**: Includes formats like `.vtp` for PolyData, `.vtu` for UnstructuredGrid, and `.vts` for StructuredGrid, among others.

Example:

```
<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="8" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="6">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 ...
        </DataArray>
      </Points>
      <Polys>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1 2 3 ...
        </DataArray>
        ...
      </Polys>
    </Piece>
  </PolyData>
</VTKFile>
```

III. Third-Party File Formats

- **Scope**: VTK interfaces seamlessly with many popular 3D graphics and scientific visualization file formats.
- **STL**: Predominantly used in 3D printing.
- **PLY**: Specializes in storing polygonal meshes.
- **OBJ**: Capable of holding complex 3D models encompassing geometry, texture, and material data.

Example (STL):

```
solid vtkGenerated
  facet normal 0 0 -1
    outer loop
      vertex 0.0 0.0 0.0
      vertex 1.0 0.0 0.0
      vertex 1.0 1.0 0.0
    endloop
  endfacet
  ...
endsolid vtkGenerated
```

Example (PLY):

```
ply
format ascii 1.0
element vertex 8
property float x
property float y
property float z
element face 6
property list uchar int vertex_indices
end_header
0.0 0.0 0.0
1.0 0.0 0.0
...
3 0 1 2
...
```

### Reading and Writing Data

There is a suite of subclasses derived from `vtkDataReader` and `vtkDataWriter`. These subclasses are specialized for handling various VTK data structures, emphasizing efficient and accurate data manipulation. The design ensures flexibility in reading and writing different types of data while maintaining the robustness of data integrity and format compatibility.

#### Subclasses for Data Reading and Writing

Each subclass under `vtkDataReader` and `vtkDataWriter` is tailored for specific data structures, facilitating precise and optimized read/write operations:

- **vtkPolyDataReader** and **vtkPolyDataWriter**: For handling polygonal data, commonly used in 3D graphics and modeling.
- **vtkStructuredPointsReader** and **vtkStructuredPointsWriter**: Optimized for structured point datasets, where data is arranged in a regular grid.
- **vtkStructuredGridReader** and **vtkStructuredGridWriter**: Suitable for structured grid data, a step above structured points in complexity, allowing for non-uniform grids.
- **vtkUnstructuredGridReader** and **vtkUnstructuredGridWriter**: Designed for unstructured grid data, which is the most flexible, accommodating irregularly spaced data points.

#### Example

Below is a Python script demonstrating how to read data from an STL file (common in 3D printing and modeling) and write it into VTK's native format.

```python
import vtk

# Initialize an STL reader and set the file to read
stl_reader = vtk.vtkSTLReader()
stl_reader.SetFileName("input.stl")

# Set up a VTK writer and connect it to the output of the STL reader
vtk_writer = vtk.vtkPolyDataWriter()
vtk_writer.SetInputConnection(stl_reader.GetOutputPort())
vtk_writer.SetFileName("output.vtk")

# Execute the writing process to convert the STL file to a VTK file
vtk_writer.Write()
```

This script shows the straightforward approach of VTK in converting data between different formats, highlighting its powerful data processing capabilities.

## Readers and Writers Comparison

A comparison of various readers and writers for different formats is provided below:

| Format    | Reader Class                    | Output Data Type            | Writer Class                  | Input Data Type          |
|-----------|---------------------------------|-----------------------------|-------------------------------|--------------------------|
| STL       | `vtkSTLReader`                  | `vtkPolyData`               | `vtkSTLWriter`                | `vtkPolyData`            |
| OBJ       | `vtkOBJReader`                  | `vtkPolyData`               | `vtkOBJWriter`                | `vtkPolyData`            |
| VTK (Legacy) | `vtkUnstructuredGridReader`   | `vtkUnstructuredGrid`       | `vtkUnstructuredGridWriter`   | `vtkUnstructuredGrid`    |
|           | `vtkStructuredGridReader`       | `vtkStructuredGrid`         | `vtkStructuredGridWriter`     | `vtkStructuredGrid`      |
|           | `vtkPolyDataReader`             | `vtkPolyData`               | `vtkPolyDataWriter`           | `vtkPolyData`            |
|           | `vtkRectilinearGridReader`      | `vtkRectilinearGrid`        | `vtkRectilinearGridWriter`    | `vtkRectilinearGrid`     |
|           | `vtkStructuredPointsReader`     | `vtkStructuredPoints`       | `vtkStructuredPointsWriter`   | `vtkStructuredPoints`    |
| VTU       | `vtkXMLUnstructuredGridReader`  | `vtkUnstructuredGrid`       | `vtkXMLUnstructuredGridWriter`| `vtkUnstructuredGrid`    |
| VTM       | `vtkXMLMultiBlockDataReader`    | `vtkMultiBlockDataSet`      | `vtkXMLMultiBlockDataWriter`  | `vtkMultiBlockDataSet`   |
| OpenFOAM  | `vtkOpenFOAMReader`             | `vtkMultiBlockDataSet`      | N/A                           | N/A                      |
| EnSight   | `vtkEnSightGoldReader`          | `vtkMultiBlockDataSet`      | N/A                           | N/A                      |
