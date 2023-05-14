## VTK Input and Output

VTK offers a comprehensive suite of tools for reading and writing a variety of data formats. This includes the native VTK file formats (legacy and XML-based), as well as numerous third-party formats.

### Common File Formats

VTK supports an extensive range of data formats, including:

- **Legacy VTK file format**: This format, ASCII or binary, supports various data structures and attributes. It comprises five sections: file version and identifier, header, type, structure, and data attributes.

- **XML-based VTK file format**: This format is more flexible and extensible than the legacy format. It is ASCII or binary and supports different data structures and attributes. Examples of file extensions include `.vtp` (PolyData), `.vtu` (UnstructuredGrid), and `.vts` (StructuredGrid).

- **Third-party file formats**: VTK can interface with several widely-used file formats in 3D graphics and scientific visualization, such as STL (commonly used for 3D printing), PLY (stores polygonal meshes), and OBJ (stores 3D models including geometry, textures, and materials).

### Reading and Writing Data

VTK uses subclasses of `vtkDataReader` and `vtkDataWriter` to read and write data objects, respectively. These subclasses cater to different data structures:

- **vtkPolyDataReader**, **vtkPolyDataWriter**
- **vtkStructuredPointsReader**, **vtkStructuredPointsWriter**
- **vtkStructuredGridReader**, **vtkStructuredGridWriter**
- **vtkUnstructuredGridReader**, **vtkUnstructuredGridWriter**

#### Example

Here is an example of reading an STL file and writing it as a VTK file:

```python
import vtk

# Create an STL reader
stl_reader = vtk.vtkSTLReader()
stl_reader.SetFileName("input.stl")

# Create a VTK writer
vtk_writer = vtk.vtkPolyDataWriter()
vtk_writer.SetInputConnection(stl_reader.GetOutputPort())
vtk_writer.SetFileName("output.vtk")

# Write the output file
vtk_writer.Write()
```

## Readers and Writers Comparison

A comparison of various readers and writers for different formats is provided below:

| Format    | Reader                     | Reader Output              | Writer                       | Writer Input             |
|-----------|----------------------------|----------------------------|------------------------------|--------------------------|
| STL       | vtkSTLReader               | vtkPolyData                | vtkSTLWriter                 | vtkPolyData              |
| OBJ       | vtkOBJReader               | vtkPolyData                | vtkOBJWriter                 | vtkPolyData              |
| VTK (old) | vtkUnstructuredGridReader  | vtkUnstructuredGrid        | vtkUnstructuredGridWriter    | vtkUnstructuredGrid      |
|           | vtkStructuredGridReader    | vtkStructuredGrid          | vtkStructuredGridWriter      | vtkStructuredGrid        |
|           | vtkPolyDataReader          | vtkPolyData                | vtkPolyDataWriter            | vtkPolyData              |
|           | vtkRectilinearGridReader   | vtkRectilinearGrid         | vtkRectilinearGridWriter     | vtkRectilinearGrid       |
|           | vtkStructuredPointsReader  | vtkStructuredPoints        | vtkStructuredPointsWriter    | vtkStructuredPoints      |
| VTU       | vtkXMLUnstructuredGridReader | vtkUnstructuredGrid      | vtkXMLUnstructuredGridWriter | vtkUnstructuredGrid      |
| VTM       | vtkXMLMultiBlockDataReader | vtkMultiBlockDataSet       | vtkXMLMultiBlockDataWriter   | vtkMultiBlockDataSet     |
| OpenFOAM  | vtkOpenFOAMReader          | vtkMultiBlockDataSet       | NA                           | NA                       |
| EnSight   | vtkEnSightGoldReader       | vtkMultiBlockDataSet       | NA                           | NA                       |
