## VTK Input and Output

### Overview
* VTK provides functionality for reading and writing various data formats
* Common file formats include:
  - Legacy VTK file format
  - XML-based VTK file format
  - Third-party file formats (e.g., STL, PLY, OBJ)

### vtkDataReader and vtkDataWriter
* Base classes for reading and writing data objects in VTK
* Subclasses include:
  - vtkPolyDataReader, vtkPolyDataWriter
  - vtkStructuredPointsReader, vtkStructuredPointsWriter
  - vtkStructuredGridReader, vtkStructuredGridWriter
  - vtkUnstructuredGridReader, vtkUnstructuredGridWriter

### Legacy VTK file format
* ASCII or binary format
* Supports various data structures and attributes
* Consists of five parts: file version and identifier, header, type, structure, and attributes

#### Example
Legacy VTK file format structure

```
vtk DataFile Version x.x

a title (up to 256 characters, terminated by a new line)
ASCII or BINARY
DATASET type (e.g., STRUCTURED_POINTS, STRUCTURED_GRID, UNSTRUCTURED_GRID, POLYDATA, RECTILINEAR_GRID)
Additional information related to the dataset type (dimensions, spacing, coordinates)
```

### XML-based VTK file format
* ASCII or binary format
* More flexible and extensible than the legacy format
* Supports various data structures and attributes
* Example file extensions: .vtp (PolyData), .vtu (UnstructuredGrid), .vts (StructuredGrid)

### Third-party file formats
* VTK can read and write several popular file formats used in 3D graphics and scientific visualization
* Examples:
  - STL: widely used for 3D printing and computer-aided manufacturing
  - PLY: stores polygonal meshes with vertex, face, and edge data
  - OBJ: stores 3D models with geometry, textures, and materials

## Example: Reading an STL file and writing it as a VTK file

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
