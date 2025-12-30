## Input and Output

Input/output is where VTK starts being a tool you can actually plug into real workflows. Your data almost never starts life inside VTK, it comes from scanners, simulators, CAD tools, or research pipelines. So the ability to read reliably, preserve attributes, and write back out in the right format is what makes your visualizations reproducible and shareable.

VTK offers a comprehensive suite of tools for reading and writing a variety of data formats. This includes the native VTK file formats (legacy and XML-based), as well as numerous third-party formats.

A useful way to think about it: file formats are not just “containers,” they encode assumptions, what kind of dataset it is, how connectivity is represented, and which attributes are supported. When a file load “works” but looks wrong, it’s usually because something about those assumptions didn’t match your expectations.

**Do:** treat file I/O as part of your pipeline design, not a last-minute step.
**Don’t:** assume every format preserves every attribute (normals, scalars, materials, units, time steps, etc.) the same way.

### Common File Formats

VTK supports an extensive range of data formats, including:

Before diving into formats, here’s the big “why”: the *best* format depends on what you care about most, human readability, storage size, speed, compatibility, or rich metadata. That’s why VTK still supports both the old legacy format and the newer XML formats, plus popular third-party mesh standards.

**Do:** pick the format that matches your goal (debugging vs production vs interoperability).
**Don’t:** default to the first one that loads, your future self will thank you for choosing intentionally.

#### I. Legacy VTK File Format

* **Nature**: ASCII or binary.
* **Features**: Supports various data structures and attributes.
* **Structure**: Composed of five main sections - file version and identifier, header, dataset type, dataset structure, and data attributes.

Legacy VTK is the “classic” option: straightforward, widely supported, and easy to inspect if you choose ASCII. It’s great for learning, debugging, and quick exports. The tradeoff is that it’s less extensible and less structured than the XML family, and it can get bulky for large datasets.

**Do:** use legacy ASCII when you want to quickly sanity-check what’s inside a file.
**Don’t:** rely on it for very large datasets or long-term pipelines where extensibility and metadata matter.

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

#### II. XML-Based VTK File Format

* **Nature**: ASCII or binary, offering enhanced flexibility and extensibility.
* **Features**: Supports a diverse range of data structures and attributes.
* **File Extensions**: Includes formats like `.vtp` for PolyData, `.vtu` for UnstructuredGrid, and `.vts` for StructuredGrid, among others.

XML-based VTK formats are the “modern default” for many workflows because they’re more descriptive and extensible, and they map cleanly to specific dataset types through file extensions. That explicit mapping is a big deal: it reduces ambiguity and makes toolchains more reliable.

**Do:** prefer XML formats when you want robust, future-proof storage and clearer dataset typing.
**Don’t:** assume “XML” means “human friendly.” These files can still be large, and binary/appended data is common for performance.

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

#### III. Third-Party File Formats

* **Scope**: VTK interfaces seamlessly with many popular 3D graphics and scientific visualization file formats.
* **STL**: Predominantly used in 3D printing.
* **PLY**: Specializes in storing polygonal meshes.
* **OBJ**: Capable of holding complex 3D models encompassing geometry, texture, and material data.

Third-party formats are where interoperability lives. You use these when your data is coming from (or going to) other ecosystems, CAD tools, 3D modeling software, printing pipelines, or simulation packages. The catch is that each format has its own “personality.” Some are geometry-only, some carry materials, some are loose about units, and some ignore advanced attributes.

**Do:** choose third-party formats when you need compatibility with external tools.
**Don’t:** assume they preserve VTK-specific richness (like arbitrary point/cell arrays) unless you’ve verified it.

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

I/O in VTK is designed to feel consistent: you pick a reader or writer class for a format, set the filename, connect the pipeline, and execute. That consistency is the real quality-of-life feature, once you learn one reader/writer pair, you’ve basically learned the pattern for dozens.

There is a suite of subclasses derived from `vtkDataReader` and `vtkDataWriter`. These subclasses are specialized for handling various VTK data structures, emphasizing efficient and accurate data manipulation. The design ensures flexibility in reading and writing different types of data while maintaining the robustness of data integrity and format compatibility.

The important “why” behind the subclass approach: different dataset types need different handling. Writing `vtkPolyData` is not the same as writing `vtkImageData`, and mixing them up tends to fail loudly, or worse, appear to work while silently dropping information.

**Do:** match the reader/writer to the dataset type you expect.
**Don’t:** ignore the output type, knowing whether you got `vtkPolyData` vs `vtkUnstructuredGrid` determines what filters you can apply next.

#### Subclasses for Data Reading and Writing

Each subclass under `vtkDataReader` and `vtkDataWriter` is tailored for specific data structures, facilitating precise and optimized read/write operations:

* **vtkPolyDataReader** and **vtkPolyDataWriter**: For handling polygonal data, commonly used in 3D graphics and modeling.
* **vtkStructuredPointsReader** and **vtkStructuredPointsWriter**: Optimized for structured point datasets, where data is arranged in a regular grid.
* **vtkStructuredGridReader** and **vtkStructuredGridWriter**: Suitable for structured grid data, a step above structured points in complexity, allowing for non-uniform grids.
* **vtkUnstructuredGridReader** and **vtkUnstructuredGridWriter**: Designed for unstructured grid data, which is the most flexible, accommodating irregularly spaced data points.

#### Example

This example is a great “real life” scenario: you get something common and portable (STL), but you want to bring it into the VTK world so you can run filters, add attributes, and keep a pipeline consistent. Converting to a VTK-native format is often the first step in building a repeatable workflow.

**Do:** treat conversion scripts like this as reusable utilities in your toolbox.
**Don’t:** forget what STL represents: it’s primarily a surface triangle mesh, so you shouldn’t expect volumetric cells or rich per-point attributes unless you add them later.

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

Tables like this are useful because they prevent a super common beginner mistake: picking a reader/writer because the file extension looks right, without checking the dataset type it produces or expects. In VTK, dataset type is everything, filters, mappers, and pipeline expectations all hang on it.

**Do:** use the “Output Data Type / Input Data Type” columns as your reality check.
**Don’t:** expect every format to be symmetric (read + write). Some formats are read-only in practice, especially for complex third-party ecosystems.

| Format       | Reader Class                   | Output Data Type       | Writer Class                   | Input Data Type        |
| ------------ | ------------------------------ | ---------------------- | ------------------------------ | ---------------------- |
| STL          | `vtkSTLReader`                 | `vtkPolyData`          | `vtkSTLWriter`                 | `vtkPolyData`          |
| OBJ          | `vtkOBJReader`                 | `vtkPolyData`          | `vtkOBJWriter`                 | `vtkPolyData`          |
| VTK (Legacy) | `vtkUnstructuredGridReader`    | `vtkUnstructuredGrid`  | `vtkUnstructuredGridWriter`    | `vtkUnstructuredGrid`  |
|              | `vtkStructuredGridReader`      | `vtkStructuredGrid`    | `vtkStructuredGridWriter`      | `vtkStructuredGrid`    |
|              | `vtkPolyDataReader`            | `vtkPolyData`          | `vtkPolyDataWriter`            | `vtkPolyData`          |
|              | `vtkRectilinearGridReader`     | `vtkRectilinearGrid`   | `vtkRectilinearGridWriter`     | `vtkRectilinearGrid`   |
|              | `vtkStructuredPointsReader`    | `vtkStructuredPoints`  | `vtkStructuredPointsWriter`    | `vtkStructuredPoints`  |
| VTU          | `vtkXMLUnstructuredGridReader` | `vtkUnstructuredGrid`  | `vtkXMLUnstructuredGridWriter` | `vtkUnstructuredGrid`  |
| VTM          | `vtkXMLMultiBlockDataReader`   | `vtkMultiBlockDataSet` | `vtkXMLMultiBlockDataWriter`   | `vtkMultiBlockDataSet` |
| OpenFOAM     | `vtkOpenFOAMReader`            | `vtkMultiBlockDataSet` | N/A                            | N/A                    |
| EnSight      | `vtkEnSightGoldReader`         | `vtkMultiBlockDataSet` | N/A                            | N/A                    |
