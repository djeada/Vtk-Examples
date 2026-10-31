# End-to-end VTK case study: build, export, and validate a mesh

This small project connects the main ideas from the notes: generate geometry, produce triangles, attach a named scalar array, export to VTP, read the result back, validate it, and configure color mapping. The workflow uses **VTK and the Python standard library only**. It runs without opening a rendering window, making it useful in CI or on a headless workstation.

## Problem and expected output

Construct a sphere and use each point's Z coordinate as a scalar named `elevation`. The sphere source defaults to radius 0.5, so the scalar range should be approximately `[-0.5, 0.5]`. Save to `surface.vtp` and confirm that the file retains the same number of points, cells, and scalar tuples.

`vtkSphereSource → vtkTriangleFilter → elevation point array → vtkXMLPolyDataWriter → vtkXMLPolyDataReader → vtkPolyDataMapper`

## Complete, runnable script

Save this example as `case_study.py` and run `python case_study.py`. It writes into the current directory.

~~~python
from pathlib import Path
import vtk


def build_surface():
    source = vtk.vtkSphereSource()
    source.SetThetaResolution(48)
    source.SetPhiResolution(24)

    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(source.GetOutputPort())
    triangles.Update()
    mesh = triangles.GetOutput()

    elevation = vtk.vtkFloatArray()
    elevation.SetName("elevation")
    elevation.SetNumberOfTuples(mesh.GetNumberOfPoints())
    for point_id in range(mesh.GetNumberOfPoints()):
        _, _, z = mesh.GetPoint(point_id)
        elevation.SetValue(point_id, z)
    mesh.GetPointData().SetScalars(elevation)
    return mesh


def save_and_validate(mesh, filename):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(filename))
    writer.SetInputData(mesh)
    if writer.Write() != 1:
        raise RuntimeError("VTP export failed")

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(filename))
    reader.Update()
    loaded = reader.GetOutput()

    if loaded.GetNumberOfPoints() != mesh.GetNumberOfPoints():
        raise RuntimeError("Point count changed during export")
    if loaded.GetNumberOfCells() != mesh.GetNumberOfCells():
        raise RuntimeError("Cell count changed during export")
    values = loaded.GetPointData().GetArray("elevation")
    if values is None or values.GetNumberOfTuples() != loaded.GetNumberOfPoints():
        raise RuntimeError("The elevation array was not preserved")
    return loaded


def make_mapper(mesh):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetScalarModeToUsePointData()
    mapper.SelectColorArray("elevation")
    mapper.SetScalarRange(-0.5, 0.5)
    return mapper


if __name__ == "__main__":
    path = Path("surface.vtp")
    original = build_surface()
    restored = save_and_validate(original, path)
    mapper = make_mapper(restored)
    print("Saved:", path.resolve())
    print("Points:", restored.GetNumberOfPoints())
    print("Cells:", restored.GetNumberOfCells())
    print(
        "Elevation range:",
        restored.GetPointData().GetArray("elevation").GetRange()
    )
~~~

The script should print a positive number of points and cells and an elevation range close to `(-0.5, 0.5)`. You can inspect the resulting file in ParaView, another VTK application, or a text editor if you explicitly choose ASCII XML output.

## Understand the design decisions

| Step | Why it is necessary | What would go wrong otherwise |
| --- | --- | --- |
| Triangle conversion | Creates triangles understood by many mesh tools | Later filters may require different cell types |
| Point-data allocation | Gives each point exactly one scalar | Wrong tuple counts or associations can break coloring |
| Explicit array name | Allows filters/mappers to select the intended field | A mapper may use an unintended active array |
| VTP export | Preserves polygonal geometry and arrays | Surface-only formats like STL drop general attributes |
| Read-back assertions | Catch lost arrays and changed topology early | A successful write alone can conceal data loss |
| Fixed scalar range | Colors carry stable numeric meaning | Auto-ranging makes side-by-side comparison unreliable |

In this script, the mesh is updated before directly attaching the elevation array. The writer uses `SetInputData()` on this resulting object. For a pipeline whose fields must be recomputed as source parameters change, implement the array computation as a filter and connect the writer to its output port instead.

## Extend the case study

1. Replace elevation with the radial distance from the origin. Predict the scalar range *before* running the script.
2. Attach a value to every **cell** rather than every point. Change the mapper's scalar mode accordingly.
3. Increase sphere resolution. Compare triangle count, VTP file size, and export time; note what each metric measures.
4. Export STL instead of VTP. Reload it and check whether the elevation array survived.
5. Add a renderer and actor using the mapper, then give users a fixed camera and color legend; see [interactivity](https://github.com/djeada/Vtk-Examples/blob/main/notes/05_interactivity.md).

## Troubleshooting checklist

- **No file:** inspect working directory and `writer.Write()` return value.
- **Missing colors:** check scalar name, point/cell association, mapper scalar selection, and scalar range.
- **Unexpected geometry:** print `GetBounds()` and ensure source units/coordinates are correct.
- **Different data after reload:** compare types, cells, points, and arrays—not just a screenshot.
- **Slow processing:** profile source, filter, write, and render stages separately as discussed in [performance](https://github.com/djeada/Vtk-Examples/blob/main/notes/07_performance_optimization_and_parallelism.md).

This is an intentionally small validation target: it exercises the complete data path without conflating file correctness with GUI setup.
