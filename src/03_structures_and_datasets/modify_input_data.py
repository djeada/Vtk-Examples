import random

import numpy as np
from pyvtk import *
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from vtkmodules.vtkIOLegacy import (
    vtkPolyDataReader,
    vtkStructuredPointsReader,
    vtkStructuredGridReader,
    vtkRectilinearGridReader,
)
from vtkmodules.vtkRenderingCore import vtkDataSetMapper

from src.io.read_vtk import read_vtk
from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "example"

width = 6
height = 4
length = 3


def f(x, y, z):
    return abs(x - width / 2) ** 2


def extend_mappers(output, mappers):
    if hasattr(output, "SetOrigin"):
        output.SetOrigin(-(1 + len(mappers)) * (width), 0, 0)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    mappers.append(mapper)


if __name__ == "__main__":
    mappers = []

    vtk_data = VtkData(
        StructuredPoints([width, height, length]),
        PointData(Scalars(np.random.randint(20, size=72))),
    )
    vtk_data.tofile(file_name, "ascii")
    modified_file_name = f"{file_name}.vtk"
    output = read_vtk(modified_file_name, vtkStructuredPointsReader)
    extend_mappers(output, mappers)

    vtk_data.point_data.append(vtk_data.structure.Scalars(f, "x*y*z"))
    modified_file_name = f"{file_name}_sp.vtk"
    vtk_data.tofile(modified_file_name)
    output = read_vtk(modified_file_name, vtkStructuredPointsReader)
    extend_mappers(output, mappers)

    points = [
        (i, j, k) for k in range(length) for j in range(height) for i in range(width)
    ]
    new_vtk_data = VtkData(StructuredGrid([width, height, length], points))
    # extend vtk_data with new data
    vtk_data.structure.points.clear()

    for i in range(len(new_vtk_data.structure.points)):
        vtk_data.structure.points.append(new_vtk_data.structure.points[i])

    modified_file_name = f"{file_name}_sg.vtk"
    vtk_data.tofile(modified_file_name)
    output = read_vtk(modified_file_name, vtkStructuredPointsReader)
    extend_mappers(output, mappers)

    vtk_data = VtkData(RectilinearGrid(range(width), range(height), range(length)))
    vtk_data.point_data.append(vtk_data.structure.Scalars(f, "x*y*z"))
    modified_file_name = f"{file_name}_rg.vtk"
    vtk_data.tofile(modified_file_name)
    output = read_vtk(modified_file_name, vtkRectilinearGridReader)
    extend_mappers(output, mappers)

    pipeline = VisualisationPipeline(mappers)
    pipeline.run()
