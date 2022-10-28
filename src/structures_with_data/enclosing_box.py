from dataclasses import dataclass
import numpy as np
import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline


@dataclass
class Point:
    x: float = 0
    y: float = 0
    z: float = 0


box = Point(5.0, 5.0, 5.0)
density = 30
n_points = 100000


def f(x, y, z):
    return np.sin(x) * np.sin(y) * z ** 2


def generate_grid():
    grid = vtk.vtkImageData()
    grid.SetDimensions(density, density, density)
    grid.SetSpacing(box.x / density, box.y / density, box.z / density)
    return grid


def generate_data():
    data = vtk.vtkFloatArray()
    data.SetNumberOfValues(n_points)

    for i in range(density):
        z = box.z / density * i - box.z / 2.0
        for j in range(density):
            y = box.y / density * j - box.y / 2.0
            for k in range(density):
                x = box.x / density * k - box.x / 2.0
                n = k + j * density + i * density * density
                data.SetValue(n, f(x, y, z))

    return data


def generate_contour(grid):
    contour = vtk.vtkContourFilter()
    contour.SetInputData(grid)
    contour.SetValue(0, 0.1)
    contour.SetValue(1, 1.0)
    contour.SetValue(2, 5.0)
    contour.SetValue(3, 10.0)
    return contour


def generate_outline(contour):
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(contour.GetOutputPort())
    return outline


if __name__ == "__main__":
    grid = generate_grid()
    data = generate_data()
    contour = generate_contour(grid)
    outline = generate_outline(contour)
    grid.GetPointData().SetScalars(data)

    contour_mapper = vtk.vtkPolyDataMapper()
    contour_mapper.SetInputConnection(contour.GetOutputPort())
    contour_mapper.UseLookupTableScalarRangeOn()

    outline_mapper = vtk.vtkPolyDataMapper()
    outline_mapper.SetInputConnection(outline.GetOutputPort())

    pipeline = VisualisationPipeline([outline_mapper, contour_mapper])
    pipeline.run()
