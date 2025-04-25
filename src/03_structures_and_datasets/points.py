"""
n VTK, a point is a location in space. Points are used to define the geometry of graphical primitives. You can think of points as the vertices of graphical primitives such as lines, polygons, and volumes.

To store and manipulate points, VTK provides the vtkPoints class. This class stores an array of 3D points, and provides methods for inserting points, retrieving points, and other point-related operations.
"""

import numpy as np
import vtk
import vtk.util.numpy_support as vtk_np

from src.common.simple_pipeline import VisualisationPipeline


def create_vtk_points():
    """
    Create a vtkPoints object and insert some points.
    """
    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(0.0, 1.0, 0.0)
    return points


def print_points_info(points):
    """
    Print information about the vtkPoints object.
    """
    print("Number of points:", points.GetNumberOfPoints())
    for i in range(points.GetNumberOfPoints()):
        point = points.GetPoint(i)
        print(f"Point {i}: {point}")


def vtk_points_to_numpy(points):
    """
    Convert vtkPoints to a NumPy array.
    """
    return vtk_np.vtk_to_numpy(points.GetData())


def numpy_to_vtk_points(numpy_array):
    """
    Convert a NumPy array to vtkPoints.
    """
    vtk_points = vtk.vtkPoints()
    vtk_points.SetData(vtk_np.numpy_to_vtk(numpy_array))
    return vtk_points


def create_and_visualize_polydata(points):
    """
    Create polydata from vtkPoints and visualize them.
    """
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    glyphFilter = vtk.vtkVertexGlyphFilter()
    glyphFilter.SetInputData(polydata)
    glyphFilter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyphFilter.GetOutputPort())

    pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
    pipeline.run()


def main():
    vtk_points = create_vtk_points()
    print_points_info(vtk_points)

    numpy_points = vtk_points_to_numpy(vtk_points)
    print("Numpy array:\n", numpy_points)

    new_numpy_points = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
    new_vtk_points = numpy_to_vtk_points(new_numpy_points)
    print_points_info(new_vtk_points)

    create_and_visualize_polydata(vtk_points)


if __name__ == "__main__":
    main()
