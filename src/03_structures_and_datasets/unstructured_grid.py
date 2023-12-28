"""
This module demonstrates the creation, manipulation, and visualization of an unstructured grid using the Visualization Toolkit (VTK). An unstructured grid is a versatile data structure in VTK that allows for representing complex geometries composed of various types of cells. In this example, a hexahedron and a tetrahedron are created, along with their associated points, and then combined into a single unstructured grid. Additionally, scalar and vector fields are attached to this grid to illustrate the concept of data association with grid points and cells.

Key Components and Workflow:

1. Hexahedron and Tetrahedron Creation:
   - Hexahedron (vtkHexahedron) and tetrahedron (vtkTetra) cells are created to represent two different 3D geometric shapes.
   - The point ids for each cell are set to define the geometry of these cells.

2. Points Creation (vtkPoints):
   - A vtkPoints object is created to store the coordinates of the vertices for both the hexahedron and tetrahedron.

3. Unstructured Grid (vtkUnstructuredGrid):
   - An unstructured grid is defined to hold the created cells.
   - This grid type is chosen for its ability to handle diverse cell types within a single data structure.

4. Scalar and Vector Fields:
   - Scalar values (random) are attached to the points of the grid, demonstrating how to associate data with grid points.
   - Vector values (random) are attached to the cells of the grid, showing data association with grid cells.

5. Visualization:
   - A vtkDataSetMapper is used to map the unstructured grid data to graphical primitives.
   - The visualization pipeline, encapsulated in the 'VisualisationPipeline' class, is then utilized to render the grid with the scalar and vector fields. The pipeline is configured to highlight the edges of the grid and enhance the visual understanding of its structure.

"""

import numpy as np
import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_hexahedron():
    """
    Create and return a hexahedron cell.
    """
    hexahedron = vtk.vtkHexahedron()
    for i in range(8):
        hexahedron.GetPointIds().SetId(i, i)
    return hexahedron


def create_tetrahedron():
    """
    Create and return a tetrahedron cell.
    """
    tetrahedron = vtk.vtkTetra()
    for i in range(4):
        tetrahedron.GetPointIds().SetId(i, i + 8)
    return tetrahedron


def create_points():
    """
    Create and return vtkPoints for the hexahedron and tetrahedron.
    """
    points = vtk.vtkPoints()
    # Hexahedron points
    for i in range(8):
        points.InsertNextPoint(i % 2, (i // 2) % 2, i // 4)
    # Tetrahedron points
    points.InsertNextPoint(2, 2, 2)
    points.InsertNextPoint(2, 3, 2)
    points.InsertNextPoint(3, 2, 2)
    points.InsertNextPoint(2, 2, 3)
    return points


def create_unstructured_grid(points, hexahedron, tetrahedron):
    """
    Create and return an unstructured grid containing the hexahedron and tetrahedron.
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())
    ugrid.InsertNextCell(tetrahedron.GetCellType(), tetrahedron.GetPointIds())
    return ugrid


def attach_fields_to_grid(ugrid, points):
    """
    Attach scalar and vector fields to the unstructured grid.
    """
    # Attach a scalar field to the points
    scalars = vtk.vtkFloatArray()
    scalars.SetName("Scalars")
    for i in range(points.GetNumberOfPoints()):
        scalars.InsertNextValue(np.random.rand())
    ugrid.GetPointData().SetScalars(scalars)

    # Attach a vector field to the cells
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetName("Vectors")
    for i in range(ugrid.GetNumberOfCells()):
        vectors.InsertNextTuple3(np.random.rand(), np.random.rand(), np.random.rand())
    ugrid.GetCellData().SetVectors(vectors)


def visualize_unstructured_grid(ugrid):
    """
    Visualize the unstructured grid with scalar and vector fields.
    """
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(ugrid)
    mapper.SetScalarModeToUsePointData()
    mapper.SetColorModeToMapScalars()
    mapper.SelectColorArray("Scalars")

    pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
    pipeline.run()


def main():
    hexahedron = create_hexahedron()
    tetrahedron = create_tetrahedron()
    points = create_points()
    ugrid = create_unstructured_grid(points, hexahedron, tetrahedron)
    attach_fields_to_grid(ugrid, points)
    visualize_unstructured_grid(ugrid)


if __name__ == "__main__":
    main()
