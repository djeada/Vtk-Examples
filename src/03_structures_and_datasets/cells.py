"""
The cells in VTK represent a topology or a type of connectivity between points.
They do not hold any positional data of their own but refer to points which are separately stored.
They can be of various types: for example, a line (connecting two points), a triangle (connecting three points),
or a tetrahedron (connecting four points).
"""

import math

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_points_and_cell(cell_type, offset):
    points = vtk.vtkPoints()
    cell = None

    if cell_type == "triangle":
        points.InsertNextPoint(offset, 0, 0)
        points.InsertNextPoint(offset + 1, 0, 0)
        points.InsertNextPoint(offset + 0.5, 1, 0)
        cell = vtk.vtkTriangle()
        for i in range(3):
            cell.GetPointIds().SetId(i, i)

    elif cell_type == "quad":
        points.InsertNextPoint(offset + 2, 0, 0)
        points.InsertNextPoint(offset + 3, 0, 0)
        points.InsertNextPoint(offset + 3, 1, 0)
        points.InsertNextPoint(offset + 2, 1, 0)
        cell = vtk.vtkQuad()
        for i in range(4):
            cell.GetPointIds().SetId(i, i)

    elif cell_type == "line":
        points.InsertNextPoint(offset + 2, 0, 0)
        points.InsertNextPoint(offset + 7, 0, 0)
        cell = vtk.vtkLine()
        for i in range(2):
            cell.GetPointIds().SetId(i, i)

    elif cell_type == "vertex":
        points.InsertNextPoint(offset + 6, 0, 0)
        cell = vtk.vtkVertex()
        cell.GetPointIds().SetId(0, 0)

    elif cell_type == "polygon":
        for i in range(5):
            angle = 2 * math.pi * i / 5
            x = math.cos(angle) + offset + 7
            y = math.sin(angle)
            points.InsertNextPoint(x, y, 0)
        cell = vtk.vtkPolygon()
        cell.GetPointIds().SetNumberOfIds(5)
        for i in range(5):
            cell.GetPointIds().SetId(i, i)

    elif cell_type == "tetra":
        points.InsertNextPoint(offset + 9, 0, 0)
        points.InsertNextPoint(offset + 10, 0, 0)
        points.InsertNextPoint(offset + 9.5, 1, 0)
        points.InsertNextPoint(offset + 9.75, 0.5, 1)
        cell = vtk.vtkTetra()
        for i in range(4):
            cell.GetPointIds().SetId(i, i)

    return points, cell


def main():
    cell_types = ["triangle", "quad", "line", "vertex", "polygon", "tetra"]
    offsets = [0, 3, 6, 9, 12, 15]

    all_points = vtk.vtkPoints()
    all_cells = vtk.vtkCellArray()

    for cell_type, offset in zip(cell_types, offsets):
        points, cell = create_points_and_cell(cell_type, offset)
        point_ids = [
            all_points.InsertNextPoint(points.GetPoint(i))
            for i in range(points.GetNumberOfPoints())
        ]
        for i, pid in enumerate(point_ids):
            cell.GetPointIds().SetId(i, pid)
        all_cells.InsertNextCell(cell)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(all_points)

    # Set appropriate data for each cell type
    polydata.SetVerts(all_cells)  # For vertices
    polydata.SetLines(all_cells)  # For lines
    polydata.SetPolys(all_cells)  # For polygons

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()


if __name__ == "__main__":
    main()
