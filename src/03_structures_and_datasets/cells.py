"""
The cells in VTK represent a topology or a type of connectivity between points. They do not hold any positional data of their own but refer to points which are separately stored. They can be of various types: for example, a line (connecting two points), a triangle (connecting three points), or a tetrahedron (connecting four points).
"""
import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_points():
    """
    Create a set of points.
    """
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(0, 1, 0)
    points.InsertNextPoint(1, 1, 0)
    return points


def create_triangle():
    """
    Create a triangle cell.
    """
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 0)
    triangle.GetPointIds().SetId(1, 1)
    triangle.GetPointIds().SetId(2, 2)
    return triangle


def create_quad():
    """
    Create a quad cell.
    """
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, 0)
    quad.GetPointIds().SetId(1, 1)
    quad.GetPointIds().SetId(2, 3)
    quad.GetPointIds().SetId(3, 2)
    return quad


def create_cell_array(triangle, quad):
    """
    Create a cell array and add the cells to it.
    """
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(triangle)
    cells.InsertNextCell(quad)
    return cells


def create_polydata(points, cells):
    """
    Create polydata to hold the points and cells.
    """
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)
    return polydata


def extract_and_print_cell_info(polydata):
    """
    Extract and print information about each cell.
    """
    for i in range(polydata.GetNumberOfCells()):
        cell = polydata.GetCell(i)
        cell_type = cell.GetCellType()
        print(f"Cell {i} is of type {cell_type}:")
        for j in range(cell.GetNumberOfPoints()):
            point_id = cell.GetPointId(j)
            x, y, z = polydata.GetPoint(point_id)
            print(f"    Point {j}: ({x}, {y}, {z})")


def main():
    points = create_points()
    triangle = create_triangle()
    quad = create_quad()
    cells = create_cell_array(triangle, quad)
    polydata = create_polydata(points, cells)

    extract_and_print_cell_info(polydata)

    # Visualizing the cells
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    # Display the geometry
    pipeline = VisualisationPipeline(mappers=[mapper], edges_visible=True)
    pipeline.run()


if __name__ == "__main__":
    main()
