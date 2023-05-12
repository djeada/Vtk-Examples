"""
The cells in VTK represent a topology or a type of connectivity between points. They do not hold any positional data of their own but refer to points which are separately stored. They can be of various types: for example, a line (connecting two points), a triangle (connecting three points), or a tetrahedron (connecting four points).
"""
import vtk

from src.simple_pipeline import VisualisationPipeline

# Create the points
points = vtk.vtkPoints()
points.InsertNextPoint(0, 0, 0)  # Point 0
points.InsertNextPoint(1, 0, 0)  # Point 1
points.InsertNextPoint(0, 1, 0)  # Point 2
points.InsertNextPoint(1, 1, 0)  # Point 3

# Create a triangle cell
triangle = vtk.vtkTriangle()
triangle.GetPointIds().SetId(0, 0)  # the first point of the triangle is point 0
triangle.GetPointIds().SetId(1, 1)  # the second point is point 1
triangle.GetPointIds().SetId(2, 2)  # the third point is point 2

# Create a quad cell
quad = vtk.vtkQuad()
quad.GetPointIds().SetId(0, 0)  # the first point of the quad is point 0
quad.GetPointIds().SetId(1, 1)  # the second point is point 1
quad.GetPointIds().SetId(2, 3)  # the third point is point 3
quad.GetPointIds().SetId(3, 2)  # the fourth point is point 2

# Create a cell array and add the cells to it
cells = vtk.vtkCellArray()
cells.InsertNextCell(triangle)
cells.InsertNextCell(quad)

# Create polydata to hold the points and cells
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetPolys(cells)

# Now let's extract the cells
for i in range(polydata.GetNumberOfCells()):
    cell = polydata.GetCell(i)
    cell_type = cell.GetCellType()
    print(f"Cell {i} is of type {cell_type}:")
    for j in range(cell.GetNumberOfPoints()):
        point_id = cell.GetPointId(j)
        x, y, z = polydata.GetPoint(point_id)
        print(f"    Point {j}: ({x}, {y}, {z})")

# Visualizing the cells
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)

# Display the squares
pipeline = VisualisationPipeline(mappers=[mapper], edges_visible=True)
pipeline.run()
