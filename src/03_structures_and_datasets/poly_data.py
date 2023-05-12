import vtk

# Create three points.
point1 = [0.0, 0.0, 0.0]
point2 = [1.0, 0.0, 0.0]
point3 = [0.5, 0.866, 0.0]

# Create a vtkPoints object and store the points in it.
points = vtk.vtkPoints()
points.InsertNextPoint(point1)
points.InsertNextPoint(point2)
points.InsertNextPoint(point3)

# Create a triangle.
triangle = vtk.vtkTriangle()
triangle.GetPointIds().SetId(0, 0)
triangle.GetPointIds().SetId(1, 1)
triangle.GetPointIds().SetId(2, 2)

# Create a cell array to store the triangle.
triangles = vtk.vtkCellArray()
triangles.InsertNextCell(triangle)

# Create a polydata object.
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetPolys(triangles)

# Create a mapper and an actor.
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)

actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a renderer and a window.
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)

window = vtk.vtkRenderWindow()
window.AddRenderer(renderer)

# Create an interactor.
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

# Start the visualization.
interactor.Initialize()
interactor.Start()
