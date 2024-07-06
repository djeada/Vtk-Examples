import vtk
import numpy as np

def naca0012(x):
    """ Calculate the NACA 0012 airfoil shape. """
    m = 0
    p = 0
    t = 0.12
    c = 1.0
    y_t = 5 * t * c * (0.2969 * np.sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c)**2 + 0.2843 * (x/c)**3 - 0.1015 * (x/c)**4)
    return y_t

# Generate airfoil points
x_coords = np.linspace(0, 1, 100)
y_coords = naca0012(x_coords)

# Create VTK points
points = vtk.vtkPoints()
for x, y in zip(x_coords, y_coords):
    points.InsertNextPoint(x, y, 0)

# Create a polygon
polygon = vtk.vtkPolygon()
polygon.GetPointIds().SetNumberOfIds(len(x_coords))
for i in range(len(x_coords)):
    polygon.GetPointIds().SetId(i, i)

# Create a cell array to store the polygon
polygons = vtk.vtkCellArray()
polygons.InsertNextCell(polygon)

# Create a polydata object
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetPolys(polygons)

# Write the polydata to a file
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName("naca0012.vtp")
writer.SetInputData(polydata)
writer.Write()

# Now, visualize the airfoil
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)

actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)
renderer.SetBackground(1, 1, 1)

renderWindow.Render()
renderWindowInteractor.Start()
