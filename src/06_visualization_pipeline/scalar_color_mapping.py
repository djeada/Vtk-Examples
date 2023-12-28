import numpy as np
import vtk

# Create a 3D grid of points, each with a scalar value
dims = (20, 20, 20)
x = np.linspace(-1.0, 1.0, dims[0])
y = np.linspace(-1.0, 1.0, dims[1])
z = np.linspace(-1.0, 1.0, dims[2])
x, y, z = np.meshgrid(x, y, z, indexing="ij")

# Create a scalar value for each point in the grid
scalars = np.sqrt(x ** 2 + y ** 2 + z ** 2)

# Create a structured grid object and assign the points and scalars
grid = vtk.vtkStructuredGrid()
grid.SetDimensions(dims)

points = vtk.vtkPoints()
points.SetNumberOfPoints(np.prod(dims))

for i in range(dims[0]):
    for j in range(dims[1]):
        for k in range(dims[2]):
            points.SetPoint(
                i * dims[1] * dims[2] + j * dims[2] + k,
                x[i, j, k],
                y[i, j, k],
                z[i, j, k],
            )
grid.SetPoints(points)

scalarArray = vtk.vtkDoubleArray()
scalarArray.SetNumberOfTuples(np.prod(dims))
for i in range(np.prod(dims)):
    scalarArray.SetTuple(i, (scalars.ravel()[i],))
scalarArray.SetName("Scalars")
grid.GetPointData().SetScalars(scalarArray)

# Create a mapper and set the scalar range to the range of our scalar data
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(grid)
mapper.SetScalarModeToUsePointData()
mapper.SetScalarRange(scalars.min(), scalars.max())

# Create an actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a renderer, add the actor to it
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)

# Create a render window, add the renderer to it
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# Create an interactor, add the render window to it
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

# Initialize the interactor and start the rendering loop
interactor.Initialize()
interactor.Start()
