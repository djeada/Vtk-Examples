import numpy as np
import vtk

# Define the number of seed points
NUMBER_OF_SEED_POINTS = 30

# Create a grid of points.
points = vtk.vtkPoints()
for i in np.linspace(-10, 10, 21):
    for j in np.linspace(-10, 10, 21):
        for k in np.linspace(-10, 10, 21):
            points.InsertNextPoint(i, j, k)

# Create vectors at each point to represent a dipole field.
vectors = vtk.vtkDoubleArray()
vectors.SetNumberOfComponents(3)
vectors.SetNumberOfTuples(points.GetNumberOfPoints())
for i in range(points.GetNumberOfPoints()):
    x, y, z = points.GetPoint(i)
    r1 = (
        np.sqrt((x - 1) ** 2 + y ** 2 + z ** 2) + 1e-2
    )  # Add a small number to avoid division by zero.
    r2 = (
        np.sqrt((x + 1) ** 2 + y ** 2 + z ** 2) + 1e-2
    )  # Add a small number to avoid division by zero.
    vectors.SetTuple3(
        i,
        (x - 1) / r1 ** 3 - (x + 1) / r2 ** 3,
        y / r1 ** 3 - y / r2 ** 3,
        z / r1 ** 3 - z / r2 ** 3,
    )  # Dipole field.

# Combine the points and vectors into a structured grid.
grid = vtk.vtkStructuredGrid()
grid.SetDimensions(21, 21, 21)
grid.SetPoints(points)
grid.GetPointData().SetVectors(vectors)

# Seed the streamlines at the points along the x-axis.
seeds = vtk.vtkPointSource()
seeds.SetCenter(0, 0, 0)
seeds.SetRadius(10)  # The radius of the spherical shell where the seeds are placed.
seeds.SetNumberOfPoints(NUMBER_OF_SEED_POINTS)

# Create the streamlines.
streamer = vtk.vtkStreamTracer()
streamer.SetInputData(grid)
streamer.SetSourceConnection(seeds.GetOutputPort())
streamer.SetMaximumPropagation(500)  # Increase this to let the streamlines go further.
streamer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
streamer.SetInitialIntegrationStep(0.1)  # Increased step size
streamer.SetIntegrationDirectionToBoth()
streamer.SetComputeVorticity(False)

# Visualize the streamlines.
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(streamer.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(0, 0, 0)

# Create a glyph source (a sphere)
glyph_source = vtk.vtkSphereSource()
glyph_source.SetRadius(0.5)

# Define the points where the sources are located
source_points = vtk.vtkPoints()
source_points.InsertNextPoint(-5, 0, 0)  # location of the first source
source_points.InsertNextPoint(5, 0, 0)  # location of the second source

# Create a polydata object
source_polydata = vtk.vtkPolyData()
source_polydata.SetPoints(source_points)

render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)
interactor.Initialize()
interactor.Start()
