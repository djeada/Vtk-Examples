import vtk

# Create a cube
cubeSource = vtk.vtkCubeSource()

# Create a mapper and set the cube source as its input
cubeMapper = vtk.vtkPolyDataMapper()
cubeMapper.SetInputConnection(cubeSource.GetOutputPort())

# Create an actor and set the mapper as its input
cubeActor = vtk.vtkActor()
cubeActor.SetMapper(cubeMapper)

# Create a renderer and add the cube actor to it
renderer = vtk.vtkRenderer()
renderer.AddActor(cubeActor)

# Create a render window and add the renderer to it
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# Create an interactor and set the render window as its input
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

# Initialize the interactor
interactor.Initialize()

# The main animation loop
for i in range(360):
    # Rotate the cube actor
    cubeActor.RotateY(1)

    # Render the scene
    renderWindow.Render()

# Start the interaction
interactor.Start()
