import vtk

# Create a sphere
sphereSource = vtk.vtkSphereSource()
sphereSource.SetCenter(0.0, 0.0, 0.0)
sphereSource.SetRadius(0.5)

# Create a mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphereSource.GetOutputPort())

# Create an actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a renderer
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)

# Create a render window
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# Create a text actor
textActor = vtk.vtkTextActor()
textActor.SetInput("This is a sphere")
textActor.GetTextProperty().SetFontSize(24)
textActor.GetTextProperty().SetColor(1.0, 0.0, 0.0)  # Red text
textActor.SetPosition(20, 20)

# Add the text actor to the renderer
renderer.AddActor2D(textActor)

# Create an interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

# Start the rendering loop
interactor.Initialize()
interactor.Start()
