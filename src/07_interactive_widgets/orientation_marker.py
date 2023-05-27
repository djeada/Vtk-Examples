import vtk

# Create two vtkPlaneSource
plane1 = vtk.vtkPlaneSource()
plane1.SetXResolution(10)
plane1.SetYResolution(10)

plane2 = vtk.vtkPlaneSource()
plane2.SetXResolution(10)
plane2.SetYResolution(10)

# Define two vtkPolyDataMappers and vtkActor
mapper1 = vtk.vtkPolyDataMapper()
mapper1.SetInputConnection(plane1.GetOutputPort())
actor1 = vtk.vtkActor()
actor1.SetMapper(mapper1)
actor1.GetProperty().SetColor(1, 0, 0)  # Set the color of plane1 to red

mapper2 = vtk.vtkPolyDataMapper()
mapper2.SetInputConnection(plane2.GetOutputPort())
actor2 = vtk.vtkActor()
actor2.SetMapper(mapper2)
actor2.GetProperty().SetColor(0, 0, 1)  # Set the color of plane2 to blue

# Define a vtkRenderer and add the two actors
renderer = vtk.vtkRenderer()
renderer.AddActor(actor1)
renderer.AddActor(actor2)

# Create a render window
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

# Create a render window interactor
render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)

# The callback function to modify the plane when the widget is manipulated
def plane_widget_callback(obj, event):
    transform = vtk.vtkTransform()
    obj.GetTransform(transform)
    matrix = transform.GetMatrix()
    obj.GetProp3D().SetUserMatrix(matrix)


# Add the custom box widget for the two planes
plane_widget1 = vtk.vtkBoxWidget()
plane_widget1.SetInteractor(render_window_interactor)
plane_widget1.SetProp3D(actor1)
plane_widget1.On()
plane_widget1.AddObserver("InteractionEvent", plane_widget_callback)

plane_widget2 = vtk.vtkBoxWidget()
plane_widget2.SetInteractor(render_window_interactor)
plane_widget2.SetProp3D(actor2)
plane_widget2.On()
plane_widget2.AddObserver("InteractionEvent", plane_widget_callback)

# Initialize and start the interactor
render_window_interactor.Initialize()
render_window.Render()
render_window_interactor.Start()
