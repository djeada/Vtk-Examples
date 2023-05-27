import vtk


def update_sphere_size(obj, event):
    global sphere1
    value = obj.GetRepresentation().GetValue()
    sphere1.SetRadius(value)
    sphere1.Modified()
    render_window.Render()


# Create a render window, renderer, and interactor
render_window = vtk.vtkRenderWindow()
renderer = vtk.vtkRenderer()
render_window.AddRenderer(renderer)
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# Create two spheres with different colors
sphere1 = vtk.vtkSphereSource()
sphere1.SetCenter(-1.0, 0.0, 0.0)
sphere1.SetRadius(0.5)
sphere1.Update()

sphere2 = vtk.vtkSphereSource()
sphere2.SetCenter(1.0, 0.0, 0.0)
sphere2.SetRadius(0.5)
sphere2.Update()

# Create two mapper and actor for the spheres
mapper1 = vtk.vtkPolyDataMapper()
mapper1.SetInputConnection(sphere1.GetOutputPort())

actor1 = vtk.vtkActor()
actor1.SetMapper(mapper1)
actor1.GetProperty().SetColor(1.0, 0.0, 0.0)

mapper2 = vtk.vtkPolyDataMapper()
mapper2.SetInputConnection(sphere2.GetOutputPort())

actor2 = vtk.vtkActor()
actor2.SetMapper(mapper2)
actor2.GetProperty().SetColor(0.0, 0.0, 1.0)

# Add the actors to the renderer
renderer.AddActor(actor1)
renderer.AddActor(actor2)

# Create a slider widget
slider_rep = vtk.vtkSliderRepresentation2D()
slider_rep.SetMinimumValue(0.1)
slider_rep.SetMaximumValue(1.0)
slider_rep.SetValue(0.5)
slider_rep.SetSliderLength(0.05)
slider_rep.SetSliderWidth(0.02)
slider_rep.SetEndCapLength(0.01)
slider_rep.SetEndCapWidth(0.02)
slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_rep.GetPoint1Coordinate().SetValue(0.2, 0.1)
slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
slider_rep.GetPoint2Coordinate().SetValue(0.8, 0.1)

# Set the background color of the slider representation
slider_rep.GetSliderProperty().SetColor(0.2, 0.2, 0.2)
slider_rep.GetTubeProperty().SetColor(0.7, 0.7, 0.7)
slider_rep.GetCapProperty().SetColor(0.2, 0.2, 0.2)

slider_widget = vtk.vtkSliderWidget()
slider_widget.SetInteractor(interactor)
slider_widget.SetRepresentation(slider_rep)
slider_widget.AddObserver(vtk.vtkCommand.InteractionEvent, update_sphere_size)

# Set up the camera and background color
camera = renderer.GetActiveCamera()
camera.SetPosition(0, 0, 5)
camera.SetFocalPoint(0, 0, 0)
renderer.SetBackground(1.0, 1.0, 1.0)

# Start the interaction
interactor.Initialize()
render_window.Render()
slider_widget.On()
interactor.Start()
