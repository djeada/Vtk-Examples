import vtk
from vtkmodules.vtkRenderingCore import vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, vtkActor, vtkPolyDataMapper
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkInteractionWidgets import vtkSliderWidget
from vtkmodules.vtkCommonCore import vtkCommand
from vtkmodules.vtkRenderingAnnotation import vtkSliderRepresentation2D

def create_sphere_actor(center: tuple, radius: float, color: tuple) -> vtkActor:
    """
    Create a sphere actor with specified center, radius, and color.

    Args:
        center (tuple): Center coordinates of the sphere (x, y, z).
        radius (float): Radius of the sphere.
        color (tuple): RGB color of the sphere.

    Returns:
        vtkActor: An actor representing a colored sphere.
    """
    sphere_source = vtkSphereSource()
    sphere_source.SetCenter(*center)
    sphere_source.SetRadius(radius)
    sphere_source.Update()

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(sphere_source.GetOutputPort())

    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    return actor, sphere_source

def setup_slider_widget(interactor: vtkRenderWindowInteractor, callback_func) -> vtkSliderWidget:
    """
    Set up a slider widget for interaction.

    Args:
        interactor (vtkRenderWindowInteractor): The render window interactor.
        callback_func: The callback function for the slider interaction.

    Returns:
        vtkSliderWidget: Configured slider widget.
    """
    slider_rep = vtkSliderRepresentation2D()
    slider_rep.SetMinimumValue(0.1)
    slider_rep.SetMaximumValue(1.0)
    slider_rep.SetValue(0.5)
    slider_rep.SetSliderLength(0.05)
    slider_rep.SetSliderWidth(0.02)
    slider_rep.SetEndCapLength(0.01)
    slider_rep.SetEndCapWidth(0.02)
    slider_rep.GetSliderProperty().SetColor(0.2, 0.2, 0.2)
    slider_rep.GetTubeProperty().SetColor(0.7, 0.7, 0.7)
    slider_rep.GetCapProperty().SetColor(0.2, 0.2, 0.2)
    slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint1Coordinate().SetValue(0.2, 0.1)
    slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint2Coordinate().SetValue(0.8, 0.1)

    slider_widget = vtkSliderWidget()
    slider_widget.SetInteractor(interactor)
    slider_widget.SetRepresentation(slider_rep)
    slider_widget.AddObserver(vtkCommand.InteractionEvent, callback_func)

    return slider_widget

def update_sphere_size(obj, event):
    """
    Callback function to update the sphere size when the slider is moved.

    Args:
        obj: The object associated with the callback.
        event: The event triggering the callback.
    """
    value = obj.GetRepresentation().GetValue()
    sphere1.SetRadius(value)
    sphere1.Modified()
    render_window.Render()

if __name__ == "__main__":
    render_window = vtkRenderWindow()
    renderer = vtkRenderer()
    render_window.AddRenderer(renderer)
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Create two spheres with different colors
    actor1, sphere1 = create_sphere_actor((-1.0, 0.0, 0.0), 0.5, (1.0, 0.0, 0.0))
    actor2, _ = create_sphere_actor((1.0, 0.0, 0.0), 0.5, (0.0, 0.0, 1.0))

    # Add the actors to the renderer
    renderer.AddActor(actor1)
    renderer.AddActor(actor2)

    # Create and configure a slider widget
    slider_widget = setup_slider_widget(interactor, update_sphere_size)

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
