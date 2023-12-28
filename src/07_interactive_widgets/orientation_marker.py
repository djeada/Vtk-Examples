import vtk
from vtkmodules.vtkCommonTransforms import vtkTransform
from vtkmodules.vtkFiltersSources import vtkPlaneSource
from vtkmodules.vtkInteractionWidgets import vtkBoxWidget
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)


def create_plane_actor(color: tuple) -> vtk.vtkActor:
    """
    Create a plane actor with specified color.

    Args:
        color (tuple): A tuple of three float values representing RGB color.

    Returns:
        vtk.vtkActor: An actor representing a colored plane.
    """
    plane = vtkPlaneSource()
    plane.SetXResolution(10)
    plane.SetYResolution(10)

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(plane.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)

    return actor


def plane_widget_callback(obj, event):
    """
    Callback function to modify the plane when the widget is manipulated.

    Args:
        obj: The object associated with the callback.
        event: The event triggering the callback.
    """
    transform = vtkTransform()
    obj.GetTransform(transform)
    matrix = transform.GetMatrix()
    obj.GetProp3D().SetUserMatrix(matrix)


def setup_interaction(
    plane_actor: vtk.vtkActor, render_window_interactor: vtk.vtkRenderWindowInteractor
):
    """
    Set up interaction widget for a plane actor.

    Args:
        plane_actor (vtk.vtkActor): The plane actor to interact with.
        render_window_interactor (vtk.vtkRenderWindowInteractor): The render window interactor.
    """
    plane_widget = vtkBoxWidget()
    plane_widget.SetInteractor(render_window_interactor)
    plane_widget.SetProp3D(plane_actor)
    plane_widget.On()
    plane_widget.AddObserver("InteractionEvent", plane_widget_callback)


# Main execution
if __name__ == "__main__":
    # Create two plane actors with different colors
    red_plane_actor = create_plane_actor((1, 0, 0))
    blue_plane_actor = create_plane_actor((0, 0, 1))

    # Define a vtkRenderer and add the two actors
    renderer = vtkRenderer()
    renderer.AddActor(red_plane_actor)
    renderer.AddActor(blue_plane_actor)

    # Create and set up render window and interactor
    render_window = vtkRenderWindow()
    render_window.AddRenderer(renderer)

    render_window_interactor = vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Set up interaction for both planes
    setup_interaction(red_plane_actor, render_window_interactor)
    setup_interaction(blue_plane_actor, render_window_interactor)

    # Initialize and start the interactor
    render_window_interactor.Initialize()
    render_window.Render()
    render_window_interactor.Start()
