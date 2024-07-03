import vtk
import numpy as np


def generate_flower_mesh(petal_count: int, points_per_petal: int) -> vtk.vtkPolyData:
    theta = np.linspace(0, 2 * np.pi, points_per_petal)
    radius = np.cos(petal_count * theta) + 2

    points = vtk.vtkPoints()
    for i in range(points_per_petal):
        x = radius[i] * np.cos(theta[i])
        y = radius[i] * np.sin(theta[i])
        points.InsertNextPoint(x, y, 0)

    points.InsertNextPoint(0, 0, 0)

    polygons = vtk.vtkCellArray()
    for i in range(points_per_petal - 1):
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(3)
        polygon.GetPointIds().SetId(0, i)
        polygon.GetPointIds().SetId(1, i + 1)
        polygon.GetPointIds().SetId(2, points_per_petal)
        polygons.InsertNextCell(polygon)

    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(3)
    polygon.GetPointIds().SetId(0, points_per_petal - 1)
    polygon.GetPointIds().SetId(1, 0)
    polygon.GetPointIds().SetId(2, points_per_petal)
    polygons.InsertNextCell(polygon)

    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points)
    poly_data.SetPolys(polygons)

    return poly_data


def update_flower_mesh(petal_count: int, points_per_petal: int, actor: vtk.vtkActor) -> None:
    new_mesh = generate_flower_mesh(petal_count, points_per_petal)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(new_mesh)
    actor.SetMapper(mapper)


def slider_callback(
    obj: vtk.vtkObject,
    event: str,
    petal_slider: vtk.vtkSliderWidget,
    points_slider: vtk.vtkSliderWidget,
    actor: vtk.vtkActor,
    render_window: vtk.vtkRenderWindow,
    text_actor: vtk.vtkTextActor
) -> None:
    petal_count = int(petal_slider.GetRepresentation().GetValue())
    points_per_petal = int(points_slider.GetRepresentation().GetValue())
    update_flower_mesh(petal_count, points_per_petal, actor)
    text_actor.SetInput(f"Petals: {petal_count}\nPoints: {points_per_petal}")
    render_window.Render()


def create_slider_widget(
    title: str,
    interactor: vtk.vtkRenderWindowInteractor,
    min_value: float,
    max_value: float,
    initial_value: float,
    position: list[float]
) -> vtk.vtkSliderWidget:
    slider_rep = vtk.vtkSliderRepresentation2D()
    slider_rep.SetMinimumValue(min_value)
    slider_rep.SetMaximumValue(max_value)
    slider_rep.SetValue(initial_value)
    slider_rep.SetTitleText(title)
    slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint1Coordinate().SetValue(position[0], position[1])
    slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint2Coordinate().SetValue(position[2], position[3])
    slider_widget = vtk.vtkSliderWidget()
    slider_widget.SetInteractor(interactor)
    slider_widget.SetRepresentation(slider_rep)
    slider_widget.SetAnimationModeToAnimate()
    slider_widget.EnabledOn()
    return slider_widget


def rotate_flower(obj: vtk.vtkObject, event: str) -> None:
    actor.RotateZ(1)
    obj.GetRenderWindow().Render()


if __name__ == "__main__":
    petal_count = 8
    points_per_petal = 80

    flower_mesh = generate_flower_mesh(petal_count, points_per_petal)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(flower_mesh)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetEdgeVisibility(True)
    actor.GetProperty().SetEdgeColor(0, 0, 0)
    actor.GetProperty().SetColor(1, 1, 1)  # Light-colored flower

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    renderer.AddActor(actor)
    renderer.SetBackground(0.1, 0.1, 0.1)  # Dark background

    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(f"Petals: {petal_count}\nPoints: {points_per_petal}")
    text_actor.GetTextProperty().SetColor(1, 1, 1)  # White text
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.SetPosition2(10, 40)
    text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    text_actor.SetPosition(0.75, 0.85)
    renderer.AddActor2D(text_actor)

    petal_slider = create_slider_widget(
        "Petals", render_window_interactor, 1, 20, petal_count,
        [0.1, 0.1, 0.4, 0.1]
    )
    points_slider = create_slider_widget(
        "Points", render_window_interactor, 10, 200, points_per_petal,
        [0.1, 0.2, 0.4, 0.2]
    )

    slider_callback_closure = lambda obj, event: slider_callback(
        obj, event, petal_slider, points_slider, actor, render_window,
        text_actor
    )
    petal_slider.AddObserver("EndInteractionEvent", slider_callback_closure)
    points_slider.AddObserver("EndInteractionEvent", slider_callback_closure)

    render_window_interactor.AddObserver("TimerEvent", rotate_flower)
    render_window_interactor.CreateRepeatingTimer(100)

    render_window.Render()
    render_window_interactor.Start()
