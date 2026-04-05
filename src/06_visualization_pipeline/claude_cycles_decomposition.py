"""
Interactive visualization of Knuth's odd-m Hamiltonian cycle on Z_m^3.

Paper:
https://www-cs-faculty.stanford.edu/~knuth/papers/claude-cycles.pdf

This example focuses on the first Hamiltonian cycle described on page 3 of
Knuth's note for odd m >= 3:
  - spheres are placed on the m x m x m lattice
  - a timer advances through the Hamiltonian cycle
  - visited vertices are dimmed
  - the current vertex is highlighted
  - a path tube shows the traversed prefix of the cycle

Controls:
  - Play / Pause / Reset buttons
  - Slider "m" to rebuild the lattice with odd sizes only
  - Slider "speed (ms)" to control the timer interval
  - Space toggles play/pause
  - R resets the animation

Examples:
  python src/07_interactive_widgets/knuth_hamiltonian_cycle.py
  python src/07_interactive_widgets/knuth_hamiltonian_cycle.py --m 7 --speed-ms 80
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.button import Button2D


BACKGROUND_COLOR = (0.08, 0.09, 0.12)
TEXT_COLOR = (0.96, 0.97, 0.98)
PANEL_COLOR = (0.05, 0.06, 0.09)
BASE_POINT_COLOR = (0.48, 0.50, 0.56)
VISITED_POINT_COLOR = (0.34, 0.78, 0.98)
CURRENT_POINT_COLOR = (1.00, 0.82, 0.22)
PATH_COLOR = (0.98, 0.54, 0.18)
BOX_COLOR = (0.36, 0.39, 0.46)

EXPECTED_M3_SEQUENCE = """
022 002 000 001 011 012 010 020 021
121 101 111 112 122 102 100 110 120
220 221 201 202 200 210 211 212 222
""".split()


def validate_odd_m(value: int) -> int:
    if value < 3 or value % 2 == 0:
        raise ValueError("m must be odd and >= 3")
    return value


def nearest_odd(value: float, minimum: int, maximum: int) -> int:
    candidate = int(round((value - minimum) / 2.0)) * 2 + minimum
    candidate = max(minimum, min(maximum, candidate))
    if candidate % 2 == 0:
        candidate = candidate + 1 if candidate < maximum else candidate - 1
    return candidate


def cycle_step(vertex: tuple[int, int, int], m: int) -> tuple[int, int, int]:
    """
    First Hamiltonian cycle from Knuth's proof for odd m.

    Let s = (i + j + k) mod m.
      - s == 0:     bump i if j == m - 1, else bump k
      - 0 < s < m-1: bump k if i == m - 1, else bump j
      - s == m-1:  bump k if i == 0, else bump j
    """
    i, j, k = vertex
    s = (i + j + k) % m

    if s == 0:
        return ((i + 1) % m, j, k) if j == m - 1 else (i, j, (k + 1) % m)
    if s == m - 1:
        return (i, j, (k + 1) % m) if i == 0 else (i, (j + 1) % m, k)
    return (i, j, (k + 1) % m) if i == m - 1 else (i, (j + 1) % m, k)


def build_cycle_order(m: int) -> list[tuple[int, int, int]]:
    validate_odd_m(m)

    start = (0, m - 1, 2)
    order = [start]
    seen = {start}
    current = start

    for _ in range(m**3 - 1):
        current = cycle_step(current, m)
        if current in seen:
            raise ValueError(f"cycle repeats early for m={m} at {current}")
        order.append(current)
        seen.add(current)

    if cycle_step(current, m) != start:
        raise ValueError(f"cycle does not close for m={m}")

    if m == 3:
        expected = [tuple(int(ch) for ch in token) for token in EXPECTED_M3_SEQUENCE]
        if order != expected:
            raise ValueError("m=3 cycle does not match Knuth's printed example")

    return order


def lattice_position(i: int, j: int, k: int, m: int, spacing: float = 1.0) -> tuple[float, float, float]:
    offset = (m - 1) * 0.5
    return (
        spacing * (i - offset),
        spacing * (j - offset),
        spacing * (k - offset),
    )


def make_text_actor(
    x: float,
    y: float,
    font_size: int = 18,
    bold: bool = True,
) -> vtk.vtkTextActor:
    actor = vtk.vtkTextActor()
    prop = actor.GetTextProperty()
    prop.SetColor(*TEXT_COLOR)
    prop.SetFontSize(font_size)
    prop.SetBold(bold)
    prop.SetLineSpacing(1.15)
    prop.SetShadow(False)
    prop.SetJustificationToLeft()
    prop.SetVerticalJustificationToTop()
    actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    actor.SetPosition(x, y)
    return actor


def make_panel_actor(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    fill_color: tuple[float, float, float],
    opacity: float = 0.88,
) -> vtk.vtkActor2D:
    points = vtk.vtkPoints()
    points.InsertNextPoint(x0, y0, 0.0)
    points.InsertNextPoint(x1, y0, 0.0)
    points.InsertNextPoint(x1, y1, 0.0)
    points.InsertNextPoint(x0, y1, 0.0)

    quad = vtk.vtkCellArray()
    quad.InsertNextCell(4)
    for point_id in range(4):
        quad.InsertCellPoint(point_id)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(quad)

    coordinate = vtk.vtkCoordinate()
    coordinate.SetCoordinateSystemToNormalizedViewport()

    mapper = vtk.vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(coordinate)

    actor = vtk.vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*fill_color)
    actor.GetProperty().SetOpacity(opacity)

    return actor


def make_slider(
    interactor: vtk.vtkRenderWindowInteractor,
    minimum: float,
    maximum: float,
    value: float,
    title: str,
    p1: tuple[float, float],
    p2: tuple[float, float],
) -> vtk.vtkSliderWidget:
    rep = vtk.vtkSliderRepresentation2D()
    rep.SetMinimumValue(minimum)
    rep.SetMaximumValue(maximum)
    rep.SetValue(value)
    rep.SetTitleText(title)
    rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    rep.GetPoint1Coordinate().SetValue(*p1)
    rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    rep.GetPoint2Coordinate().SetValue(*p2)
    rep.SetSliderLength(0.02)
    rep.SetSliderWidth(0.03)
    rep.SetTubeWidth(0.005)
    rep.SetLabelFormat("%.0f")
    rep.SetTitleHeight(0.018)
    rep.SetLabelHeight(0.014)

    rep.GetTitleProperty().SetColor(*TEXT_COLOR)
    rep.GetTitleProperty().SetBold(True)
    rep.GetLabelProperty().SetColor(*TEXT_COLOR)
    rep.GetTubeProperty().SetColor(0.58, 0.62, 0.70)
    rep.GetSliderProperty().SetColor(0.96, 0.97, 0.98)
    rep.GetCapProperty().SetColor(0.96, 0.97, 0.98)

    widget = vtk.vtkSliderWidget()
    widget.SetInteractor(interactor)
    widget.SetRepresentation(rep)
    widget.SetAnimationModeToAnimate()
    widget.EnabledOn()
    return widget


def make_glyph_actor(
    polydata: vtk.vtkPolyData,
    radius: float,
    color: tuple[float, float, float] | None = None,
    opacity: float = 1.0,
) -> tuple[vtk.vtkGlyph3D, vtk.vtkActor]:
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(18)
    sphere.SetPhiResolution(18)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToDataScalingOff()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(opacity)
    if color is not None:
        actor.GetProperty().SetColor(*color)
    return glyph, actor


class KnuthHamiltonianCycleApp:
    def __init__(self, initial_m: int = 5, initial_speed_ms: int = 120, max_m: int = 15):
        self.max_m = max(3, max_m if max_m % 2 == 1 else max_m - 1)
        self.m = nearest_odd(initial_m, 3, self.max_m)
        self.speed_ms = max(20, int(initial_speed_ms))

        self.order: list[tuple[int, int, int]] = []
        self.vertex_positions: dict[tuple[int, int, int], tuple[float, float, float]] = {}
        self.completed_steps = 0
        self.complete = False
        self.playing = False
        self.timer_id: int | None = None

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(*BACKGROUND_COLOR)

        self.render_window = vtk.vtkRenderWindow()
        self.render_window.SetSize(1280, 900)
        self.render_window.AddRenderer(self.renderer)

        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)

        self.base_polydata = vtk.vtkPolyData()
        self.visited_polydata = vtk.vtkPolyData()
        self.current_polydata = vtk.vtkPolyData()
        self.path_polydata = vtk.vtkPolyData()

        _, self.base_actor = make_glyph_actor(
            self.base_polydata,
            radius=0.12,
            color=BASE_POINT_COLOR,
            opacity=0.18,
        )
        _, self.visited_actor = make_glyph_actor(
            self.visited_polydata,
            radius=0.12,
            color=VISITED_POINT_COLOR,
            opacity=0.96,
        )
        _, self.current_actor = make_glyph_actor(
            self.current_polydata,
            radius=0.12,
            color=CURRENT_POINT_COLOR,
            opacity=1.0,
        )

        path_tube = vtk.vtkTubeFilter()
        path_tube.SetInputData(self.path_polydata)
        path_tube.SetRadius(0.03)
        path_tube.SetNumberOfSides(16)
        path_tube.CappingOn()

        path_mapper = vtk.vtkPolyDataMapper()
        path_mapper.SetInputConnection(path_tube.GetOutputPort())

        self.path_actor = vtk.vtkActor()
        self.path_actor.SetMapper(path_mapper)
        self.path_actor.GetProperty().SetColor(*PATH_COLOR)
        self.path_actor.GetProperty().SetOpacity(0.88)

        self.box_actor = vtk.vtkActor()
        self.info_panel_actor = make_panel_actor(0.015, 0.58, 0.29, 0.98, PANEL_COLOR)
        self.controls_panel_actor = make_panel_actor(0.015, 0.02, 0.31, 0.28, PANEL_COLOR)
        self.title_actor = make_text_actor(0.035, 0.95, font_size=17)
        self.progress_actor = make_text_actor(0.035, 0.865, font_size=18)
        self.coords_actor = make_text_actor(0.035, 0.805, font_size=17, bold=False)
        self.status_actor = make_text_actor(0.035, 0.735, font_size=11, bold=False)
        self.controls_actor = make_text_actor(0.035, 0.255, font_size=13, bold=False)

        self.renderer.AddActor2D(self.info_panel_actor)
        self.renderer.AddActor2D(self.controls_panel_actor)
        self.renderer.AddActor(self.box_actor)
        self.renderer.AddActor(self.base_actor)
        self.renderer.AddActor(self.visited_actor)
        self.renderer.AddActor(self.current_actor)
        self.renderer.AddActor(self.path_actor)
        self.renderer.AddActor2D(self.title_actor)
        self.renderer.AddActor2D(self.progress_actor)
        self.renderer.AddActor2D(self.coords_actor)
        self.renderer.AddActor2D(self.status_actor)
        self.renderer.AddActor2D(self.controls_actor)

        self.m_slider = make_slider(
            self.interactor,
            3,
            self.max_m,
            self.m,
            "m",
            (0.05, 0.19),
            (0.28, 0.19),
        )
        self.speed_slider = make_slider(
            self.interactor,
            20,
            500,
            self.speed_ms,
            "speed ms",
            (0.05, 0.12),
            (0.28, 0.12),
        )

        self.m_slider.AddObserver("EndInteractionEvent", self.on_m_slider)
        self.speed_slider.AddObserver("EndInteractionEvent", self.on_speed_slider)
        self.interactor.AddObserver("TimerEvent", self.on_timer)
        self.interactor.AddObserver("KeyPressEvent", self.on_key_press)

        button_style = {
            "background_color": (34, 39, 49),
            "text_color": TEXT_COLOR,
            "font_size": 15,
        }
        self.play_button = Button2D(self.interactor, "Play", (42, 28), (86, 36), **button_style)
        self.pause_button = Button2D(self.interactor, "Pause", (138, 28), (86, 36), **button_style)
        self.reset_button = Button2D(self.interactor, "Reset", (234, 28), (86, 36), **button_style)

        self.play_button.add_click_callback(self.play)
        self.pause_button.add_click_callback(self.pause_and_refresh)
        self.reset_button.add_click_callback(self.reset_animation)

        self.rebuild_scene(reset_camera=True)

    def build_points_polydata(self) -> vtk.vtkPolyData:
        points = vtk.vtkPoints()

        self.vertex_positions = {}
        for i in range(self.m):
            for j in range(self.m):
                for k in range(self.m):
                    vertex = (i, j, k)
                    position = lattice_position(i, j, k, self.m)
                    points.InsertNextPoint(*position)
                    self.vertex_positions[vertex] = position

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        return polydata

    def rebuild_scene(self, reset_camera: bool = False) -> None:
        self.pause()
        self.order = build_cycle_order(self.m)
        self.completed_steps = 0
        self.complete = False

        new_base = self.build_points_polydata()
        self.base_polydata.DeepCopy(new_base)
        self.base_polydata.Modified()

        self.update_box()
        self.update_dynamic_geometry()
        self.update_text()

        if reset_camera:
            self.reset_camera()

        self.render_window.Render()

    def update_box(self) -> None:
        cube = vtk.vtkCubeSource()
        side = (self.m - 1) + 0.9
        cube.SetXLength(side)
        cube.SetYLength(side)
        cube.SetZLength(side)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cube.GetOutputPort())
        self.box_actor.SetMapper(mapper)
        self.box_actor.GetProperty().SetRepresentationToWireframe()
        self.box_actor.GetProperty().SetColor(*BOX_COLOR)
        self.box_actor.GetProperty().SetOpacity(0.18)
        self.box_actor.GetProperty().SetLineWidth(1.0)

    def update_dynamic_geometry(self) -> None:
        self.update_visited_points()
        self.update_current_point()
        self.update_path()

    def update_visited_points(self) -> None:
        points = vtk.vtkPoints()
        for index in range(self.completed_steps):
            points.InsertNextPoint(*self.vertex_positions[self.order[index]])
        self.visited_polydata.SetPoints(points)
        self.visited_polydata.Modified()

    def update_current_point(self) -> None:
        points = vtk.vtkPoints()
        if self.complete:
            current_vertex = self.order[0]
        else:
            current_vertex = self.order[self.completed_steps]
        points.InsertNextPoint(*self.vertex_positions[current_vertex])
        self.current_polydata.SetPoints(points)
        self.current_polydata.Modified()

    def update_path(self) -> None:
        path_points = vtk.vtkPoints()
        path_vertices = self.order[: self.completed_steps + 1]
        if self.complete:
            path_vertices = path_vertices + [self.order[0]]

        for vertex in path_vertices:
            path_points.InsertNextPoint(*self.vertex_positions[vertex])

        lines = vtk.vtkCellArray()
        if path_points.GetNumberOfPoints() >= 2:
            polyline = vtk.vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(path_points.GetNumberOfPoints())
            for point_id in range(path_points.GetNumberOfPoints()):
                polyline.GetPointIds().SetId(point_id, point_id)
            lines.InsertNextCell(polyline)

        self.path_polydata.SetPoints(path_points)
        self.path_polydata.SetLines(lines)
        self.path_polydata.Modified()

    def current_vertex(self) -> tuple[int, int, int]:
        if self.complete:
            return self.order[0]
        return self.order[self.completed_steps]

    def current_index(self) -> int:
        if self.complete:
            return 0
        return self.completed_steps

    def update_text(self) -> None:
        current = self.current_vertex()
        total = len(self.order)
        current_index = self.current_index()
        state = "playing" if self.playing else "complete" if self.complete else "paused"
        progress = f"index {current_index}"
        coords = f"({current[0]}, {current[1]}, {current[2]})"
        status = (
            f"Vertices: {total}\n"
            f"Step: {self.completed_steps + 1} / {total}\n"
            f"m: {self.m}\n"
            f"State: {state}\n"
            f"Speed: {self.speed_ms} ms"
        )
        controls = "Space: play/pause\nR: reset"
        self.title_actor.SetInput("Knuth cycle on Z_m^3")
        self.progress_actor.SetInput(progress)
        self.coords_actor.SetInput(coords)
        self.status_actor.SetInput(status)
        self.controls_actor.SetInput(controls)

    def reset_camera(self) -> None:
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.Azimuth(35)
        camera.Elevation(28)
        camera.ParallelProjectionOn()
        camera.SetParallelScale(max(3.5, (self.m - 1) * 1.2))
        self.renderer.ResetCameraClippingRange()

    def restart_timer(self) -> None:
        if self.timer_id is not None:
            self.interactor.DestroyTimer(self.timer_id)
        self.timer_id = self.interactor.CreateRepeatingTimer(self.speed_ms)

    def play(self) -> None:
        if self.complete:
            self.reset_animation()
        if not self.playing:
            self.playing = True
            self.restart_timer()
            self.update_text()
            self.render_window.Render()

    def pause(self) -> None:
        if self.timer_id is not None:
            self.interactor.DestroyTimer(self.timer_id)
            self.timer_id = None
        self.playing = False

    def pause_and_refresh(self) -> None:
        self.pause()
        self.update_text()
        self.render_window.Render()

    def reset_animation(self) -> None:
        self.pause()
        self.completed_steps = 0
        self.complete = False
        self.update_dynamic_geometry()
        self.update_text()
        self.render_window.Render()

    def advance(self) -> None:
        if self.complete:
            return

        if self.completed_steps < len(self.order) - 1:
            self.completed_steps += 1
        else:
            self.complete = True
            self.pause()

        self.update_dynamic_geometry()
        self.update_text()
        self.render_window.Render()

    def on_timer(self, _obj: vtk.vtkObject, _event: str) -> None:
        if self.playing:
            self.advance()

    def on_m_slider(self, obj: vtk.vtkObject, _event: str) -> None:
        slider_value = obj.GetRepresentation().GetValue()
        new_m = nearest_odd(slider_value, 3, self.max_m)
        obj.GetRepresentation().SetValue(new_m)
        if new_m != self.m:
            was_playing = self.playing
            self.m = new_m
            self.rebuild_scene(reset_camera=True)
            if was_playing:
                self.play()

    def on_speed_slider(self, obj: vtk.vtkObject, _event: str) -> None:
        new_speed = max(20, int(round(obj.GetRepresentation().GetValue())))
        obj.GetRepresentation().SetValue(new_speed)
        if new_speed != self.speed_ms:
            self.speed_ms = new_speed
            if self.playing:
                self.restart_timer()
            self.update_text()
            self.render_window.Render()

    def on_key_press(self, obj: vtk.vtkObject, _event: str) -> None:
        key = obj.GetKeySym().lower()
        if key == "space":
            if self.playing:
                self.pause()
            else:
                self.play()
            self.update_text()
            self.render_window.Render()
        elif key == "r":
            self.reset_animation()

    def run(self) -> None:
        self.render_window.Render()
        self.interactor.Initialize()
        self.interactor.Start()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=5, help="Odd lattice size m >= 3")
    parser.add_argument(
        "--speed-ms",
        type=int,
        default=120,
        help="Animation timer interval in milliseconds",
    )
    parser.add_argument(
        "--max-m",
        type=int,
        default=15,
        help="Largest odd m allowed by the slider",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    app = KnuthHamiltonianCycleApp(
        initial_m=validate_odd_m(args.m),
        initial_speed_ms=max(20, args.speed_ms),
        max_m=args.max_m,
    )
    app.run()


if __name__ == "__main__":
    main()
