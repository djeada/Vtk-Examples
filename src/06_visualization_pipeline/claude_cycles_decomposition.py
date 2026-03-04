"""
Interactive explorer for Claude's cycle decomposition from Knuth's 2026 note.

Paper:
https://www-cs-faculty.stanford.edu/~knuth/papers/claude-cycles.pdf

This visualizes the directed base graph G_m on Z_m^3:
  vertices: (i, j, k) with 0 <= i, j, k < m
  arcs from each vertex: +e0, +e1, +e2 (mod m)

Claude's local permutation rule assigns one outgoing direction per cycle index
c in {0,1,2}. For odd m > 1, each chosen edge set is Hamiltonian.

Interactive controls:
  - Slider "m (odd)": change modulus m
  - Discrete cycle choice: keys 1/2/3 (or Left/Right arrows)
  - Toggle Graph/Table mode: T
  - In Graph mode: click a vertex to inspect outgoing arcs and chosen cycle edge
  - In Table mode: browse all edges and which cycle uses each edge
    (PageUp/PageDown or [ / ])
  - Base graph is always visible in Graph mode
  - P: toggle vertex points (Graph mode)
  - H: toggle explanatory overlays
  - R: reset camera
  - +/-: decrease/increase m by 2

Examples:
  python src/06_visualization_pipeline/claude_cycles_decomposition.py --m 5
  python src/06_visualization_pipeline/claude_cycles_decomposition.py --m 7 --cycle 2
  python src/06_visualization_pipeline/claude_cycles_decomposition.py --offscreen --output /tmp/claude_m5.png
"""

import argparse
import math

import vtk


CYCLE_COLORS = [
    (0.92, 0.28, 0.20),  # cycle 0
    (0.15, 0.58, 0.92),  # cycle 1
    (0.18, 0.74, 0.35),  # cycle 2
]

DIRECTION_COLORS = {
    "0": (0.90, 0.30, 0.28),  # +e0 (i)
    "1": (0.26, 0.61, 0.96),  # +e1 (j)
    "2": (0.25, 0.80, 0.42),  # +e2 (k)
}
SELECTED_EDGE_COLOR = (1.00, 0.85, 0.20)

BACKGROUND_COLOR = (0.07, 0.08, 0.11)
POINT_COLOR = (0.90, 0.90, 0.92)
BASE_GRAPH_COLOR = (0.74, 0.76, 0.80)
TEXT_COLOR = (0.97, 0.97, 0.98)
TABLE_ROWS_PER_PAGE = 22


def claude_permutation(i, j, k, m):
    """
    Return local permutation d in {'0','1','2'}^3 for vertex (i,j,k).

    Rule from Knuth's C-style code:
      s = (i + j + k) mod m
      if s == 0:      d = "012" if j == m - 1 else "210"
      elif s == m-1:  d = "120" if i > 0 else "210"
      else:           d = "201" if i == m - 1 else "102"
    """
    s = (i + j + k) % m
    if s == 0:
        return "012" if j == m - 1 else "210"
    if s == m - 1:
        return "120" if i > 0 else "210"
    return "201" if i == m - 1 else "102"


def apply_direction(vertex, direction, m):
    """Take one directed step in G_m according to direction in {'0','1','2'}."""
    i, j, k = vertex
    if direction == "0":
        return ((i + 1) % m, j, k)
    if direction == "1":
        return (i, (j + 1) % m, k)
    return (i, j, (k + 1) % m)


def build_successors(m):
    """Build successor map for each cycle choice c=0,1,2."""
    successors = [dict(), dict(), dict()]
    for i in range(m):
        for j in range(m):
            for k in range(m):
                perm = claude_permutation(i, j, k, m)
                source = (i, j, k)
                for cycle_index in range(3):
                    direction = perm[cycle_index]
                    successors[cycle_index][source] = apply_direction(source, direction, m)
    return successors


def cycle_stats(successor_map, m):
    """
    Diagnostics for one directed edge selection.

    Hamiltonian means:
    1) indegree is exactly 1 at every vertex, and
    2) one traversal visits all m^3 vertices then returns to start.
    """
    total_vertices = m**3
    indegrees = {vertex: 0 for vertex in successor_map}
    for target in successor_map.values():
        indegrees[target] += 1

    permutation = all(degree == 1 for degree in indegrees.values())
    min_indegree = min(indegrees.values())
    max_indegree = max(indegrees.values())

    start = next(iter(successor_map))
    visited = set()
    current = start
    while current not in visited:
        visited.add(current)
        current = successor_map[current]

    full_cycle = current == start and len(visited) == total_vertices
    return {
        "hamiltonian": permutation and full_cycle,
        "cycle_length": len(visited),
        "expected_length": total_vertices,
        "min_indegree": min_indegree,
        "max_indegree": max_indegree,
    }


def embed_vertex(i, j, k, m):
    """
    Embed Z_m^3 in R^3 for visual continuity under modulo wrap.

    This is a smooth trigonometric embedding, not graph distance.
    """
    tau = 2.0 * math.pi
    a = tau * i / m
    b = tau * j / m
    c = tau * k / m

    x = math.cos(a) + 0.55 * math.cos(b) + 0.30 * math.cos(c)
    y = math.sin(a) + 0.55 * math.sin(b) + 0.30 * math.sin(c)
    z = 0.75 * math.sin(a - b) + 0.45 * math.sin(b - c) + 0.30 * math.sin(c - a)
    return x, y, z


def build_points_and_ids(m):
    """Return vtkPoints, vertex->point_id map, and point_id->vertex list."""
    points = vtk.vtkPoints()
    vertex_to_id = {}
    id_to_vertex = []
    for i in range(m):
        for j in range(m):
            for k in range(m):
                pid = points.InsertNextPoint(*embed_vertex(i, j, k, m))
                vertex = (i, j, k)
                vertex_to_id[vertex] = pid
                id_to_vertex.append(vertex)
    return points, vertex_to_id, id_to_vertex


def build_cycle_polydata(points, vertex_to_id, successor_map):
    """PolyData with one directed edge per vertex for selected cycle."""
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    lines = vtk.vtkCellArray()
    for source, target in successor_map.items():
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, vertex_to_id[source])
        line.GetPointIds().SetId(1, vertex_to_id[target])
        lines.InsertNextCell(line)
    polydata.SetLines(lines)
    return polydata


def build_base_graph_polydata(points, vertex_to_id, m):
    """PolyData of all base graph arcs (+e0,+e1,+e2 at each vertex)."""
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    lines = vtk.vtkCellArray()
    for i in range(m):
        for j in range(m):
            for k in range(m):
                source = (i, j, k)
                targets = (
                    ((i + 1) % m, j, k),
                    (i, (j + 1) % m, k),
                    (i, j, (k + 1) % m),
                )
                for target in targets:
                    line = vtk.vtkLine()
                    line.GetPointIds().SetId(0, vertex_to_id[source])
                    line.GetPointIds().SetId(1, vertex_to_id[target])
                    lines.InsertNextCell(line)
    polydata.SetLines(lines)
    return polydata


def make_tube_actor(polydata, color, radius=0.02, opacity=0.95):
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(polydata)
    tube.SetRadius(radius)
    tube.SetNumberOfSides(14)
    tube.CappingOn()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    return actor


def make_segment_actor(p0, p1, color, radius=0.03, opacity=0.95):
    """Create a tube-rendered segment actor between two 3D points."""
    line_source = vtk.vtkLineSource()
    line_source.SetPoint1(*p0)
    line_source.SetPoint2(*p1)
    line_source.Update()

    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(line_source.GetOutputPort())
    tube.SetRadius(radius)
    tube.SetNumberOfSides(16)
    tube.CappingOn()

    tube_mapper = vtk.vtkPolyDataMapper()
    tube_mapper.SetInputConnection(tube.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(tube_mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    return actor


def make_vertex_marker_actor(center, radius=0.06):
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(*center)
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(24)
    sphere.SetPhiResolution(24)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*SELECTED_EDGE_COLOR)
    actor.GetProperty().SetOpacity(1.0)
    return actor


def make_points_actor(points):
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    verts = vtk.vtkCellArray()
    for pid in range(points.GetNumberOfPoints()):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(pid)
    polydata.SetVerts(verts)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*POINT_COLOR)
    actor.GetProperty().SetPointSize(7.0)
    actor.GetProperty().SetOpacity(0.70)
    return actor


def make_text_actor(text, x, y, size=18):
    actor = vtk.vtkTextActor()
    actor.SetInput(text)
    prop = actor.GetTextProperty()
    prop.SetFontSize(size)
    prop.SetColor(*TEXT_COLOR)
    prop.SetBold(True)
    actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    actor.SetPosition(x, y)
    return actor


def make_slider(interactor, minimum, maximum, value, title, p1, p2):
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
    rep.SetTubeWidth(0.004)
    rep.SetLabelFormat("%.0f")

    title_prop = rep.GetTitleProperty()
    title_prop.SetColor(*TEXT_COLOR)
    title_prop.SetFontSize(16)
    label_prop = rep.GetLabelProperty()
    label_prop.SetColor(*TEXT_COLOR)
    label_prop.SetFontSize(14)
    tube_prop = rep.GetTubeProperty()
    tube_prop.SetColor(0.75, 0.75, 0.78)
    slider_prop = rep.GetSliderProperty()
    slider_prop.SetColor(0.95, 0.95, 0.96)

    widget = vtk.vtkSliderWidget()
    widget.SetInteractor(interactor)
    widget.SetRepresentation(rep)
    widget.SetAnimationModeToAnimate()
    widget.EnabledOn()
    return widget


def format_vertex(vertex):
    return f"({vertex[0]:>2},{vertex[1]:>2},{vertex[2]:>2})"


class ClaudeCyclesExplorer:
    """Interactive VTK app with graph and table modes."""

    def __init__(self, initial_m=5, max_m=13, initial_cycle=0):
        self.max_m = max(3, max_m if max_m % 2 == 1 else max_m - 1)
        self.m = self._normalize_m(initial_m)
        self.cycle_index = max(0, min(2, initial_cycle))

        self.view_mode = "graph"
        self.table_offset = 0
        self.edge_rows = []
        self.cycle_edge_counts = [0, 0, 0]

        self.show_points = True
        self.show_help = True
        self.selected_vertex = None

        self.successors = None
        self.cycle_diagnostics = None
        self.points = None
        self.vertex_to_id = None
        self.id_to_vertex = None

        self.base_actor = None
        self.points_actor = None
        self.cycle_actor = None
        self.inspect_actors = []

        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(*BACKGROUND_COLOR)

        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetSize(1600, 920)

        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)
        self.interactor.AddObserver("KeyPressEvent", self._on_key_press)
        self.interactor.AddObserver("LeftButtonPressEvent", self._on_left_click)
        self.interactor.AddObserver("MouseWheelForwardEvent", self._on_mouse_wheel_forward)
        self.interactor.AddObserver("MouseWheelBackwardEvent", self._on_mouse_wheel_backward)

        self.camera = vtk.vtkCamera()
        self.camera.SetPosition(0.0, -6.0, 2.8)
        self.camera.SetFocalPoint(0.0, 0.0, 0.0)
        self.camera.SetViewUp(0.0, 0.0, 1.0)
        self.renderer.SetActiveCamera(self.camera)

        self.title_actor = make_text_actor("", 0.02, 0.94, size=24)
        self.status_actor = make_text_actor("", 0.02, 0.89, size=18)
        self.def_actor = make_text_actor("", 0.02, 0.84, size=16)

        self.cycle_panel_title_actor = make_text_actor("", 0.73, 0.82, size=18)
        self.cycle_button_actors = [
            make_text_actor("", 0.73, 0.77, size=17),
            make_text_actor("", 0.73, 0.73, size=17),
            make_text_actor("", 0.73, 0.69, size=17),
        ]
        self.cycle_button_regions = [
            (0.71, 0.99, 0.76, 0.80),
            (0.71, 0.99, 0.72, 0.76),
            (0.71, 0.99, 0.68, 0.72),
        ]

        self.rule_actor = make_text_actor("", 0.02, 0.17, size=15)
        self.help_actor = make_text_actor("", 0.02, 0.03, size=15)
        self.inspect_actor = make_text_actor("", 0.69, 0.03, size=16)

        self.table_actor = make_text_actor("", 0.18, 0.78, size=15)
        self.table_actor.GetTextProperty().SetFontFamilyToCourier()
        self.table_actor.GetTextProperty().SetJustificationToLeft()
        self.table_actor.GetTextProperty().SetVerticalJustificationToTop()
        self.table_footer_actor = make_text_actor("", 0.18, 0.05, size=15)

        self.renderer.AddActor2D(self.title_actor)
        self.renderer.AddActor2D(self.status_actor)
        self.renderer.AddActor2D(self.def_actor)
        self.renderer.AddActor2D(self.cycle_panel_title_actor)
        for actor in self.cycle_button_actors:
            self.renderer.AddActor2D(actor)
        self.renderer.AddActor2D(self.rule_actor)
        self.renderer.AddActor2D(self.help_actor)
        self.renderer.AddActor2D(self.inspect_actor)
        self.renderer.AddActor2D(self.table_actor)
        self.renderer.AddActor2D(self.table_footer_actor)

        self.m_slider = make_slider(
            self.interactor,
            3,
            self.max_m,
            self.m,
            "m (odd)",
            (0.70, 0.61),
            (0.96, 0.61),
        )
        self.m_slider.AddObserver("InteractionEvent", self._on_m_slider)

        self.point_picker = vtk.vtkPointPicker()
        self.point_picker.SetTolerance(0.03)
        self.point_picker.PickFromListOn()

        self.rebuild_scene(reset_camera=True)

    def _normalize_m(self, raw_value):
        m = int(round(raw_value))
        m = max(3, min(self.max_m, m))
        if m % 2 == 0:
            if m == self.max_m:
                m -= 1
            else:
                m += 1
        return max(3, m)

    def _set_m(self, new_m, reset_camera=False):
        normalized = self._normalize_m(new_m)
        if normalized == self.m:
            return
        self.m = normalized
        self.table_offset = 0
        self.m_slider.GetRepresentation().SetValue(self.m)
        self.selected_vertex = None
        self.rebuild_scene(reset_camera=reset_camera)

    def _set_cycle(self, cycle):
        cycle = max(0, min(2, int(cycle)))
        if cycle == self.cycle_index:
            return
        self.cycle_index = cycle
        self._refresh_cycle_actor()
        self._refresh_inspection_actors()
        self._update_text()
        self._apply_view_mode_visibility()
        self.render_window.Render()

    def _toggle_view_mode(self):
        self.view_mode = "table" if self.view_mode == "graph" else "graph"
        if self.view_mode == "graph":
            # Defensive guard: if actors were dropped by backend/window events, rebuild them.
            if self.base_actor is None or self.cycle_actor is None or self.points_actor is None:
                self.rebuild_scene(reset_camera=False)
                return
        self._apply_view_mode_visibility()
        self._update_text()
        self.renderer.ResetCameraClippingRange()
        self.render_window.Render()

    def _max_table_offset(self):
        return max(0, len(self.edge_rows) - TABLE_ROWS_PER_PAGE)

    def _scroll_table(self, delta_rows):
        self.table_offset = max(
            0,
            min(self.table_offset + delta_rows, self._max_table_offset()),
        )
        self._update_text()
        self.render_window.Render()

    def _on_m_slider(self, obj, _event):
        proposed = obj.GetRepresentation().GetValue()
        normalized = self._normalize_m(proposed)
        obj.GetRepresentation().SetValue(normalized)
        if normalized != self.m:
            self.m = normalized
            self.table_offset = 0
            self.selected_vertex = None
            self.rebuild_scene(reset_camera=False)

    def _on_key_press(self, _obj, _event):
        key = self.interactor.GetKeySym().lower()

        if key in ("t", "v"):
            self._toggle_view_mode()
            return
        if key in ("1", "2", "3"):
            self._set_cycle(int(key) - 1)
            return
        if key == "left":
            self._set_cycle(self.cycle_index - 1)
            return
        if key == "right":
            self._set_cycle(self.cycle_index + 1)
            return
        if key in ("plus", "equal", "kp_add"):
            self._set_m(self.m + 2, reset_camera=False)
            return
        if key in ("minus", "underscore", "kp_subtract"):
            self._set_m(self.m - 2, reset_camera=False)
            return
        if key in ("bracketright", "pagedown"):
            if self.view_mode == "table":
                self._scroll_table(TABLE_ROWS_PER_PAGE)
            return
        if key in ("bracketleft", "pageup"):
            if self.view_mode == "table":
                self._scroll_table(-TABLE_ROWS_PER_PAGE)
            return
        if key == "down":
            if self.view_mode == "table":
                self._scroll_table(1)
            return
        if key == "up":
            if self.view_mode == "table":
                self._scroll_table(-1)
            return

        if key == "p":
            self.show_points = not self.show_points
            self._apply_view_mode_visibility()
        elif key == "h":
            self.show_help = not self.show_help
            self._apply_view_mode_visibility()
        elif key == "r":
            self.renderer.ResetCamera()
            self.renderer.ResetCameraClippingRange()
        else:
            return

        self._update_text()
        self.render_window.Render()

    def _on_mouse_wheel_forward(self, _obj, _event):
        if self.view_mode == "table":
            self._scroll_table(-3)

    def _on_mouse_wheel_backward(self, _obj, _event):
        if self.view_mode == "table":
            self._scroll_table(3)

    def _on_left_click(self, _obj, _event):
        x, y = self.interactor.GetEventPosition()
        width, height = self.render_window.GetSize()
        if width <= 0 or height <= 0:
            return

        xn = x / float(width)
        yn = y / float(height)

        for idx, (x0, x1, y0, y1) in enumerate(self.cycle_button_regions):
            if x0 <= xn <= x1 and y0 <= yn <= y1:
                self._set_cycle(idx)
                return

        if self.view_mode != "graph" or self.points_actor is None:
            return

        self.point_picker.InitializePickList()
        self.point_picker.AddPickList(self.points_actor)
        picked = self.point_picker.Pick(x, y, 0, self.renderer)
        if not picked:
            return

        point_id = self.point_picker.GetPointId()
        if point_id < 0 or point_id >= len(self.id_to_vertex):
            return

        self.selected_vertex = self.id_to_vertex[point_id]
        self._refresh_inspection_actors()
        self._update_text()
        self.render_window.Render()

    def _clear_graph_actors(self):
        if self.base_actor is not None:
            self.renderer.RemoveActor(self.base_actor)
            self.base_actor = None
        if self.points_actor is not None:
            self.renderer.RemoveActor(self.points_actor)
            self.points_actor = None
        if self.cycle_actor is not None:
            self.renderer.RemoveActor(self.cycle_actor)
            self.cycle_actor = None

    def _clear_inspection_actors(self):
        for actor in self.inspect_actors:
            self.renderer.RemoveActor(actor)
        self.inspect_actors = []

    def _refresh_cycle_actor(self):
        if self.cycle_actor is not None:
            self.renderer.RemoveActor(self.cycle_actor)
            self.cycle_actor = None

        cycle_polydata = build_cycle_polydata(
            self.points,
            self.vertex_to_id,
            self.successors[self.cycle_index],
        )
        self.cycle_actor = make_tube_actor(
            cycle_polydata,
            CYCLE_COLORS[self.cycle_index],
            radius=0.026,
            opacity=0.96,
        )
        self.renderer.AddActor(self.cycle_actor)

    def _refresh_inspection_actors(self):
        self._clear_inspection_actors()

        if self.selected_vertex is None:
            return

        source = self.selected_vertex
        source_pid = self.vertex_to_id[source]
        source_point = self.points.GetPoint(source_pid)

        outgoing = {
            "0": ((source[0] + 1) % self.m, source[1], source[2]),
            "1": (source[0], (source[1] + 1) % self.m, source[2]),
            "2": (source[0], source[1], (source[2] + 1) % self.m),
        }

        for direction, target_vertex in outgoing.items():
            target_point = self.points.GetPoint(self.vertex_to_id[target_vertex])
            direction_actor = make_segment_actor(
                source_point,
                target_point,
                DIRECTION_COLORS[direction],
                radius=0.018,
                opacity=0.88,
            )
            self.inspect_actors.append(direction_actor)
            self.renderer.AddActor(direction_actor)

        chosen_direction = claude_permutation(*source, self.m)[self.cycle_index]
        chosen_target = outgoing[chosen_direction]
        chosen_target_point = self.points.GetPoint(self.vertex_to_id[chosen_target])

        chosen_actor = make_segment_actor(
            source_point,
            chosen_target_point,
            SELECTED_EDGE_COLOR,
            radius=0.035,
            opacity=1.0,
        )
        marker_actor = make_vertex_marker_actor(source_point, radius=0.07)

        self.inspect_actors.append(chosen_actor)
        self.inspect_actors.append(marker_actor)
        self.renderer.AddActor(chosen_actor)
        self.renderer.AddActor(marker_actor)

    def _build_edge_rows(self):
        rows = []
        counts = [0, 0, 0]
        for i in range(self.m):
            for j in range(self.m):
                for k in range(self.m):
                    source = (i, j, k)
                    perm = claude_permutation(i, j, k, self.m)
                    direction_to_cycle = {perm[c]: c for c in range(3)}
                    for direction in ("0", "1", "2"):
                        target = apply_direction(source, direction, self.m)
                        cycle_id = direction_to_cycle[direction]
                        counts[cycle_id] += 1
                        rows.append((source, target, direction, cycle_id))
        self.edge_rows = rows
        self.cycle_edge_counts = counts
        self.table_offset = min(self.table_offset, self._max_table_offset())

    def _cycle_menu_entries(self):
        entries = []
        for idx in range(3):
            state = (
                "Hamiltonian"
                if self.cycle_diagnostics[idx]["hamiltonian"]
                else "Not Hamiltonian"
            )
            mark = ">>" if idx == self.cycle_index else "  "
            entries.append(f"{mark} [{idx + 1}] Cycle {idx}  {state}")
        return entries

    def _rule_text(self):
        return (
            "Local rule (d is a permutation of directions 0/1/2):\n"
            "s=(i+j+k) mod m\n"
            "if s=0:      d=012 if j=m-1 else 210\n"
            "elif s=m-1:  d=120 if i>0   else 210\n"
            "else:        d=201 if i=m-1 else 102\n"
            "Cycle c follows direction d[c] at every vertex."
        )

    def _inspect_text(self):
        if self.selected_vertex is None:
            return (
                "Vertex inspector:\n"
                "Click any vertex point to inspect\n"
                "+e0 (red), +e1 (blue), +e2 (green),\n"
                "selected cycle edge is highlighted in gold."
            )

        i, j, k = self.selected_vertex
        perm = claude_permutation(i, j, k, self.m)
        chosen = perm[self.cycle_index]
        axis = {"0": "i", "1": "j", "2": "k"}[chosen]
        return (
            f"Vertex inspector:\n"
            f"v=({i},{j},{k}), s={(i + j + k) % self.m}\n"
            f"d(v)={perm}; current cycle={self.cycle_index}\n"
            f"chosen direction={chosen} (increment {axis})"
        )

    def _update_table_text(self):
        total_edges = len(self.edge_rows)
        start = max(0, min(self.table_offset, self._max_table_offset()))
        end = min(start + TABLE_ROWS_PER_PAGE, total_edges)

        header = " idx FROM (i,j,k)     TO (i,j,k)       dir cycle"
        separator = "-" * len(header)
        lines = [header, separator]

        for idx in range(start, end):
            source, target, direction, cycle_id = self.edge_rows[idx]
            lines.append(
                f"{idx:>4} {format_vertex(source):<14} {format_vertex(target):<14} +e{direction}   C{cycle_id}"
            )

        if start == end:
            lines.append("(no edges)")

        self.table_actor.SetInput("\n".join(lines))
        self.table_footer_actor.SetInput(
            f"Rows {start}..{max(start, end - 1)} of {max(0, total_edges - 1)} | "
            "scroll: mouse wheel / Up-Down / [ ] / PageUp-PageDown"
        )

    def _apply_view_mode_visibility(self):
        graph_mode = self.view_mode == "graph"

        if self.base_actor is not None:
            self.base_actor.SetVisibility(1 if graph_mode else 0)
        if self.cycle_actor is not None:
            self.cycle_actor.SetVisibility(1 if graph_mode else 0)
        if self.points_actor is not None:
            self.points_actor.SetVisibility(1 if (graph_mode and self.show_points) else 0)
        for actor in self.inspect_actors:
            actor.SetVisibility(1 if graph_mode else 0)

        cycle_ui_vis = 1 if graph_mode else 0
        self.cycle_panel_title_actor.SetVisibility(cycle_ui_vis)
        for actor in self.cycle_button_actors:
            actor.SetVisibility(cycle_ui_vis)
        if graph_mode:
            self.m_slider.EnabledOn()
        else:
            self.m_slider.EnabledOff()

        if graph_mode:
            text_vis = 1 if self.show_help else 0
            self.rule_actor.SetVisibility(text_vis)
            self.help_actor.SetVisibility(text_vis)
            self.inspect_actor.SetVisibility(1)
            self.table_actor.SetVisibility(0)
            self.table_footer_actor.SetVisibility(0)
        else:
            self.rule_actor.SetVisibility(0)
            self.help_actor.SetVisibility(1 if self.show_help else 0)
            self.inspect_actor.SetVisibility(0)
            self.table_actor.SetVisibility(1)
            self.table_footer_actor.SetVisibility(1)

    def _update_text(self):
        total_vertices = self.m**3
        total_arcs = 3 * total_vertices
        selected = self.cycle_diagnostics[self.cycle_index]

        if self.view_mode == "graph":
            self.title_actor.SetInput(
                f"Claude Cycles Explorer | GRAPH | m={self.m} | viewing cycle {self.cycle_index}"
            )
            self.status_actor.SetInput(
                f"Selected cycle length={selected['cycle_length']}/{selected['expected_length']} | "
                f"indegree=[{selected['min_indegree']},{selected['max_indegree']}]"
            )
            self.def_actor.SetInput(
                f"Base graph G_m: {total_vertices} vertices, {total_arcs} directed arcs. "
                "Vertices v=(i,j,k), arcs: +e0,+e1,+e2 (mod m)."
            )

            self.cycle_panel_title_actor.SetInput("Cycle selection (click or keys 1/2/3):")
            entries = self._cycle_menu_entries()
            for idx, actor in enumerate(self.cycle_button_actors):
                actor.SetInput(entries[idx])
                if idx == self.cycle_index:
                    actor.GetTextProperty().SetColor(*SELECTED_EDGE_COLOR)
                else:
                    actor.GetTextProperty().SetColor(*TEXT_COLOR)

            self.rule_actor.SetInput(self._rule_text())
            points_state = "ON" if self.show_points else "OFF"
            self.help_actor.SetInput(
                "Mode: GRAPH | Controls: m slider (right), click cycle list or keys [1/2/3], Left/Right cycle, click vertex inspect, "
                "T table mode, P points, H help, R reset, +/- m by 2 "
                f"(points={points_state})"
            )
            self.inspect_actor.SetInput(self._inspect_text())
        else:
            self.title_actor.SetInput(
                f"Claude Cycles Explorer | TABLE | m={self.m} | complete edge partition"
            )
            self.status_actor.SetInput(
                f"All directed edges listed with assigned cycle C0/C1/C2 | "
                f"C0={self.cycle_edge_counts[0]}, C1={self.cycle_edge_counts[1]}, C2={self.cycle_edge_counts[2]}"
            )
            self.def_actor.SetInput(
                "Each row is one base-graph edge: FROM -> TO, direction (+e0/+e1/+e2), and cycle that uses it."
            )
            self.help_actor.SetInput(
                "Mode: TABLE | Controls: T graph mode, mouse wheel scroll, Up/Down line scroll, [ ] or PageUp/PageDown page scroll, H help"
            )
            self._update_table_text()

        self._apply_view_mode_visibility()

    def rebuild_scene(self, reset_camera=False):
        self._clear_graph_actors()
        self._clear_inspection_actors()

        self.successors = build_successors(self.m)
        self.cycle_diagnostics = [cycle_stats(self.successors[c], self.m) for c in range(3)]
        self.points, self.vertex_to_id, self.id_to_vertex = build_points_and_ids(self.m)
        self._build_edge_rows()

        base_polydata = build_base_graph_polydata(self.points, self.vertex_to_id, self.m)
        self.base_actor = make_tube_actor(
            base_polydata,
            BASE_GRAPH_COLOR,
            radius=0.007,
            opacity=0.16,
        )
        self.points_actor = make_points_actor(self.points)

        self.renderer.AddActor(self.base_actor)
        self.renderer.AddActor(self.points_actor)

        self._refresh_cycle_actor()
        self._refresh_inspection_actors()
        self._update_text()

        if reset_camera:
            self.renderer.ResetCamera()
        self.renderer.ResetCameraClippingRange()
        self.render_window.Render()

    def start(self):
        self.interactor.Initialize()
        self.render_window.Render()
        self.interactor.Start()


def save_screenshot(render_window, output_path):
    window_to_image = vtk.vtkWindowToImageFilter()
    window_to_image.SetInput(render_window)
    window_to_image.SetInputBufferTypeToRGB()
    window_to_image.ReadFrontBufferOff()
    window_to_image.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(output_path)
    writer.SetInputConnection(window_to_image.GetOutputPort())
    writer.Write()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Interactive visualization of Knuth's Claude cycles decomposition."
    )
    parser.add_argument(
        "--m",
        type=int,
        default=5,
        help="Initial odd modulus m (default: 5).",
    )
    parser.add_argument(
        "--max-m",
        type=int,
        default=13,
        help="Maximum odd m allowed by slider (default: 13).",
    )
    parser.add_argument(
        "--cycle",
        type=int,
        default=0,
        help="Initial cycle index (0,1,2).",
    )
    parser.add_argument(
        "--offscreen",
        action="store_true",
        help="Render once without interaction (use with --output).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="",
        help="Optional PNG output path.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    app = ClaudeCyclesExplorer(
        initial_m=args.m,
        max_m=args.max_m,
        initial_cycle=args.cycle,
    )

    if args.offscreen:
        if app.m != args.m:
            print(f"Adjusted m from {args.m} to odd value {app.m}.")
        if args.output:
            save_screenshot(app.render_window, args.output)
            print(f"Saved screenshot: {args.output}")
        else:
            print("Offscreen render completed (no --output path given).")
        return

    total_vertices = app.m**3
    total_arcs = 3 * total_vertices
    print(f"G_m with m={app.m}: vertices={total_vertices}, directed_arcs={total_arcs}")
    for cycle_index in range(3):
        stats = app.cycle_diagnostics[cycle_index]
        state = "Hamiltonian" if stats["hamiltonian"] else "Not Hamiltonian"
        print(
            f"Cycle {cycle_index}: {state}; "
            f"length={stats['cycle_length']}/{stats['expected_length']}; "
            f"indegree=[{stats['min_indegree']},{stats['max_indegree']}]"
        )
    print("Controls: m slider; [1/2/3] cycle; T graph/table; click vertex inspect in graph mode.")

    app.start()


if __name__ == "__main__":
    main()
