"""
Ahmed Body Visualization: improved geometry processing and field visualization.

Features:
- Robust STL loading with cleaning, triangulation, optional smoothing, and normals.
- Multiple scalar fields (synthetic pressure, curvature, height, length).
- Diverging/linear colormaps with scalar bar.
- Feature edges and wireframe overlays for surface topology cues.
- Interactive key controls for field cycling and view toggles.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import vtk


FIELD_PRESSURE = "StaticPressure"
FIELD_CURVATURE = "GaussianCurvature"
FIELD_HEIGHT = "Height"
FIELD_X = "NormalizedX"

FIELD_ORDER = [FIELD_PRESSURE, FIELD_CURVATURE, FIELD_HEIGHT, FIELD_X]
FIELD_META = {
    FIELD_PRESSURE: {"label": "Static Pressure (synthetic)", "diverging": True},
    FIELD_CURVATURE: {"label": "Gaussian Curvature", "diverging": True},
    FIELD_HEIGHT: {"label": "Normalized Height", "diverging": False},
    FIELD_X: {"label": "Normalized Length", "diverging": False},
}


def find_repo_root(start: Path) -> Path | None:
    """Walk upward to find the repository root containing the Ahmed STL."""
    for parent in (start,) + tuple(start.parents):
        if (parent / "data" / "stls" / "ahmed_body.stl").exists():
            return parent
    return None


def default_stl_path() -> Path:
    repo_root = find_repo_root(Path(__file__).resolve().parent)
    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[2]
    return repo_root / "data" / "stls" / "ahmed_body.stl"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ahmed Body visualization demo.")
    parser.add_argument(
        "--stl",
        type=str,
        default=str(default_stl_path()),
        help="Path to the Ahmed body STL file.",
    )
    parser.add_argument(
        "--field",
        type=str,
        choices=FIELD_ORDER,
        default=FIELD_PRESSURE,
        help="Initial scalar field to display.",
    )
    parser.add_argument(
        "--smooth-iterations",
        type=int,
        default=12,
        help="Number of smoothing iterations (0 disables smoothing).",
    )
    parser.add_argument(
        "--feature-angle",
        type=float,
        default=40.0,
        help="Feature angle for edge extraction in degrees.",
    )
    parser.add_argument(
        "--no-edges",
        action="store_true",
        help="Disable feature edge overlay.",
    )
    parser.add_argument(
        "--wireframe",
        action="store_true",
        help="Enable wireframe overlay on startup.",
    )
    return parser.parse_args()


def resolve_stl_path(path_str: str) -> Path:
    path = Path(path_str)
    if path.exists():
        return path

    if not path.is_absolute():
        repo_root = find_repo_root(Path(__file__).resolve().parent)
        if repo_root is not None:
            candidate = repo_root / path
            if candidate.exists():
                return candidate

        script_root = Path(__file__).resolve().parents[2]
        candidate = script_root / path
        if candidate.exists():
            return candidate

    raise FileNotFoundError(f"STL file not found: {path_str}")


def load_and_preprocess_geometry(path: Path, smooth_iterations: int) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(path))

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(reader.GetOutputPort())

    triangulate = vtk.vtkTriangleFilter()
    triangulate.SetInputConnection(cleaner.GetOutputPort())

    normals_input = triangulate
    if smooth_iterations > 0:
        smoother = vtk.vtkWindowedSincPolyDataFilter()
        smoother.SetInputConnection(triangulate.GetOutputPort())
        smoother.SetNumberOfIterations(smooth_iterations)
        smoother.BoundarySmoothingOff()
        smoother.FeatureEdgeSmoothingOff()
        smoother.NonManifoldSmoothingOn()
        smoother.NormalizeCoordinatesOn()
        normals_input = smoother

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(normals_input.GetOutputPort())
    normals.SetComputePointNormals(True)
    normals.SetComputeCellNormals(False)
    normals.SetAutoOrientNormals(True)
    normals.SetConsistency(True)
    normals.SetFeatureAngle(60.0)
    normals.SplittingOn()
    normals.Update()

    output = normals.GetOutput()
    if output.GetNumberOfPoints() == 0:
        raise ValueError("Loaded STL contains no points.")
    return output


def add_array(
    polydata: vtk.vtkPolyData, name: str, values: list[float]
) -> vtk.vtkDoubleArray:
    array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples(len(values))
    for idx, value in enumerate(values):
        array.SetTuple1(idx, value)
    polydata.GetPointData().AddArray(array)
    return array


def compute_percentile_range(
    values: list[float], lower: float = 2.0, upper: float = 98.0
) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    if len(values) < 10:
        return min(values), max(values)
    values.sort()
    low_idx = max(int(len(values) * lower / 100.0), 0)
    high_idx = min(int(len(values) * upper / 100.0), len(values) - 1)
    return values[low_idx], values[high_idx]


def compute_robust_range(
    array: vtk.vtkDoubleArray, sample_limit: int = 50000
) -> tuple[float, float]:
    num_values = array.GetNumberOfTuples()
    if num_values == 0:
        return 0.0, 0.0
    step = max(num_values // sample_limit, 1)
    sampled = [array.GetTuple1(i) for i in range(0, num_values, step)]
    return compute_percentile_range(sampled, 2.0, 98.0)


def add_scalar_fields(polydata: vtk.vtkPolyData) -> dict[str, tuple[float, float]]:
    bounds = polydata.GetBounds()
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    length = max(x_max - x_min, 1e-6)
    width = max(y_max - y_min, 1e-6)
    height = max(z_max - z_min, 1e-6)

    pressure_values = []
    height_values = []
    x_values = []

    num_points = polydata.GetNumberOfPoints()
    for i in range(num_points):
        x, y, z = polydata.GetPoint(i)
        nx = (x - x_min) / length
        ny = (y - y_min) / width
        nz = (z - z_min) / height

        front_peak = math.exp(-(((nx - 0.08) / 0.08) ** 2))
        wake_suction = math.exp(-(((nx - 0.9) / 0.12) ** 2))
        side_bias = math.cos(math.pi * (ny - 0.5))
        roof_bias = 1.0 - abs(nz - 0.55) * 1.4

        pressure = (
            240.0 * front_peak
            - 180.0 * wake_suction
            + 45.0 * side_bias
            + 35.0 * roof_bias
        )
        pressure_values.append(pressure)
        height_values.append(nz)
        x_values.append(nx)

    pressure_array = add_array(polydata, FIELD_PRESSURE, pressure_values)
    height_array = add_array(polydata, FIELD_HEIGHT, height_values)
    x_array = add_array(polydata, FIELD_X, x_values)

    curvature_filter = vtk.vtkCurvatures()
    curvature_filter.SetInputData(polydata)
    curvature_filter.SetCurvatureTypeToGaussian()
    curvature_filter.Update()
    curvature_array_source = curvature_filter.GetOutput().GetPointData().GetScalars()
    curvature_array = vtk.vtkDoubleArray()
    curvature_array.DeepCopy(curvature_array_source)
    curvature_array.SetName(FIELD_CURVATURE)
    polydata.GetPointData().AddArray(curvature_array)

    ranges = {
        FIELD_PRESSURE: pressure_array.GetRange(),
        FIELD_HEIGHT: height_array.GetRange(),
        FIELD_X: x_array.GetRange(),
        FIELD_CURVATURE: compute_robust_range(curvature_array),
    }

    return ranges


def make_diverging_lut(
    value_range: tuple[float, float]
) -> vtk.vtkColorTransferFunction:
    min_val, max_val = value_range
    max_abs = max(abs(min_val), abs(max_val))
    color_tf = vtk.vtkColorTransferFunction()
    color_tf.SetColorSpaceToDiverging()
    color_tf.AddRGBPoint(-max_abs, 0.22, 0.36, 0.76)
    color_tf.AddRGBPoint(0.0, 0.95, 0.95, 0.95)
    color_tf.AddRGBPoint(max_abs, 0.78, 0.18, 0.18)
    return color_tf


def make_sequential_lut() -> vtk.vtkLookupTable:
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetHueRange(0.62, 0.08)
    lut.SetSaturationRange(0.7, 0.9)
    lut.SetValueRange(0.85, 0.95)
    lut.Build()
    return lut


def configure_light_kit(light_kit: vtk.vtkLightKit) -> None:
    """Configure a light kit with compatibility for older VTK versions."""
    light_kit.SetKeyLightIntensity(0.9)

    if hasattr(light_kit, "SetFillLightIntensity"):
        light_kit.SetFillLightIntensity(0.4)
    elif hasattr(light_kit, "SetKeyToFillRatio"):
        light_kit.SetKeyToFillRatio(2.5)

    if hasattr(light_kit, "SetBackLightIntensity"):
        light_kit.SetBackLightIntensity(0.3)
    elif hasattr(light_kit, "SetKeyToBackRatio"):
        light_kit.SetKeyToBackRatio(3.5)

    if hasattr(light_kit, "SetHeadLightIntensity"):
        light_kit.SetHeadLightIntensity(0.6)
    elif hasattr(light_kit, "SetKeyToHeadRatio"):
        light_kit.SetKeyToHeadRatio(1.6)


class AhmedBodyViewer:
    def __init__(
        self,
        polydata: vtk.vtkPolyData,
        field_ranges: dict[str, tuple[float, float]],
        initial_field: str,
        feature_angle: float,
        show_edges: bool,
        show_wireframe: bool,
    ) -> None:
        self.polydata = polydata
        self.field_ranges = field_ranges
        self.field_order = FIELD_ORDER
        self.active_field = initial_field

        self.luts = {}
        for field_name in self.field_order:
            field_range = field_ranges.get(field_name, (0.0, 1.0))
            if FIELD_META[field_name]["diverging"]:
                self.luts[field_name] = make_diverging_lut(field_range)
            else:
                self.luts[field_name] = make_sequential_lut()

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputData(polydata)
        self.mapper.SetScalarModeToUsePointData()
        self.mapper.SetColorModeToMapScalars()
        self.mapper.ScalarVisibilityOn()

        self.surface_actor = vtk.vtkActor()
        self.surface_actor.SetMapper(self.mapper)
        self.surface_actor.GetProperty().SetDiffuse(0.85)
        self.surface_actor.GetProperty().SetSpecular(0.3)
        self.surface_actor.GetProperty().SetSpecularPower(25.0)

        self.wireframe_actor = vtk.vtkActor()
        wireframe_mapper = vtk.vtkPolyDataMapper()
        wireframe_mapper.SetInputData(polydata)
        wireframe_mapper.ScalarVisibilityOff()
        self.wireframe_actor.SetMapper(wireframe_mapper)
        self.wireframe_actor.GetProperty().SetRepresentationToWireframe()
        self.wireframe_actor.GetProperty().SetColor(0.1, 0.1, 0.1)
        self.wireframe_actor.GetProperty().SetOpacity(0.25)
        self.wireframe_actor.SetVisibility(show_wireframe)

        edge_filter = vtk.vtkFeatureEdges()
        edge_filter.SetInputData(polydata)
        edge_filter.SetFeatureAngle(feature_angle)
        edge_filter.BoundaryEdgesOn()
        edge_filter.FeatureEdgesOn()
        edge_filter.ManifoldEdgesOff()
        edge_filter.NonManifoldEdgesOff()

        edge_mapper = vtk.vtkPolyDataMapper()
        edge_mapper.SetInputConnection(edge_filter.GetOutputPort())
        edge_mapper.ScalarVisibilityOff()

        self.edge_actor = vtk.vtkActor()
        self.edge_actor.SetMapper(edge_mapper)
        self.edge_actor.GetProperty().SetColor(0.05, 0.05, 0.05)
        self.edge_actor.GetProperty().SetLineWidth(1.8)
        self.edge_actor.SetVisibility(show_edges)

        self.scalar_bar = vtk.vtkScalarBarActor()
        self.scalar_bar.SetNumberOfLabels(5)
        self.scalar_bar.GetLabelTextProperty().SetColor(0.12, 0.12, 0.12)
        self.scalar_bar.GetTitleTextProperty().SetColor(0.12, 0.12, 0.12)

        self.corner_annotation = vtk.vtkCornerAnnotation()
        self.corner_annotation.GetTextProperty().SetColor(0.12, 0.12, 0.12)
        self.corner_annotation.SetLinearFontScaleFactor(2.0)
        self.corner_annotation.SetNonlinearFontScaleFactor(1.0)
        self.corner_annotation.SetMaximumFontSize(14)
        self.corner_annotation.SetText(2, "Ahmed Body Visualization")
        self.corner_annotation.SetText(
            0,
            "Keys: [K] field  [1-4] select  [G] edges  [V] wireframe  [B] bar  [R] reset  [H] help",
        )
        self.corner_annotation.SetText(
            1,
            f"Points: {polydata.GetNumberOfPoints()}\nCells: {polydata.GetNumberOfCells()}",
        )

        self.renderer = vtk.vtkRenderer()
        self.renderer.GradientBackgroundOn()
        self.renderer.SetBackground(0.95, 0.95, 0.98)
        self.renderer.SetBackground2(0.72, 0.78, 0.85)
        self.renderer.AddActor(self.surface_actor)
        self.renderer.AddActor(self.wireframe_actor)
        self.renderer.AddActor(self.edge_actor)
        self.renderer.AddViewProp(self.scalar_bar)
        self.renderer.AddViewProp(self.corner_annotation)

        light_kit = vtk.vtkLightKit()
        configure_light_kit(light_kit)
        light_kit.AddLightsToRenderer(self.renderer)

        self.render_window = vtk.vtkRenderWindow()
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetSize(1000, 700)
        self.render_window.SetWindowName("Ahmed Body Visualization")

        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)
        self.interactor_style = AhmedBodyInteractorStyle(self)
        if hasattr(self.interactor_style, "SetDefaultRenderer"):
            self.interactor_style.SetDefaultRenderer(self.renderer)
        self.interactor.SetInteractorStyle(self.interactor_style)

        self.axes = vtk.vtkAxesActor()
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(self.axes)
        self.axes_widget.SetInteractor(self.interactor)
        self.axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
        self.axes_widget.SetEnabled(1)
        self.axes_widget.InteractiveOff()

        self.set_active_field(initial_field)

    def set_active_field(self, field_name: str) -> None:
        self.active_field = field_name
        self.polydata.GetPointData().SetActiveScalars(field_name)

        field_range = self.field_ranges.get(field_name, (0.0, 1.0))
        if FIELD_META[field_name]["diverging"]:
            max_abs = max(abs(field_range[0]), abs(field_range[1]))
            field_range = (-max_abs, max_abs)

        lookup = self.luts[field_name]
        if isinstance(lookup, vtk.vtkLookupTable):
            lookup.SetRange(field_range)
            lookup.Build()
        self.mapper.SetLookupTable(lookup)
        self.mapper.SetScalarRange(field_range)
        self.scalar_bar.SetLookupTable(lookup)
        self.scalar_bar.SetTitle(FIELD_META[field_name]["label"])

        range_text = f"Range: {field_range[0]:.3g} to {field_range[1]:.3g}"
        self.corner_annotation.SetText(
            3, f"{FIELD_META[field_name]['label']}\n{range_text}"
        )

    def reset_camera(self) -> None:
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.Azimuth(25)
        camera.Elevation(18)
        camera.Dolly(1.15)
        self.renderer.ResetCameraClippingRange()

    def handle_key(self, key: str) -> bool:
        if key == "k":
            idx = self.field_order.index(self.active_field)
            new_field = self.field_order[(idx + 1) % len(self.field_order)]
            self.set_active_field(new_field)
        elif key in ("1", "2", "3", "4"):
            index = int(key) - 1
            if 0 <= index < len(self.field_order):
                self.set_active_field(self.field_order[index])
        elif key == "g":
            self.edge_actor.SetVisibility(not self.edge_actor.GetVisibility())
        elif key == "v":
            self.wireframe_actor.SetVisibility(not self.wireframe_actor.GetVisibility())
        elif key == "b":
            self.scalar_bar.SetVisibility(not self.scalar_bar.GetVisibility())
        elif key == "r":
            self.reset_camera()
        elif key == "h":
            print(
                "Keys: k=cycle field, 1-4=select field, g=edges, v=wireframe, "
                "b=scalar bar, r=reset camera"
            )
        else:
            return False

        self.render_window.Render()
        return True

    def start(self) -> None:
        self.reset_camera()
        self.render_window.Render()
        self.interactor.Initialize()
        self.interactor.Start()


class AhmedBodyInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
    def __init__(self, viewer: AhmedBodyViewer):
        super().__init__()
        self.viewer = viewer
        self.AddObserver("KeyPressEvent", self._on_key_press)

    def _normalized_key(self, interactor: vtk.vtkRenderWindowInteractor) -> str:
        key_sym = interactor.GetKeySym()
        if key_sym:
            return key_sym.lower()
        key_code = interactor.GetKeyCode()
        if isinstance(key_code, int):
            if key_code == 0:
                return ""
            key_code = chr(key_code)
        if not key_code:
            return ""
        key_str = str(key_code)
        if key_str == " ":
            return "space"
        return key_str.lower()

    def _on_key_press(self, obj: vtk.vtkObject, event: str) -> None:
        interactor = self.GetInteractor()
        if interactor is None:
            return
        key = self._normalized_key(interactor)
        if key and self.viewer.handle_key(key):
            interactor.AbortFlagOn()


def main() -> None:
    args = parse_args()
    stl_path = resolve_stl_path(args.stl)
    polydata = load_and_preprocess_geometry(stl_path, args.smooth_iterations)
    field_ranges = add_scalar_fields(polydata)
    polydata.GetPointData().SetActiveScalars(args.field)

    viewer = AhmedBodyViewer(
        polydata=polydata,
        field_ranges=field_ranges,
        initial_field=args.field,
        feature_angle=args.feature_angle,
        show_edges=not args.no_edges,
        show_wireframe=args.wireframe,
    )
    viewer.start()


if __name__ == "__main__":
    main()
