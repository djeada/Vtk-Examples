"""Scientific VTK scenes for reels, with explicit display transforms."""

from contextlib import redirect_stdout
import io
import math
from pathlib import Path
import runpy
from unittest.mock import patch

import numpy as np
from PIL import Image
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from src.common.field_visualization import field_points, quad_field
from src.reels.models import Model, example_module


BACKGROUND = (9, 15, 25)


def color_samples(name="inferno"):
    from matplotlib import colormaps

    # Trim only the darkest end for legibility; the scale remains monotonic.
    return (colormaps[name](np.linspace(0.13, 1, 256))[:, :3] * 255).astype(np.uint8)


def lookup_table(values, colors):
    low, high = float(np.nanmin(values)), float(np.nanmax(values))
    table = vtk.vtkLookupTable()
    table.SetNumberOfTableValues(256)
    table.SetRange(low, high if high > low else low + 1)
    table.SetNanColor(0.13, 0.19, 0.27, 1)
    for index, color in enumerate(colors):
        table.SetTableValue(index, *(color / 255), 1.0)
    table.Build()
    return table


def _actor(data, lookup=None, color=(1, 1, 1), opacity=1):
    mapper = vtk.vtkPolyDataMapper()
    if isinstance(data, vtk.vtkAlgorithmOutput):
        mapper.SetInputConnection(data)
    else:
        mapper.SetInputData(data)
    if lookup is None:
        mapper.ScalarVisibilityOff()
    else:
        mapper.SetLookupTable(lookup)
        mapper.UseLookupTableScalarRangeOn()
        mapper.InterpolateScalarsBeforeMappingOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().LightingOff()
    return actor


class SceneRenderer:
    """One reusable offscreen scene; captions are composed separately."""

    def __init__(self, model, width, height, generic_renderer=None):
        self.model = model
        # Recorded LBM snapshots are periodic over all NY rows, including the
        # interval from the last row back to the first.
        self.periodic_y = model.kind == "flow" and model.steps is not None
        self.y_period = (
            len(model.y) * (model.y[1] - model.y[0]) if self.periodic_y else None
        )
        renderers = generic_renderer or [vtk.vtkRenderer()]
        self.renderer = renderers[0]
        self.renderer.SetBackground(*(value / 255 for value in BACKGROUND))
        self.renderer.SetBackground2(0.025, 0.045, 0.065)
        self.renderer.GradientBackgroundOn()
        self.renderer.SetUseFXAA(True)
        self.window = vtk.vtkRenderWindow()
        self.window.SetOffScreenRendering(1)
        self.window.SetMultiSamples(0)
        self.window.SetSize(width, height)
        self.window.AddRenderer(self.renderer)
        for renderer in renderers[1:]:
            self.window.AddRenderer(renderer)
            renderer.GetActiveCamera().ParallelProjectionOn()
            renderer.ResetCamera()
        self.capture = vtk.vtkWindowToImageFilter()
        self.capture.SetInput(self.window)
        self.capture.SetInputBufferTypeToRGB()
        self.capture.ReadFrontBufferOff()
        self.capture.ShouldRerenderOff()
        self.colors = color_samples("viridis" if model.kind == "flow" else "inferno")
        self.table = lookup_table(model.fields, self.colors)
        self.mesh = None
        self.scalars = None
        self.phase = None
        self.legend = (self.table.GetRange(), model.unit)
        self.particles = None
        self._resources = []
        if generic_renderer is not None:
            # Existing example captions are replaced by the reel's own layout.
            # Keep scientific legends and every original viewport.
            for renderer in renderers:
                if len(renderers) == 1:
                    renderer.SetViewport(0, 0, 1, 1)
                for item in list(renderer.GetActors2D()):
                    if not item.IsA("vtkScalarBarActor"):
                        renderer.RemoveActor2D(item)
            self.legend = None
        elif model.kind == "rod":
            self._rod()
        else:
            self._field()
        self._frame_camera()

    def _rod(self):
        model = self.model
        module = example_module(model.key)
        self.mesh, _ = module.create_cylinder_rod(
            model.x, model.fields[0, 0], radius=0.065
        )
        rod = _actor(self.mesh, self.table)
        self.renderer.AddActor(rod)
        self._outline(self.mesh)
        for end in (0, -1):
            cap = module.create_end_cap(
                model.x[end], model.fields[0, 0, end], radius=0.065
            )
            self.renderer.AddActor(_actor(cap, self.table))
        # Arrow signs derive from the actual conductive gradient or prescribed u.
        if model.key == "heat_convection_solver_1d":
            arrows = np.ones_like(model.x)
        else:
            conductivity = (
                50 * (1 + 0.001 * model.fields[0, 0])
                if model.key == "enhanced_heat_transfer_solver_1d"
                else 50
            )
            arrows = -conductivity * np.gradient(model.fields[0, 0], model.x)
        points = vtk.vtkPoints()
        vectors = vtk.vtkFloatArray()
        vectors.SetNumberOfComponents(3)
        for index in np.linspace(5, len(model.x) - 6, 9).astype(int):
            points.InsertNextPoint(model.x[index], -0.16, 0.06)
            vectors.InsertNextTuple3(float(np.sign(arrows[index])), 0, 0)
        self._glyphs(points, vectors, 0.065)

    def _field(self):
        model = self.model
        values = model.fields[0]
        height = 0.36 * (model.x[-1] - model.x[0]) if model.kind == "surface" else 0
        # Expanded channel y is a declared visualization transform, not a mesh change.
        display_y = model.y * 5 if model.kind == "channel" else model.y
        if self.periodic_y:
            display_y = np.r_[display_y, model.y[0] + self.y_period]
            values = np.vstack((values, values[0]))
        points, scalars = field_points(
            model.x, display_y, np.nan_to_num(values), height=height
        )
        scalars = numpy_to_vtk(values.ravel(), deep=True)
        scalars.SetName(model.field_name)
        self.scalars = scalars
        self.mesh = quad_field(len(model.x), len(display_y), points, scalars)
        self.renderer.AddActor(_actor(self.mesh, self.table))
        self._outline(self.mesh)
        if model.kind == "surface":
            contour = vtk.vtkContourFilter()
            contour.SetInputData(self.mesh)
            low, high = self.table.GetRange()
            for index, value in enumerate(np.linspace(low, high, 13)[1:-1]):
                contour.SetValue(index, float(value))
            self.renderer.AddActor(
                _actor(contour.GetOutputPort(), color=(1, 0.78, 0.49), opacity=0.5)
            )
            self._resources.append(contour)
        if model.solid is not None:
            rows, cols = np.where(model.solid)
            cube = vtk.vtkCubeSource()
            cube.SetBounds(
                cols.min() - 0.5,
                cols.max() + 0.5,
                rows.min() - 0.5,
                rows.max() + 0.5,
                0.05,
                2,
            )
            self.renderer.AddActor(
                _actor(cube.GetOutputPort(), color=(0.16, 0.22, 0.31))
            )
            self._resources.append(cube)
        if model.vectors is not None:
            points = vtk.vtkPoints()
            vectors = vtk.vtkFloatArray()
            vectors.SetNumberOfComponents(3)
            stride = max(1, len(model.x) // 9)
            for j in range(2, len(model.y) - 2, max(1, len(model.y) // 7)):
                for i in range(2, len(model.x) - 2, stride):
                    if model.solid is not None and model.solid[j, i]:
                        continue
                    u, v = model.vectors[0, :, j, i]
                    speed = np.hypot(u, v)
                    if speed < 1e-10:
                        continue
                    points.InsertNextPoint(
                        model.x[i],
                        display_y[j],
                        height + 0.025 * (model.x[-1] - model.x[0]),
                    )
                    vectors.InsertNextTuple3(u / speed, v / speed, 0)
            if model.kind in ("surface", "channel"):
                self._glyphs(points, vectors, 0.045 * (model.x[-1] - model.x[0]))
            if model.kind == "flow":
                self._prepare_particles()
                self._streamlines(model.vectors[-1])

    def _outline(self, mesh):
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(mesh)
        edges.BoundaryEdgesOn()
        edges.FeatureEdgesOff()
        edges.NonManifoldEdgesOff()
        edges.ManifoldEdgesOff()
        actor = _actor(edges.GetOutputPort(), color=(0.53, 0.66, 0.82), opacity=0.7)
        actor.GetProperty().SetLineWidth(1.4)
        self.renderer.AddActor(actor)
        self._resources.append(edges)

    def _glyphs(self, points, vectors, scale):
        data = vtk.vtkPolyData()
        data.SetPoints(points)
        data.GetPointData().SetVectors(vectors)
        arrow = vtk.vtkArrowSource()
        arrow.SetTipResolution(12)
        arrow.SetShaftResolution(10)
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(data)
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetVectorModeToUseVector()
        glyph.SetScaleModeToDataScalingOff()
        glyph.SetScaleFactor(scale)
        glyph.OrientOn()
        self.renderer.AddActor(
            _actor(glyph.GetOutputPort(), color=(0.78, 0.92, 1), opacity=0.65)
        )
        self._resources.extend((arrow, glyph, data))

    def _streamlines(self, velocity):
        model = self.model
        image = vtk.vtkImageData()
        image.SetDimensions(len(model.x), len(model.y), 1)
        image.SetOrigin(model.x[0], model.y[0], 0)
        image.SetSpacing(model.x[1] - model.x[0], model.y[1] - model.y[0], 1)
        values = np.c_[
            velocity[0].ravel(), velocity[1].ravel(), np.zeros(velocity[0].size)
        ]
        vectors = numpy_to_vtk(values, deep=True)
        vectors.SetName("velocity")
        image.GetPointData().SetVectors(vectors)
        seeds = vtk.vtkPoints()
        for py in np.linspace(
            model.y[0] + 0.08 * np.ptp(model.y),
            model.y[-1] - 0.08 * np.ptp(model.y),
            13,
        ):
            seeds.InsertNextPoint(model.x[0] + 0.06 * np.ptp(model.x), py, 0)
        if model.steps is None:
            for px in np.linspace(0.2, 0.8, 6):
                seeds.InsertNextPoint(px, float(model.y[-1]) * 0.9, 0)
        data = vtk.vtkPolyData()
        data.SetPoints(seeds)
        trace = vtk.vtkStreamTracer()
        trace.SetInputData(image)
        trace.SetSourceData(data)
        trace.SetIntegratorTypeToRungeKutta45()
        trace.SetIntegrationDirectionToBoth()
        trace.SetMaximumPropagation(float(np.ptp(model.x)) * 8)
        trace.SetMaximumNumberOfSteps(3000)
        trace.SetComputeVorticity(False)
        actor = _actor(trace.GetOutputPort(), color=(0.65, 0.96, 0.91), opacity=0.3)
        actor.GetProperty().SetLineWidth(1.3)
        # A steady cavity may use paths from its solved field. Transient flow
        # uses current-velocity particles instead of future-frame streamlines.
        if model.steps is None:
            self.renderer.AddActor(actor)
        self._resources.extend((image, data, trace))

    def _prepare_particles(self):
        model = self.model
        rng = np.random.default_rng(8)
        count = 100 if model.steps is not None else 80
        self.particles = np.c_[
            rng.uniform(model.x[0], model.x[-1], count),
            rng.uniform(
                model.y[0],
                model.y[0] + self.y_period if self.periodic_y else model.y[-1],
                count,
            ),
        ]
        if model.solid is not None:
            # Seed only fluid sites; a visible tracer must never start in solid.
            for _ in range(100):
                ix = np.rint(self.particles[:, 0]).astype(int)
                iy = np.rint(
                    (self.particles[:, 1] - model.y[0]) / (model.y[1] - model.y[0])
                ).astype(int) % len(model.y)
                invalid = model.solid[iy, ix]
                if not invalid.any():
                    break
                self.particles[invalid, 0] = rng.uniform(
                    model.x[0], model.x[-1], invalid.sum()
                )
                self.particles[invalid, 1] = rng.uniform(
                    model.y[0],
                    model.y[0] + self.y_period if self.periodic_y else model.y[-1],
                    invalid.sum(),
                )
            else:
                raise ValueError("Could not seed tracers outside the solid mask")
        self.initial_particles = self.particles.copy()
        self.particle_points = vtk.vtkPoints()
        self.particle_points.SetDataTypeToDouble()
        for px, py in self.particles:
            self.particle_points.InsertNextPoint(px, py, 0.018 * np.ptp(model.x))
        data = vtk.vtkPolyData()
        data.SetPoints(self.particle_points)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(0.004 * np.ptp(model.x))
        sphere.SetThetaResolution(8)
        sphere.SetPhiResolution(8)
        glyph = vtk.vtkGlyph3DMapper()
        glyph.SetInputData(data)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        glyph.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.SetMapper(glyph)
        actor.GetProperty().SetColor(0.9, 1, 0.96)
        actor.GetProperty().LightingOff()
        self.renderer.AddActor(actor)
        self._resources.extend((data, sphere, glyph))

    def _velocity_at(self, positions, velocity):
        model = self.model
        px = np.clip(
            (positions[:, 0] - model.x[0]) / (model.x[1] - model.x[0]),
            0,
            len(model.x) - 1.001,
        )
        py = (positions[:, 1] - model.y[0]) / (model.y[1] - model.y[0])
        py = (
            py % len(model.y)
            if self.periodic_y
            else np.clip(py, 0, len(model.y) - 1.001)
        )
        ix, iy = px.astype(int), py.astype(int)
        fx, fy = px - ix, py - iy
        next_y = (iy + 1) % len(model.y) if self.periodic_y else iy + 1
        result = (
            velocity[:, iy, ix] * (1 - fx) * (1 - fy)
            + velocity[:, iy, ix + 1] * fx * (1 - fy)
            + velocity[:, next_y, ix] * (1 - fx) * fy
            + velocity[:, next_y, ix + 1] * fx * fy
        )
        return result.T

    def _update_particles(self, velocity, delta):
        if self.particles is None or delta <= 0:
            return
        model = self.model
        # Substep actual grid-space advection to avoid jumping through solids.
        maximum = float(np.max(np.hypot(velocity[0], velocity[1])))
        substeps = max(
            1, int(math.ceil(delta * maximum / (0.45 * (model.x[1] - model.x[0]))))
        )
        step = delta / substeps
        for _ in range(substeps):
            trial = self.particles + self._velocity_at(self.particles, velocity) * step
            if self.periodic_y:
                trial[:, 1] = (trial[:, 1] - model.y[0]) % self.y_period + model.y[0]
            invalid = (trial[:, 0] < model.x[0]) | (trial[:, 0] > model.x[-1])
            if not self.periodic_y:
                invalid |= (trial[:, 1] < model.y[0]) | (trial[:, 1] > model.y[-1])
            if model.solid is not None:
                ix = np.clip(np.rint(trial[:, 0]).astype(int), 0, len(model.x) - 1)
                iy = np.rint(
                    (trial[:, 1] - model.y[0]) / (model.y[1] - model.y[0])
                ).astype(int)
                iy = (
                    iy % len(model.y)
                    if self.periodic_y
                    else np.clip(iy, 0, len(model.y) - 1)
                )
                invalid |= model.solid[iy, ix]
            trial[invalid] = self.initial_particles[invalid]
            self.particles = trial
        for index, (px, py) in enumerate(self.particles):
            self.particle_points.SetPoint(index, px, py, 0.018 * np.ptp(model.x))
        self.particle_points.Modified()

    def _frame_camera(self):
        model = self.model
        camera = self.renderer.GetActiveCamera()
        camera.ParallelProjectionOn()
        if model.kind == "rod":
            camera.SetPosition(0.5, -1.9, 1.1)
            camera.SetFocalPoint(0.5, 0, 0)
            camera.SetViewUp(0, 0, 1)
        elif model.kind in ("flow", "channel"):
            center = ((model.x[0] + model.x[-1]) / 2, (model.y[0] + model.y[-1]) / 2, 0)
            if model.kind == "channel":
                center = (center[0], 2.5 * model.y[-1], 0)
            camera.SetPosition(center[0], center[1], max(np.ptp(model.x), 1) * 3)
            camera.SetFocalPoint(*center)
            camera.SetViewUp(0, 1, 0)
        elif model.kind == "surface":
            camera.SetPosition(2, -2.5, 2.6)
            camera.SetFocalPoint(0.5, 0.5, 0.15)
            camera.SetViewUp(0, 0, 1)
        self.renderer.ResetCamera()
        camera.Zoom(
            1.20
            if model.kind == "flow" and model.steps is None
            else 1.40
            if model.kind in ("flow", "channel")
            else 1.16
            if model.kind == "surface"
            else 1.10
        )
        self.camera_position = np.array(camera.GetPosition())
        self.camera_focal = np.array(camera.GetFocalPoint())
        self.camera_scale = camera.GetParallelScale()

    def render(self, phase):
        """phase is playback position, not fictitious simulation time."""
        phase = float(np.clip(phase, 0, 1))
        if self.phase is not None and phase < self.phase and self.particles is not None:
            self.particles = self.initial_particles.copy()
            for index, (px, py) in enumerate(self.particles):
                self.particle_points.SetPoint(
                    index, px, py, 0.018 * np.ptp(self.model.x)
                )
            self.particle_points.Modified()
            self.phase = None
        model = self.model
        camera = self.renderer.GetActiveCamera()
        angle = math.radians(7 * math.sin(phase * math.pi * 2))
        offset = self.camera_position - self.camera_focal
        if model.kind in ("surface", "rod", "generic"):
            c, s = math.cos(angle), math.sin(angle)
            offset = np.array(
                [
                    c * offset[0] - s * offset[1],
                    s * offset[0] + c * offset[1],
                    offset[2],
                ]
            )
        camera.SetPosition(*(self.camera_focal + offset))
        camera.SetParallelScale(
            self.camera_scale * (1 - 0.035 * math.sin(math.pi * phase))
        )
        if model.vectors is not None and model.kind == "flow":
            position = phase * (len(model.fields) - 1)
            lower = int(position)
            upper = min(lower + 1, len(model.fields) - 1)
            fraction = position - lower
            velocity = (
                model.vectors[lower] * (1 - fraction) + model.vectors[upper] * fraction
            )
            values = np.hypot(velocity[0], velocity[1])
            if model.solid is not None:
                values[model.solid] = np.nan
            if self.periodic_y:
                values = np.vstack((values, values[0]))
            vtk_to_numpy(self.scalars)[:] = values.ravel()
            self.scalars.Modified()
            self.mesh.Modified()
            if self.phase is not None and phase >= self.phase:
                delta = (
                    (phase - self.phase) * (model.steps[-1] - model.steps[0])
                    if model.steps is not None
                    else (phase - self.phase) * 12
                )
                self._update_particles(velocity, float(delta))
            self.phase = phase
        self.renderer.ResetCameraClippingRange()
        self.window.Render()
        self.capture.Modified()
        self.capture.Update()
        image = self.capture.GetOutput()
        width, height, _ = image.GetDimensions()
        pixels = vtk_to_numpy(image.GetPointData().GetScalars()).reshape(
            height, width, 3
        )
        return Image.fromarray(np.flipud(pixels).copy())

    def close(self):
        self.window.Finalize()


class _SceneCaptured(BaseException):
    def __init__(self, renderer):
        self.renderer = renderer


def capture_example(path, width, height):
    """Capture a selected native VTK example before its interactive event loop.

    This intentionally executes that example, just as the existing launcher
    does. Converters and examples with no native VTK render window are rejected.
    """
    real_window, real_interactor = vtk.vtkRenderWindow, vtk.vtkRenderWindowInteractor
    windows = []

    def make_window(*args, **kwargs):
        window = real_window(*args, **kwargs)
        window.SetOffScreenRendering(1)
        window.SetSize(width, height)
        windows.append(window)
        return window

    def make_interactor(*args, **kwargs):
        interactor = real_interactor(*args, **kwargs)

        def capture():
            renderers = interactor.GetRenderWindow().GetRenderers()
            if renderers.GetNumberOfItems() == 0:
                raise RuntimeError("The selected example produced no VTK renderer.")
            raise _SceneCaptured(list(renderers))

        interactor.Start = capture
        return interactor

    import matplotlib.pyplot as plt
    import sys

    original_path = list(sys.path)
    original_arguments = sys.argv[:]
    sys.argv = [str(path)]
    sys.path.insert(0, str(Path(path).parent))
    try:
        with (
            patch.object(vtk, "vtkRenderWindow", side_effect=make_window),
            patch.object(vtk, "vtkRenderWindowInteractor", side_effect=make_interactor),
            patch.object(plt, "show"),
            redirect_stdout(io.StringIO()),
        ):
            try:
                runpy.run_path(str(path), run_name="__main__")
            except _SceneCaptured as captured:
                renderers = captured.renderer
                for window in windows:
                    for renderer in renderers:
                        window.RemoveRenderer(renderer)
                # Derive display metadata from the real scene, not a fake field.
                actors = [
                    actor for renderer in renderers for actor in renderer.GetActors()
                ]
                for actor in actors:
                    if actor.GetMapper():
                        actor.GetMapper().Update()
                count = sum(
                    actor.GetMapper().GetInput().GetNumberOfPoints()
                    for actor in actors
                    if actor.GetMapper() and actor.GetMapper().GetInput()
                )
                model = Model(
                    Path(path).stem,
                    "generic",
                    np.array([0.0, 1.0]),
                    np.array([0.0, 1.0]),
                    np.zeros((1, 2, 2)),
                    "",
                    "RENDERED GEOMETRY",
                    {"source": str(path), "viewports": len(renderers)},
                    "Original example scene · camera motion is presentation",
                    [
                        ("Scene actors", str(len(actors))),
                        ("Dataset points", f"{count:,}"),
                    ],
                )
                return model, renderers
        raise ValueError(
            "This example has no native VTK scene. Use a rendering example or add a scene adapter."
        )
    finally:
        sys.path[:] = original_path
        sys.argv[:] = original_arguments
        for window in windows:
            window.Finalize()
