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


def _actor(data, lookup=None, color=(1, 1, 1), opacity=1, lit=False):
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
    prop = actor.GetProperty()
    prop.SetColor(*color)
    prop.SetOpacity(opacity)
    if lit:
        # Soft studio shading for geometry whose 3D form carries meaning.
        prop.SetInterpolationToPhong()
        prop.SetAmbient(0.38)
        prop.SetDiffuse(0.72)
        prop.SetSpecular(0.32)
        prop.SetSpecularPower(28)
    else:
        prop.LightingOff()
    return actor


def _smootherstep(value):
    value = min(1.0, max(0.0, value))
    return value * value * value * (value * (value * 6 - 15) + 10)


# Camera keyframes (azimuth°, elevation°, dolly) at each chapter's midpoint;
# the first and last entries are the opening and closing poses.
CAMERA_KEYS = (
    (-12, 8, 0.92),
    (-7, 4, 1.0),
    (6, -2, 1.03),
    (-5, 5, 1.06),
    (8, 0, 1.02),
    (-6, -3, 1.05),
    (2, 3, 1.0),
    (0, 5, 0.96),
)
# Initial view direction (from focal point toward camera), roll, and motion scale.
VIEWS = {
    "rod": ((0.0, -1.0, 0.62), 14, 0.8),
    "channel": ((0.0, -1.0, 0.95), 9, 0.6),
    "flow": ((0.0, -0.62, 1.0), 0, 0.6),
    "cavity": ((0.0, -0.32, 1.0), 0, 0.6),
    "surface": ((1.5, -3.0, 2.45), 0, 1.0),
}


class SceneRenderer:
    """One reusable offscreen scene; captions are composed separately.

    Frames are rendered at ``supersample`` times the requested size on a
    transparent background, then box-filtered down. :meth:`render` returns a
    premultiplied RGBA image, so the layout can float the scene over its own
    backdrop and derive a glow from it.
    """

    def __init__(self, model, width, height, generic_renderer=None, supersample=2, chapters=6):
        self.model = model
        self.width, self.height = int(width), int(height)
        self.ss = max(1, int(supersample))
        self.chapters = max(1, int(chapters))
        # Recorded LBM snapshots are periodic over all NY rows, including the
        # interval from the last row back to the first.
        self.periodic_y = model.kind == "flow" and model.steps is not None
        self.y_period = (
            len(model.y) * (model.y[1] - model.y[0]) if self.periodic_y else None
        )
        renderers = generic_renderer or [vtk.vtkRenderer()]
        self.renderer = renderers[0]
        self.generic = generic_renderer is not None
        if not self.generic or len(renderers) == 1:
            self.renderer.GradientBackgroundOff()
            self.renderer.SetBackground(0, 0, 0)
            self.renderer.SetBackgroundAlpha(0.0)
        self.window = vtk.vtkRenderWindow()
        self.window.SetOffScreenRendering(1)
        self.window.SetAlphaBitPlanes(1)
        self.window.SetMultiSamples(0)
        self.window.SetSize(self.width * self.ss, self.height * self.ss)
        self.window.AddRenderer(self.renderer)
        for renderer in renderers[1:]:
            self.window.AddRenderer(renderer)
            renderer.GetActiveCamera().ParallelProjectionOn()
            renderer.ResetCamera()
        self.capture = vtk.vtkWindowToImageFilter()
        self.capture.SetInput(self.window)
        self.capture.SetInputBufferTypeToRGBA()
        self.capture.ReadFrontBufferOff()
        self.capture.ShouldRerenderOff()
        self.colors = color_samples("viridis" if model.kind == "flow" else "inferno")
        self.table = lookup_table(model.fields, self.colors)
        self.mesh = None
        self.scalars = None
        self.phase = None
        self.legend = (self.table.GetRange(), model.unit)
        self.particles = None
        self.fit_points = []
        self._resources = []
        if self.generic:
            # Existing example captions are replaced by the reel's own layout.
            # Keep scientific legends and every original viewport.
            for renderer in renderers:
                if len(renderers) == 1:
                    renderer.SetViewport(0, 0, 1, 1)
                for item in list(renderer.GetActors2D()):
                    if not item.IsA("vtkScalarBarActor"):
                        renderer.RemoveActor2D(item)
            self.legend = None
        else:
            kit = vtk.vtkLightKit()
            kit.SetKeyLightIntensity(0.9)
            kit.AddLightsToRenderer(self.renderer)
            if model.kind == "rod":
                self._rod()
            else:
                self._field()
        self._frame_camera()

    @property
    def view(self):
        if self.model.kind == "flow" and self.model.steps is None:
            return "cavity"
        return self.model.kind

    def _rod(self):
        model = self.model
        module = example_module(model.key)
        # Rod radius is illustrative; a slightly thicker rod reads on a phone.
        radius = 0.072
        self.mesh, _ = module.create_cylinder_rod(
            model.x, model.fields[0, 0], radius=radius, resolution=48
        )
        self.renderer.AddActor(_actor(self.mesh, self.table, lit=True))
        self._outline(self.mesh)
        for end in (0, -1):
            cap = module.create_end_cap(
                model.x[end], model.fields[0, 0, end], radius=radius, resolution=48
            )
            self.renderer.AddActor(_actor(cap, self.table, lit=True))
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
            points.InsertNextPoint(model.x[index], -0.17, 0.0)
            vectors.InsertNextTuple3(float(np.sign(arrows[index])), 0, 0)
        self._glyphs(points, vectors, 0.07)
        # Arrow tips reach past the rod ends; keep them inside the frame.
        self.fit_points.append(
            np.array([[model.x[0] - 0.08, -0.2, 0], [model.x[-1] + 0.08, -0.2, 0]])
        )

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
        self.surface_height = height
        # Animated surfaces share one height scale across every snapshot.
        self.animated_surface = model.kind == "surface" and len(model.fields) > 1
        scalars = numpy_to_vtk(values.ravel(), deep=True)
        scalars.SetName(model.field_name)
        self.scalars = scalars
        self.mesh = quad_field(len(model.x), len(display_y), points, scalars)
        if self.animated_surface:
            self._set_surface(values)
            corners = [
                [px, py, height]
                for px in (model.x[0], model.x[-1])
                for py in (model.y[0], model.y[-1])
            ]
            self.fit_points.append(np.array(corners))
        if model.kind == "surface":
            normals = vtk.vtkPolyDataNormals()
            normals.SetInputData(self.mesh)
            normals.SplittingOff()
            self.renderer.AddActor(_actor(normals.GetOutputPort(), self.table, lit=True))
            self._resources.append(normals)
        else:
            self.renderer.AddActor(_actor(self.mesh, self.table))
        self._outline(self.mesh)
        span = float(np.ptp(model.x))
        if model.kind == "surface":
            contour = vtk.vtkContourFilter()
            contour.SetInputData(self.mesh)
            low, high = self.table.GetRange()
            for index, value in enumerate(np.linspace(low, high, 13)[1:-1]):
                contour.SetValue(index, float(value))
            lines = _actor(contour.GetOutputPort(), color=(1, 0.86, 0.62), opacity=0.55)
            lines.GetProperty().SetLineWidth(1.3 * self.ss)
            self.renderer.AddActor(lines)
            self._resources.append(contour)
        if model.solid is not None:
            rows, cols = np.where(model.solid)
            cube = vtk.vtkCubeSource()
            top = 0.05 * span
            cube.SetBounds(
                cols.min() - 0.5,
                cols.max() + 0.5,
                rows.min() - 0.5,
                rows.max() + 0.5,
                0.0,
                top,
            )
            block = _actor(cube.GetOutputPort(), color=(0.2, 0.27, 0.38), lit=True)
            self.renderer.AddActor(block)
            self._resources.append(cube)
            self.solid_box = (cols.min() - 0.5, cols.max() + 0.5, rows.min() - 0.5, rows.max() + 0.5, top)
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
                        height + 0.025 * span,
                    )
                    vectors.InsertNextTuple3(u / speed, v / speed, 0)
            if model.kind in ("surface", "channel"):
                self._glyphs(points, vectors, 0.045 * span)
            if model.kind == "flow":
                self._prepare_particles()
                self._streamlines(model.vectors[-1])

    def _set_surface(self, values):
        low, high = float(np.nanmin(self.model.fields)), float(np.nanmax(self.model.fields))
        span = high - low if high > low else 1.0
        vtk_to_numpy(self.scalars)[:] = values.ravel()
        points = vtk_to_numpy(self.mesh.GetPoints().GetData())
        points[:, 2] = (values.ravel() - low) / span * self.surface_height
        self.mesh.GetPoints().Modified()
        self.scalars.Modified()
        self.mesh.Modified()

    def _outline(self, mesh):
        edges = vtk.vtkFeatureEdges()
        edges.SetInputData(mesh)
        edges.BoundaryEdgesOn()
        edges.FeatureEdgesOff()
        edges.NonManifoldEdgesOff()
        edges.ManifoldEdgesOff()
        actor = _actor(edges.GetOutputPort(), color=(0.62, 0.76, 0.92), opacity=0.55)
        actor.GetProperty().SetLineWidth(1.2 * self.ss)
        self.renderer.AddActor(actor)
        self._resources.append(edges)

    def _glyphs(self, points, vectors, scale):
        data = vtk.vtkPolyData()
        data.SetPoints(points)
        data.GetPointData().SetVectors(vectors)
        arrow = vtk.vtkArrowSource()
        arrow.SetTipResolution(24)
        arrow.SetShaftResolution(16)
        arrow.SetTipRadius(0.12)
        arrow.SetShaftRadius(0.035)
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(data)
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetVectorModeToUseVector()
        glyph.SetScaleModeToDataScalingOff()
        glyph.SetScaleFactor(scale)
        glyph.OrientOn()
        self.renderer.AddActor(
            _actor(glyph.GetOutputPort(), color=(0.86, 0.95, 1), opacity=0.85)
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
        actor = _actor(trace.GetOutputPort(), color=(0.72, 0.98, 0.93), opacity=0.38)
        actor.GetProperty().SetLineWidth(1.2 * self.ss)
        # A steady cavity may use paths from its solved field. Transient flow
        # uses current-velocity particles instead of future-frame streamlines.
        if model.steps is None:
            self.renderer.AddActor(actor)
        self._resources.extend((image, data, trace))

    def _particle_z(self):
        return 0.012 * np.ptp(self.model.x)

    def _prepare_particles(self):
        model = self.model
        rng = np.random.default_rng(8)
        count = 140 if model.steps is not None else 110
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
        self.particle_points.SetNumberOfPoints(count)
        data = vtk.vtkPolyData()
        data.SetPoints(self.particle_points)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(0.0042 * np.ptp(model.x))
        sphere.SetThetaResolution(12)
        sphere.SetPhiResolution(12)
        glyph = vtk.vtkGlyph3DMapper()
        glyph.SetInputData(data)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.ScalingOff()
        glyph.ScalarVisibilityOff()
        actor = vtk.vtkActor()
        actor.SetMapper(glyph)
        actor.GetProperty().SetColor(1, 0.98, 0.9)
        actor.GetProperty().LightingOff()
        self.renderer.AddActor(actor)

        # Fading trails: one polyline per tracer through its recent positions.
        # Tracers move slowly in lattice units, so sample the trail sparsely.
        self.trail_length = 28
        self.trail_stride = 5
        self.trail_tick = 0
        self.inlet_rng = np.random.default_rng(9)
        self.trail_points = vtk.vtkPoints()
        self.trail_points.SetDataTypeToDouble()
        self.trail_points.SetNumberOfPoints(count * self.trail_length)
        lines = vtk.vtkCellArray()
        for index in range(count):
            lines.InsertNextCell(self.trail_length)
            for offset in range(self.trail_length):
                lines.InsertCellPoint(index * self.trail_length + offset)
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(4)
        fade = np.linspace(0, 1, self.trail_length) ** 1.6
        rgba = np.tile(
            np.c_[np.full((self.trail_length, 3), (235, 250, 255)), 200 * fade],
            (count, 1),
        ).astype(np.uint8)
        colors.DeepCopy(numpy_to_vtk(rgba, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR))
        colors.SetName("trail")
        trails = vtk.vtkPolyData()
        trails.SetPoints(self.trail_points)
        trails.SetLines(lines)
        trails.GetPointData().SetScalars(colors)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(trails)
        mapper.SetColorModeToDirectScalars()
        mapper.SetScalarModeToUsePointData()
        trail_actor = vtk.vtkActor()
        trail_actor.SetMapper(mapper)
        trail_actor.GetProperty().LightingOff()
        trail_actor.GetProperty().SetLineWidth(1.6 * self.ss)
        trail_actor.GetProperty().SetOpacity(0.999)  # route through the blended pass
        self.renderer.AddActor(trail_actor)
        self._resources.extend((data, sphere, glyph, trails, colors))
        self._reset_trails()

    def _reset_trails(self, mask=None):
        if self.particles is None:
            return
        if mask is None:
            self.history = np.repeat(self.particles[None], self.trail_length, axis=0)
        else:
            self.history[:, mask] = self.particles[mask]
        self._sync_particles()

    def _sync_particles(self):
        z = self._particle_z()
        heads = vtk_to_numpy(self.particle_points.GetData())
        heads[:, :2] = self.particles
        heads[:, 2] = z
        self.particle_points.Modified()
        trail = vtk_to_numpy(self.trail_points.GetData())
        trail[:, :2] = self.history.transpose(1, 0, 2).reshape(-1, 2)
        trail[:, 2] = z * 0.98
        self.trail_points.Modified()

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
        start = self.particles.copy()
        respawned = np.zeros(len(self.particles), dtype=bool)
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
            if model.steps is not None and invalid.any():
                # Transient channel flow: replacement tracers enter at the inlet.
                trial[invalid, 0] = model.x[0] + self.inlet_rng.uniform(0, 1.5, invalid.sum())
                trial[invalid, 1] = model.y[0] + self.inlet_rng.uniform(0, self.y_period, invalid.sum())
            respawned |= invalid
            self.particles = trial
        # Respawns and periodic wraps start a fresh trail instead of a streak.
        jump = np.abs(self.particles - start).max(axis=1) > 0.2 * np.ptp(model.x)
        self.trail_tick += 1
        if self.trail_tick % self.trail_stride == 0:
            self.history = np.roll(self.history, -1, axis=0)
        self.history[-1] = self.particles
        self._reset_trails(respawned | jump)

    def _fit_points(self):
        clouds = list(self.fit_points)
        if self.mesh is not None:
            points = vtk_to_numpy(self.mesh.GetPoints().GetData())
            clouds.append(points[:: max(1, len(points) // 6000)])
        if getattr(self, "solid_box", None) is not None:
            x0, x1, y0, y1, top = self.solid_box
            clouds.append(np.array([[x, y, top] for x in (x0, x1) for y in (y0, y1)]))
        if not clouds:
            bounds = self.renderer.ComputeVisiblePropBounds()
            clouds.append(
                np.array(
                    [
                        [x, y, z]
                        for x in bounds[0:2]
                        for y in bounds[2:4]
                        for z in bounds[4:6]
                    ]
                )
            )
        return np.vstack(clouds)

    def _projected(self, points):
        camera = self.renderer.GetActiveCamera()
        matrix = camera.GetCompositeProjectionTransformMatrix(
            self.width / self.height, -1, 1
        )
        m = np.array([[matrix.GetElement(r, c) for c in range(4)] for r in range(4)])
        clip = np.c_[points, np.ones(len(points))] @ m.T
        return clip[:, :2] / clip[:, 3:4]

    def _fit(self, fill=0.9):
        """Frame the content tightly in this scene's aspect ratio."""
        camera = self.renderer.GetActiveCamera()
        points = self._fit_points()
        aspect = self.width / self.height
        for _ in range(8):
            self.renderer.ResetCameraClippingRange()
            ndc = self._projected(points)
            low, high = ndc.min(axis=0), ndc.max(axis=0)
            center, half = (low + high) / 2, (high - low) / 2
            scale = max(half[0], half[1]) / fill
            if camera.GetParallelProjection():
                height = camera.GetParallelScale()
            else:
                height = camera.GetDistance() * math.tan(
                    math.radians(camera.GetViewAngle() / 2)
                )
            up = np.array(camera.GetViewUp())
            forward = np.array(camera.GetDirectionOfProjection())
            right = np.cross(forward, up)
            shift = right * center[0] * height * aspect + up * center[1] * height
            camera.SetFocalPoint(*(np.array(camera.GetFocalPoint()) + shift))
            camera.SetPosition(*(np.array(camera.GetPosition()) + shift))
            if camera.GetParallelProjection():
                camera.SetParallelScale(height * scale)
            else:
                camera.Dolly(1 / scale)
            if abs(scale - 1) < 1e-3 and np.abs(center).max() < 1e-3:
                break

    def _frame_camera(self):
        model = self.model
        camera = self.renderer.GetActiveCamera()
        self.motion = 1.0
        if self.generic:
            self.renderer.ResetCamera()
        else:
            direction, roll, self.motion = VIEWS[self.view]
            bounds = np.array(self.mesh.GetBounds()).reshape(3, 2)
            center = bounds.mean(axis=1)
            size = float(np.ptp(bounds, axis=1).max())
            direction = np.array(direction) / np.linalg.norm(direction)
            camera.ParallelProjectionOff()
            camera.SetViewAngle(26)
            camera.SetFocalPoint(*center)
            camera.SetPosition(*(center + direction * size * 3))
            camera.SetViewUp(0, 0, 1)
            camera.OrthogonalizeViewUp()
            camera.Roll(roll)
        self._unfitted = (
            camera.GetPosition(),
            camera.GetFocalPoint(),
            camera.GetViewUp(),
            camera.GetParallelScale(),
        )
        self._fit_base()

    def _fit_base(self):
        camera = self.renderer.GetActiveCamera()
        position, focal, up, scale = self._unfitted
        camera.SetPosition(*position)
        camera.SetFocalPoint(*focal)
        camera.SetViewUp(*up)
        camera.SetParallelScale(scale)
        if not self.generic or self.window.GetRenderers().GetNumberOfItems() == 1:
            self._fit(self.FILL)
        self.camera_position = np.array(camera.GetPosition())
        self.camera_focal = np.array(camera.GetFocalPoint())
        self.camera_up = np.array(camera.GetViewUp())
        self.camera_scale = camera.GetParallelScale()

    FILL = 0.86

    def preferred_height(self, minimum, maximum):
        """Box height that lets the framed content fill the full scene width."""
        if self.generic:
            return self.height
        ndc = self._projected(self._fit_points())
        span = ndc.max(axis=0) - ndc.min(axis=0)
        aspect = (span[0] * self.width) / max(1e-9, span[1] * self.height)
        return int(np.clip(round(self.width / aspect * 1.04), minimum, maximum))

    def set_size(self, width, height):
        self.width, self.height = int(width), int(height)
        self.window.SetSize(self.width * self.ss, self.height * self.ss)
        self._fit_base()

    def _pose(self, phase):
        keys = CAMERA_KEYS
        times = [0.0] + [(index + 0.5) / self.chapters for index in range(self.chapters)] + [1.0]
        values = [keys[0]] + [keys[1 + index % (len(keys) - 2)] for index in range(self.chapters)] + [keys[-1]]
        for index in range(len(times) - 1):
            if phase <= times[index + 1] or index == len(times) - 2:
                span = max(1e-9, times[index + 1] - times[index])
                blend = _smootherstep((phase - times[index]) / span)
                return [
                    a + (b - a) * blend
                    for a, b in zip(values[index], values[index + 1])
                ]

    def _move_camera(self, phase):
        azimuth, elevation, dolly = self._pose(phase)
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(*self.camera_position)
        camera.SetFocalPoint(*self.camera_focal)
        camera.SetViewUp(*self.camera_up)
        camera.Azimuth(azimuth * self.motion)
        camera.Elevation(elevation * self.motion)
        camera.OrthogonalizeViewUp()
        dolly = 1 + (dolly - 1) * self.motion
        if camera.GetParallelProjection():
            camera.SetParallelScale(self.camera_scale / dolly)
        else:
            camera.Dolly(dolly)

    def render(self, phase):
        """phase is playback position, not fictitious simulation time."""
        phase = float(np.clip(phase, 0, 1))
        if self.phase is not None and phase < self.phase and self.particles is not None:
            self.particles = self.initial_particles.copy()
            self._reset_trails()
            self.phase = None
        model = self.model
        self._move_camera(phase)
        if getattr(self, "animated_surface", False):
            position = phase * (len(model.fields) - 1)
            lower = int(position)
            upper = min(lower + 1, len(model.fields) - 1)
            fraction = position - lower
            self._set_surface(
                model.fields[lower] * (1 - fraction) + model.fields[upper] * fraction
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
            height, width, 4
        )
        frame = Image.fromarray(np.flipud(pixels).copy(), "RGBA")
        return frame.reduce(self.ss) if self.ss > 1 else frame

    def detail(self, scene):
        """A real crop around the obstacle, located with VTK's world/display map."""
        if self.model.solid is None:
            return None
        rows, cols = np.where(self.model.solid)
        halfx = (cols.max() - cols.min() + 1) * 3.2
        halfy = (rows.max() - rows.min() + 1) * 1.65
        cx, cy = (cols.max() + cols.min()) / 2, (rows.max() + rows.min()) / 2
        x0, x1 = self.model.x[0] + 1, self.model.x[-1] - 1
        halfx = min(halfx, (x1 - x0) / 2)
        cx = min(max(cx, x0 + halfx), x1 - halfx)
        pixels = []
        for x in (cx - halfx, cx + halfx):
            for y in (cy - halfy, cy + halfy):
                self.renderer.SetWorldPoint(float(x), float(y), 0, 1)
                self.renderer.WorldToDisplay()
                px, py, _ = self.renderer.GetDisplayPoint()
                pixels.append((px / self.ss, scene.height - py / self.ss))
        left = max(0, min(p[0] for p in pixels))
        right = min(scene.width, max(p[0] for p in pixels))
        top = max(0, min(p[1] for p in pixels))
        bottom = min(scene.height, max(p[1] for p in pixels))
        return (
            scene.crop((round(left), round(top), round(right), round(bottom)))
            if right > left and bottom > top
            else None
        )

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
