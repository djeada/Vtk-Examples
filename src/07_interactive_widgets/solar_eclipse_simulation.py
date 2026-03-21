"""
Interactive Solar Eclipse Simulation

A desktop VTK application that simulates solar eclipses in 3D with full
interactive control. The simulation models Sun, Earth, and Moon with orbital
mechanics, shadow cone geometry (umbra and penumbra), and real-time eclipse
detection.

Physical Model:
---------------
The simulation uses a simplified but physically plausible model:
- Sun is fixed at the origin as both a light source and emissive sphere
- Earth orbits the Sun in the ecliptic (X-Y) plane
- Moon orbits Earth with configurable inclination (~5 degrees real value)
- All positions are computed deterministically from simulation time
- Body radii are exaggerated relative to orbital distances for clarity

Orbital Mechanics:
------------------
Positions are computed analytically from time using circular orbits:
    Earth:  x = Re * cos(ωe * t),  y = Re * sin(ωe * t)
    Moon:   offset from Earth with inclination applied via rotation matrix

    where Re is Earth's orbital radius, ωe is Earth's angular velocity,
    and t is simulation time in days.

Eclipse Geometry:
-----------------
Solar eclipses occur when the Moon passes between the Sun and Earth near
the ecliptic plane. The shadow geometry uses similar triangles:

    Umbra length:   Lu = Rm * Dsm / (Rs - Rm)
    Penumbra angle: computed from external tangent lines

where Rs = Sun radius, Rm = Moon radius, Dsm = Sun-Moon distance.

The eclipse type is determined by whether the umbra/penumbra cones
intersect the Earth sphere:
    - Total:   umbra cone reaches and intersects Earth
    - Partial: only penumbra intersects Earth
    - None:    no shadow intersection

Interactive Controls:
---------------------
Six VTK slider widgets control simulation parameters:
    1. Time         — scrub through one full lunar orbit (0-29.53 days)
    2. Speed        — animation playback speed multiplier
    3. Scale        — body radius exaggeration factor
    4. Inclination  — Moon orbital tilt (0-10 degrees)
    5. Shadow Alpha — shadow cone transparency
    6. Earth Spin   — Earth axial rotation speed

Keyboard shortcuts provide camera presets and toggle controls:
    1-5     Camera presets (solar system, Earth, Moon, alignment, shadow)
    O       Toggle orbit path visibility
    L       Toggle label visibility
    R       Toggle sun ray visualization
    Space   Pause/resume animation
    Q/Esc   Quit

Architecture:
-------------
The code follows a clean Model-View-Controller separation:
    A. Simulation Model  — physics, time, state (CelestialBody, OrbitModel,
                           EclipseModel, SimulationClock)
    B. Rendering Layer   — VTK actors, mappers, pipelines (BodyView,
                           OrbitView, ShadowView, SceneController)
    C. Interaction Layer — widgets and callbacks (SliderFactory,
                           KeyboardController)

References:
-----------
1. Meeus, J. (1991). "Astronomical Algorithms". Willmann-Bell.
2. Espenak, F. & Meeus, J. (2006). "Five Millennium Canon of Solar Eclipses".
   NASA Technical Publication TP-2006-214141.
3. VTK Documentation: https://vtk.org/doc/nightly/html/
"""

import math
from dataclasses import dataclass, field

import numpy as np
import vtk


# =============================================================================
# Configuration Constants
# =============================================================================

# Scaled distances (not real — optimized for visualization clarity)
# The key physical constraint for eclipses is that the Moon and Sun subtend
# nearly equal angular sizes from Earth. We enforce:
#   Rs / (De - Dme) ≈ Rm / Dme
# so the umbra cone barely reaches Earth, enabling total eclipses.
SUN_RADIUS = 8.0
EARTH_RADIUS = 2.0
MOON_RADIUS = 1.2

EARTH_ORBIT_RADIUS = 60.0
MOON_ORBIT_RADIUS = 8.0

# Orbital periods in days
EARTH_ORBITAL_PERIOD = 365.25
MOON_ORBITAL_PERIOD = 29.53

# Derived angular velocities (radians per day)
EARTH_ANGULAR_VELOCITY = 2.0 * math.pi / EARTH_ORBITAL_PERIOD
MOON_ANGULAR_VELOCITY = 2.0 * math.pi / MOON_ORBITAL_PERIOD

# Earth axial tilt (radians)
EARTH_AXIAL_TILT = math.radians(23.44)

# Moon orbital inclination to ecliptic (radians)
MOON_INCLINATION_DEFAULT = math.radians(5.145)

# VTK sphere resolution
SPHERE_THETA = 64
SPHERE_PHI = 64

# Orbit path resolution (number of line segments)
ORBIT_PATH_RESOLUTION = 360

# Window dimensions
WINDOW_WIDTH = 1600
WINDOW_HEIGHT = 900

# Animation timer interval in milliseconds
TIMER_INTERVAL_MS = 33  # ~30 FPS

# Shadow cone resolution
SHADOW_CONE_RESOLUTION = 64

# Number of sun rays to draw
NUM_SUN_RAYS = 12


# =============================================================================
# Simulation Model
# =============================================================================


@dataclass
class CelestialBodyState:
    """State of a celestial body at a given instant."""

    name: str
    radius: float
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    rotation_axis: np.ndarray = field(
        default_factory=lambda: np.array([0.0, 0.0, 1.0])
    )
    rotation_angle: float = 0.0


@dataclass
class OrbitParams:
    """Orbital parameters for a celestial body."""

    radius: float
    angular_speed: float
    inclination: float = 0.0
    phase_offset: float = 0.0


@dataclass
class EclipseState:
    """Current eclipse geometry and classification."""

    aligned: bool = False
    penumbra_hits_earth: bool = False
    umbra_hits_earth: bool = False
    eclipse_type: str = "None"
    shadow_axis: np.ndarray = field(default_factory=lambda: np.zeros(3))
    umbra_length: float = 0.0
    penumbra_half_angle: float = 0.0
    umbra_half_angle: float = 0.0
    alignment_angle: float = 0.0
    moon_earth_sun_angle: float = 0.0
    centerline_offset: float = 0.0
    penumbra_radius_at_earth: float = 0.0
    umbra_radius_at_earth: float = 0.0


class SimulationClock:
    """
    Single source of truth for simulation time.

    All positions are recomputed from this time value to ensure deterministic
    behavior and stable slider scrubbing with no accumulated drift.
    """

    def __init__(self):
        self.time = 0.0
        self.speed = 1.0
        self.paused = False

    def set_time(self, t: float) -> None:
        self.time = t

    def advance(self, dt: float) -> None:
        if not self.paused:
            self.time += dt * self.speed

    def set_speed(self, s: float) -> None:
        self.speed = s

    def toggle_pause(self) -> None:
        self.paused = not self.paused


class OrbitModel:
    """Computes orbital positions from time and parameters."""

    @staticmethod
    def compute_position(
        center: np.ndarray, orbit: OrbitParams, time: float
    ) -> np.ndarray:
        """
        Compute the position of a body on a circular orbit at the given time.

        The orbit lies in the X-Y plane, optionally tilted by the inclination
        angle around the X-axis. The position is offset from the center point.

        Args:
            center: The point the body orbits around.
            orbit: Orbital parameters (radius, speed, inclination, phase).
            time: Current simulation time in days.

        Returns:
            3D position as numpy array.
        """
        angle = orbit.angular_speed * time + orbit.phase_offset

        # Position in the orbital plane (before inclination)
        x = orbit.radius * math.cos(angle)
        y = orbit.radius * math.sin(angle)
        z = 0.0

        # Apply inclination (rotation about X-axis)
        if abs(orbit.inclination) > 1e-10:
            cos_i = math.cos(orbit.inclination)
            sin_i = math.sin(orbit.inclination)
            y_rot = y * cos_i - z * sin_i
            z_rot = y * sin_i + z * cos_i
            y = y_rot
            z = z_rot

        return center + np.array([x, y, z])


class EclipseModel:
    """
    Computes eclipse geometry from body positions and radii.

    Uses similar-triangle geometry to build umbra and penumbra cones from the
    Sun through the Moon, then tests intersection with the Earth sphere.
    """

    @staticmethod
    def compute(
        sun_pos: np.ndarray,
        sun_radius: float,
        earth_pos: np.ndarray,
        earth_radius: float,
        moon_pos: np.ndarray,
        moon_radius: float,
    ) -> EclipseState:
        state = EclipseState()

        sun_to_moon = moon_pos - sun_pos
        sun_moon_dist = np.linalg.norm(sun_to_moon)
        if sun_moon_dist < 1e-10:
            return state

        shadow_dir = sun_to_moon / sun_moon_dist

        # Moon to Earth vector
        moon_to_earth = earth_pos - moon_pos
        moon_earth_dist = np.linalg.norm(moon_to_earth)

        # Alignment angle: angle between shadow axis and Moon-to-Earth direction
        if moon_earth_dist > 1e-10:
            cos_angle = np.clip(
                np.dot(shadow_dir, moon_to_earth / moon_earth_dist), -1.0, 1.0
            )
            state.alignment_angle = math.degrees(math.acos(cos_angle))
        else:
            state.alignment_angle = 0.0

        # Moon-Earth-Sun angle (how close Moon is to being between Sun and Earth)
        sun_to_earth = earth_pos - sun_pos
        sun_earth_dist = np.linalg.norm(sun_to_earth)
        earth_to_moon_dir = moon_pos - earth_pos
        earth_to_sun_dir = sun_pos - earth_pos
        if moon_earth_dist > 1e-10 and sun_earth_dist > 1e-10:
            cos_mes = np.clip(
                np.dot(
                    earth_to_moon_dir / moon_earth_dist,
                    earth_to_sun_dir / np.linalg.norm(earth_to_sun_dir),
                ),
                -1.0,
                1.0,
            )
            state.moon_earth_sun_angle = math.degrees(math.acos(cos_mes))

        state.shadow_axis = shadow_dir

        # Umbra cone: converging shadow (similar triangles)
        # The umbra tip is beyond the Moon where the apparent Sun disk
        # is fully blocked.
        if sun_radius > moon_radius:
            state.umbra_length = moon_radius * sun_moon_dist / (sun_radius - moon_radius)
        else:
            state.umbra_length = 1e6  # Moon >= Sun means infinite umbra

        # Umbra half-angle
        if state.umbra_length > 1e-10:
            state.umbra_half_angle = math.atan2(moon_radius, state.umbra_length)

        # Penumbra half-angle (diverging cone from external tangents)
        if sun_moon_dist > 1e-10:
            state.penumbra_half_angle = math.atan2(
                sun_radius + moon_radius, sun_moon_dist
            )

        # Project Earth center onto shadow axis to find closest approach
        moon_to_earth_vec = earth_pos - moon_pos
        proj_dist = np.dot(moon_to_earth_vec, shadow_dir)

        if proj_dist > 0:
            # Earth is in the forward shadow direction
            closest_point = moon_pos + proj_dist * shadow_dir
            lateral_dist = np.linalg.norm(earth_pos - closest_point)
            state.centerline_offset = lateral_dist

            # Umbra radius at Earth's distance along shadow axis
            if state.umbra_length > 1e-10 and proj_dist < state.umbra_length:
                umbra_radius_at_earth = moon_radius * (
                    1.0 - proj_dist / state.umbra_length
                )
            else:
                umbra_radius_at_earth = 0.0
            state.umbra_radius_at_earth = umbra_radius_at_earth

            # Penumbra radius at Earth's distance
            penumbra_radius_at_earth = (
                moon_radius + proj_dist * math.tan(state.penumbra_half_angle)
            )
            state.penumbra_radius_at_earth = penumbra_radius_at_earth

            # Check intersection with Earth sphere
            if lateral_dist < penumbra_radius_at_earth + earth_radius:
                state.penumbra_hits_earth = True
                state.aligned = True

            if (
                umbra_radius_at_earth > 0
                and lateral_dist < umbra_radius_at_earth + earth_radius * 0.3
            ):
                state.umbra_hits_earth = True

            # Classify eclipse
            if state.umbra_hits_earth:
                state.eclipse_type = "Total Solar Eclipse"
            elif state.penumbra_hits_earth:
                state.eclipse_type = "Partial Solar Eclipse"
            else:
                state.eclipse_type = "No Eclipse"
        else:
            state.eclipse_type = "No Eclipse"

        return state


# =============================================================================
# Rendering Layer — VTK Actors and Views
# =============================================================================


class BodyView:
    """Creates and manages VTK actors for a celestial body sphere."""

    def __init__(
        self,
        name: str,
        radius: float,
        color: tuple[float, float, float],
        label_color: tuple[float, float, float] = (1.0, 1.0, 1.0),
        emissive: bool = False,
        opacity: float = 1.0,
        glow_color: tuple[float, float, float] | None = None,
        glow_scale: float = 1.18,
        glow_opacity: float = 0.18,
    ):
        self.name = name
        self.base_radius = radius
        self.scale_factor = 1.0
        self.glow_scale = glow_scale

        self.source = vtk.vtkSphereSource()
        self.source.SetRadius(radius)
        self.source.SetThetaResolution(SPHERE_THETA)
        self.source.SetPhiResolution(SPHERE_PHI)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.source.GetOutputPort())

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        prop = self.actor.GetProperty()
        prop.SetColor(*color)
        prop.SetOpacity(opacity)

        if emissive:
            prop.SetAmbient(1.0)
            prop.SetDiffuse(0.0)
            prop.SetSpecular(0.0)
        else:
            prop.SetAmbient(0.28)
            prop.SetDiffuse(0.95)
            prop.SetSpecular(0.5)
            prop.SetSpecularPower(38.0)
            prop.SetInterpolationToPhong()

        self.glow_actor = None
        if glow_color is not None:
            self.glow_source = vtk.vtkSphereSource()
            self.glow_source.SetRadius(radius * glow_scale)
            self.glow_source.SetThetaResolution(SPHERE_THETA)
            self.glow_source.SetPhiResolution(SPHERE_PHI)

            self.glow_mapper = vtk.vtkPolyDataMapper()
            self.glow_mapper.SetInputConnection(self.glow_source.GetOutputPort())

            self.glow_actor = vtk.vtkActor()
            self.glow_actor.SetMapper(self.glow_mapper)
            glow_prop = self.glow_actor.GetProperty()
            glow_prop.SetColor(*glow_color)
            glow_prop.SetOpacity(glow_opacity)
            glow_prop.SetAmbient(1.0)
            glow_prop.SetDiffuse(0.0)
            glow_prop.SetSpecular(0.0)

        # Label — 2D screen-space text (immune to scene lighting)
        self.label = vtk.vtkTextActor()
        self.label.SetInput(name)
        tp = self.label.GetTextProperty()
        tp.SetFontSize(15)
        tp.SetColor(*label_color)
        tp.SetJustificationToCentered()
        tp.SetVerticalJustificationToBottom()
        tp.SetBold(True)
        tp.SetShadow(True)
        tp.SetShadowOffset(1, 1)
        tp.SetFontFamilyToArial()
        self.label_world_pos = np.zeros(3)

    def set_position(self, pos: np.ndarray, label_z_offset: float = 0.0) -> None:
        self.actor.SetPosition(*pos)
        if self.glow_actor is not None:
            self.glow_actor.SetPosition(*pos)
        offset = self.base_radius * self.scale_factor * 1.6 + label_z_offset
        self.label_world_pos = np.array([pos[0], pos[1], pos[2] + offset])

    def project_label(self, renderer) -> None:
        """Project the 3D label position to 2D screen coordinates."""
        coord = vtk.vtkCoordinate()
        coord.SetCoordinateSystemToWorld()
        coord.SetValue(*self.label_world_pos)
        disp = coord.GetComputedDisplayValue(renderer)
        self.label.SetPosition(disp[0], disp[1])

    def set_scale(self, factor: float) -> None:
        self.scale_factor = factor
        self.source.SetRadius(self.base_radius * factor)
        if self.glow_actor is not None:
            self.glow_source.SetRadius(self.base_radius * factor * self.glow_scale)

    def set_rotation(self, axis: np.ndarray, angle_deg: float) -> None:
        self.actor.SetOrientation(0, 0, 0)
        if np.linalg.norm(axis) > 1e-10:
            self.actor.RotateWXYZ(angle_deg, *axis)
        if self.glow_actor is not None:
            self.glow_actor.SetOrientation(0, 0, 0)
            if np.linalg.norm(axis) > 1e-10:
                self.glow_actor.RotateWXYZ(angle_deg, *axis)


class OrbitView:
    """Creates and manages VTK actors for an orbital path."""

    def __init__(self, color: tuple[float, float, float], opacity: float = 0.4):
        self.points = vtk.vtkPoints()
        self.lines = vtk.vtkCellArray()
        self.polydata = vtk.vtkPolyData()
        self.color = color
        self.opacity = opacity

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputData(self.polydata)

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        prop = self.actor.GetProperty()
        prop.SetColor(*color)
        prop.SetOpacity(opacity)
        prop.SetLineWidth(1.5)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)

    def update_path(
        self, center: np.ndarray, orbit: OrbitParams, n_points: int = ORBIT_PATH_RESOLUTION
    ) -> None:
        """Recompute the full orbit path polyline."""
        self.points = vtk.vtkPoints()
        self.lines = vtk.vtkCellArray()

        period = 2.0 * math.pi / orbit.angular_speed if orbit.angular_speed > 0 else 1.0
        dt = period / n_points

        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(n_points + 1)

        for i in range(n_points + 1):
            t = i * dt
            pos = OrbitModel.compute_position(center, orbit, t)
            self.points.InsertNextPoint(*pos)
            polyline.GetPointIds().SetId(i, i)

        self.lines.InsertNextCell(polyline)
        self.polydata.SetPoints(self.points)
        self.polydata.SetLines(self.lines)
        self.polydata.Modified()


class ShadowView:
    """
    Creates and manages VTK actors for umbra and penumbra shadow cones.

    The cones are built as custom polydata meshes from the Moon position
    extending along the shadow axis.
    """

    def __init__(self):
        # Umbra cone (dark, converging)
        self.umbra_polydata = vtk.vtkPolyData()
        self.umbra_mapper = vtk.vtkPolyDataMapper()
        self.umbra_mapper.SetInputData(self.umbra_polydata)
        self.umbra_actor = vtk.vtkActor()
        self.umbra_actor.SetMapper(self.umbra_mapper)
        prop = self.umbra_actor.GetProperty()
        prop.SetColor(0.05, 0.08, 0.16)
        prop.SetOpacity(0.5)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)

        # Penumbra cone (lighter, diverging)
        self.penumbra_polydata = vtk.vtkPolyData()
        self.penumbra_mapper = vtk.vtkPolyDataMapper()
        self.penumbra_mapper.SetInputData(self.penumbra_polydata)
        self.penumbra_actor = vtk.vtkActor()
        self.penumbra_actor.SetMapper(self.penumbra_mapper)
        prop = self.penumbra_actor.GetProperty()
        prop.SetColor(0.85, 0.62, 0.16)
        prop.SetOpacity(0.2)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)

    def set_opacity(self, alpha: float) -> None:
        self.umbra_actor.GetProperty().SetOpacity(alpha)
        self.penumbra_actor.GetProperty().SetOpacity(alpha * 0.5)

    def update(
        self,
        moon_pos: np.ndarray,
        moon_radius: float,
        shadow_dir: np.ndarray,
        eclipse: EclipseState,
    ) -> None:
        """Rebuild shadow cone meshes from current eclipse geometry."""
        if not eclipse.aligned or np.linalg.norm(shadow_dir) < 1e-10:
            self._clear_polydata(self.umbra_polydata)
            self._clear_polydata(self.penumbra_polydata)
            return

        # Build umbra cone (converging from Moon radius to tip)
        umbra_length = min(eclipse.umbra_length, 80.0)
        if umbra_length > 0.1:
            self._build_cone(
                self.umbra_polydata,
                moon_pos,
                shadow_dir,
                moon_radius,
                0.01,
                umbra_length,
            )
        else:
            self._clear_polydata(self.umbra_polydata)

        # Build penumbra cone (diverging from Moon)
        penumbra_length = min(umbra_length * 1.5, 80.0)
        if eclipse.penumbra_half_angle > 0:
            end_radius = moon_radius + penumbra_length * math.tan(
                eclipse.penumbra_half_angle
            )
            self._build_cone(
                self.penumbra_polydata,
                moon_pos,
                shadow_dir,
                moon_radius * 1.2,
                end_radius,
                penumbra_length,
            )
        else:
            self._clear_polydata(self.penumbra_polydata)

    @staticmethod
    def _build_cone(
        polydata: vtk.vtkPolyData,
        origin: np.ndarray,
        direction: np.ndarray,
        start_radius: float,
        end_radius: float,
        length: float,
    ) -> None:
        """
        Build a truncated cone mesh as VTK polydata.

        The cone starts at origin with start_radius and extends along direction
        for the given length, ending at end_radius.
        """
        n = SHADOW_CONE_RESOLUTION
        direction = direction / np.linalg.norm(direction)

        # Build a local coordinate frame perpendicular to direction
        if abs(direction[0]) < 0.9:
            up = np.array([1.0, 0.0, 0.0])
        else:
            up = np.array([0.0, 1.0, 0.0])
        u = np.cross(direction, up)
        u = u / np.linalg.norm(u)
        v = np.cross(direction, u)

        points = vtk.vtkPoints()
        triangles = vtk.vtkCellArray()

        end_point = origin + direction * length

        # Generate ring vertices at start and end
        for i in range(n):
            angle = 2.0 * math.pi * i / n
            cos_a = math.cos(angle)
            sin_a = math.sin(angle)

            # Start ring
            p_start = origin + start_radius * (cos_a * u + sin_a * v)
            points.InsertNextPoint(*p_start)

            # End ring
            p_end = end_point + end_radius * (cos_a * u + sin_a * v)
            points.InsertNextPoint(*p_end)

        # Connect rings with triangles
        for i in range(n):
            i0 = 2 * i
            i1 = 2 * i + 1
            i2 = 2 * ((i + 1) % n)
            i3 = 2 * ((i + 1) % n) + 1

            tri1 = vtk.vtkTriangle()
            tri1.GetPointIds().SetId(0, i0)
            tri1.GetPointIds().SetId(1, i1)
            tri1.GetPointIds().SetId(2, i2)
            triangles.InsertNextCell(tri1)

            tri2 = vtk.vtkTriangle()
            tri2.GetPointIds().SetId(0, i2)
            tri2.GetPointIds().SetId(1, i1)
            tri2.GetPointIds().SetId(2, i3)
            triangles.InsertNextCell(tri2)

        polydata.SetPoints(points)
        polydata.SetPolys(triangles)
        polydata.Modified()

    @staticmethod
    def _clear_polydata(polydata: vtk.vtkPolyData) -> None:
        polydata.SetPoints(vtk.vtkPoints())
        polydata.SetPolys(vtk.vtkCellArray())
        polydata.Modified()


class SunRayView:
    """Visualizes a bundle of rays from the Sun through the Moon region."""

    def __init__(self):
        self.polydata = vtk.vtkPolyData()
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputData(self.polydata)
        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        prop = self.actor.GetProperty()
        prop.SetColor(1.0, 1.0, 0.4)
        prop.SetOpacity(0.25)
        prop.SetLineWidth(1.0)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)
        self.actor.SetVisibility(False)

    def update(
        self,
        sun_pos: np.ndarray,
        sun_radius: float,
        moon_pos: np.ndarray,
        earth_pos: np.ndarray,
    ) -> None:
        """Draw rays from Sun surface toward Moon/Earth region."""
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()

        sun_to_moon = moon_pos - sun_pos
        dist = np.linalg.norm(sun_to_moon)
        if dist < 1e-10:
            self.polydata.SetPoints(points)
            self.polydata.SetLines(lines)
            self.polydata.Modified()
            return

        direction = sun_to_moon / dist

        # Build local frame
        if abs(direction[0]) < 0.9:
            up = np.array([1.0, 0.0, 0.0])
        else:
            up = np.array([0.0, 1.0, 0.0])
        u = np.cross(direction, up)
        u = u / np.linalg.norm(u)
        v = np.cross(direction, u)

        ray_length = np.linalg.norm(earth_pos - sun_pos) * 1.2
        idx = 0

        for i in range(NUM_SUN_RAYS):
            angle = 2.0 * math.pi * i / NUM_SUN_RAYS
            cos_a = math.cos(angle)
            sin_a = math.sin(angle)

            start = sun_pos + sun_radius * 0.95 * (cos_a * u + sin_a * v)
            end = start + direction * ray_length

            points.InsertNextPoint(*start)
            points.InsertNextPoint(*end)

            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, idx)
            line.GetPointIds().SetId(1, idx + 1)
            lines.InsertNextCell(line)
            idx += 2

        # Center ray
        start = sun_pos
        end = sun_pos + direction * ray_length
        points.InsertNextPoint(*start)
        points.InsertNextPoint(*end)
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, idx)
        line.GetPointIds().SetId(1, idx + 1)
        lines.InsertNextCell(line)

        self.polydata.SetPoints(points)
        self.polydata.SetLines(lines)
        self.polydata.Modified()


class EclipseFootprintView:
    """Visualizes the eclipse footprint on Earth with a filled disk and ring."""

    def __init__(self):
        self.disk_source = vtk.vtkDiskSource()
        self.disk_source.SetInnerRadius(0.0)
        self.disk_source.SetOuterRadius(1.0)
        self.disk_source.SetRadialResolution(1)
        self.disk_source.SetCircumferentialResolution(96)

        self.disk_mapper = vtk.vtkPolyDataMapper()
        self.disk_mapper.SetInputConnection(self.disk_source.GetOutputPort())

        self.disk_actor = vtk.vtkActor()
        self.disk_actor.SetMapper(self.disk_mapper)
        prop = self.disk_actor.GetProperty()
        prop.SetColor(0.05, 0.08, 0.12)
        prop.SetOpacity(0.85)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)
        self.disk_actor.SetVisibility(False)

        self.ring_source = vtk.vtkDiskSource()
        self.ring_source.SetInnerRadius(0.85)
        self.ring_source.SetOuterRadius(1.0)
        self.ring_source.SetRadialResolution(1)
        self.ring_source.SetCircumferentialResolution(96)

        self.ring_mapper = vtk.vtkPolyDataMapper()
        self.ring_mapper.SetInputConnection(self.ring_source.GetOutputPort())

        self.ring_actor = vtk.vtkActor()
        self.ring_actor.SetMapper(self.ring_mapper)
        prop = self.ring_actor.GetProperty()
        prop.SetColor(1.0, 0.85, 0.25)
        prop.SetOpacity(0.9)
        prop.SetAmbient(1.0)
        prop.SetDiffuse(0.0)
        self.ring_actor.SetVisibility(False)

    @staticmethod
    def _orient_disk(actor: vtk.vtkActor, footprint_pos: np.ndarray, shadow_dir: np.ndarray) -> None:
        actor.SetPosition(*footprint_pos)
        actor.SetOrientation(0, 0, 0)

        z_axis = np.array([0.0, 0.0, 1.0])
        target = -shadow_dir
        dot = float(np.clip(np.dot(z_axis, target), -1.0, 1.0))
        rot_axis = np.cross(z_axis, target)
        rot_axis_len = np.linalg.norm(rot_axis)

        if rot_axis_len > 1e-10:
            rot_axis = rot_axis / rot_axis_len
            angle_deg = math.degrees(math.acos(dot))
            actor.RotateWXYZ(angle_deg, *rot_axis)
        elif dot < 0.0:
            actor.RotateWXYZ(180.0, 1.0, 0.0, 0.0)

    def update(
        self,
        earth_pos: np.ndarray,
        earth_radius: float,
        shadow_dir: np.ndarray,
        eclipse: EclipseState,
    ) -> None:
        """Position and orient the eclipse footprint on Earth's sunward face."""
        if not eclipse.aligned:
            self.disk_actor.SetVisibility(False)
            self.ring_actor.SetVisibility(False)
            return

        # Place disk on Earth surface facing the shadow
        offset = earth_radius * 1.01
        footprint_pos = earth_pos - shadow_dir * offset
        self._orient_disk(self.disk_actor, footprint_pos, shadow_dir)
        self._orient_disk(self.ring_actor, footprint_pos, shadow_dir)

        # Scale footprint based on eclipse type
        if eclipse.umbra_hits_earth:
            self.disk_source.SetOuterRadius(earth_radius * 0.23)
            self.disk_actor.GetProperty().SetColor(0.04, 0.06, 0.1)
            self.disk_actor.GetProperty().SetOpacity(0.9)
            self.ring_source.SetInnerRadius(earth_radius * 0.27)
            self.ring_source.SetOuterRadius(earth_radius * 0.34)
            self.ring_actor.GetProperty().SetColor(1.0, 0.88, 0.3)
            self.ring_actor.GetProperty().SetOpacity(0.95)
            self.disk_actor.SetVisibility(True)
            self.ring_actor.SetVisibility(True)
        elif eclipse.penumbra_hits_earth:
            self.disk_source.SetOuterRadius(earth_radius * 0.52)
            self.disk_actor.GetProperty().SetColor(0.7, 0.48, 0.1)
            self.disk_actor.GetProperty().SetOpacity(0.24)
            self.ring_source.SetInnerRadius(earth_radius * 0.58)
            self.ring_source.SetOuterRadius(earth_radius * 0.7)
            self.ring_actor.GetProperty().SetColor(1.0, 0.82, 0.26)
            self.ring_actor.GetProperty().SetOpacity(0.85)
            self.disk_actor.SetVisibility(True)
            self.ring_actor.SetVisibility(True)
        else:
            self.disk_actor.SetVisibility(False)
            self.ring_actor.SetVisibility(False)


class EclipseInsetView:
    """Shows the eclipse as seen from Earth in a dedicated inset renderer."""

    def __init__(self):
        self.sun_source = vtk.vtkRegularPolygonSource()
        self.sun_source.SetNumberOfSides(128)
        self.sun_source.GeneratePolygonOn()
        self.sun_source.SetRadius(1.0)

        self.sun_mapper = vtk.vtkPolyDataMapper()
        self.sun_mapper.SetInputConnection(self.sun_source.GetOutputPort())

        self.sun_actor = vtk.vtkActor()
        self.sun_actor.SetMapper(self.sun_mapper)
        sun_prop = self.sun_actor.GetProperty()
        sun_prop.SetColor(1.0, 0.86, 0.18)
        sun_prop.SetAmbient(1.0)
        sun_prop.SetDiffuse(0.0)

        self.sun_outline_actor = vtk.vtkActor()
        self.sun_outline_actor.SetMapper(self.sun_mapper)
        outline_prop = self.sun_outline_actor.GetProperty()
        outline_prop.SetRepresentationToWireframe()
        outline_prop.SetColor(1.0, 0.97, 0.75)
        outline_prop.SetOpacity(0.85)
        outline_prop.SetLineWidth(2.0)
        outline_prop.SetAmbient(1.0)
        outline_prop.SetDiffuse(0.0)

        self.moon_source = vtk.vtkRegularPolygonSource()
        self.moon_source.SetNumberOfSides(128)
        self.moon_source.GeneratePolygonOn()
        self.moon_source.SetRadius(1.0)

        self.moon_mapper = vtk.vtkPolyDataMapper()
        self.moon_mapper.SetInputConnection(self.moon_source.GetOutputPort())

        self.moon_actor = vtk.vtkActor()
        self.moon_actor.SetMapper(self.moon_mapper)
        moon_prop = self.moon_actor.GetProperty()
        moon_prop.SetColor(0.55, 0.59, 0.68)
        moon_prop.SetAmbient(1.0)
        moon_prop.SetDiffuse(0.0)

        self.crosshair_polydata = vtk.vtkPolyData()
        crosshair_points = vtk.vtkPoints()
        crosshair_lines = vtk.vtkCellArray()
        segments = [
            (-0.95, 0.0, 0.0, 0.95, 0.0, 0.0),
            (0.0, -0.95, 0.0, 0.0, 0.95, 0.0),
        ]
        for idx, (x1, y1, z1, x2, y2, z2) in enumerate(segments):
            p0 = idx * 2
            crosshair_points.InsertNextPoint(x1, y1, z1)
            crosshair_points.InsertNextPoint(x2, y2, z2)
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, p0)
            line.GetPointIds().SetId(1, p0 + 1)
            crosshair_lines.InsertNextCell(line)
        self.crosshair_polydata.SetPoints(crosshair_points)
        self.crosshair_polydata.SetLines(crosshair_lines)

        self.crosshair_mapper = vtk.vtkPolyDataMapper()
        self.crosshair_mapper.SetInputData(self.crosshair_polydata)

        self.crosshair_actor = vtk.vtkActor()
        self.crosshair_actor.SetMapper(self.crosshair_mapper)
        cross_prop = self.crosshair_actor.GetProperty()
        cross_prop.SetColor(0.75, 0.8, 0.95)
        cross_prop.SetOpacity(0.18)
        cross_prop.SetLineWidth(1.0)
        cross_prop.SetAmbient(1.0)
        cross_prop.SetDiffuse(0.0)

    def add_to_renderer(self, renderer: vtk.vtkRenderer) -> None:
        renderer.AddActor(self.crosshair_actor)
        renderer.AddActor(self.sun_actor)
        renderer.AddActor(self.sun_outline_actor)
        renderer.AddActor(self.moon_actor)

        camera = renderer.GetActiveCamera()
        camera.SetPosition(0.0, 0.0, 8.0)
        camera.SetFocalPoint(0.0, 0.0, 0.0)
        camera.SetViewUp(0.0, 1.0, 0.0)
        camera.ParallelProjectionOn()
        camera.SetParallelScale(1.15)

    def update(
        self,
        sun_pos: np.ndarray,
        sun_radius: float,
        earth_pos: np.ndarray,
        moon_pos: np.ndarray,
        moon_radius: float,
    ) -> None:
        earth_to_sun = sun_pos - earth_pos
        earth_to_moon = moon_pos - earth_pos
        sun_dist = np.linalg.norm(earth_to_sun)
        moon_dist = np.linalg.norm(earth_to_moon)

        if sun_dist < 1e-10 or moon_dist < 1e-10:
            return

        sun_dir = earth_to_sun / sun_dist
        moon_dir = earth_to_moon / moon_dist

        ref = np.array([0.0, 0.0, 1.0])
        if np.linalg.norm(np.cross(ref, sun_dir)) < 1e-6:
            ref = np.array([0.0, 1.0, 0.0])

        x_axis = np.cross(ref, sun_dir)
        x_axis = x_axis / np.linalg.norm(x_axis)
        y_axis = np.cross(sun_dir, x_axis)
        y_axis = y_axis / np.linalg.norm(y_axis)

        moon_plane = moon_dir - np.dot(moon_dir, sun_dir) * sun_dir
        moon_x = float(np.dot(moon_plane, x_axis))
        moon_y = float(np.dot(moon_plane, y_axis))

        sun_ang = math.asin(min(0.999999, sun_radius / sun_dist))
        moon_ang = math.asin(min(0.999999, moon_radius / moon_dist))

        max_extent = max(
            sun_ang,
            moon_ang,
            abs(moon_x) + moon_ang,
            abs(moon_y) + moon_ang,
            1e-3,
        )
        scale = 0.78 / max_extent

        self.sun_actor.SetPosition(0.0, 0.0, 0.0)
        self.sun_actor.SetScale(sun_ang * scale, sun_ang * scale, 1.0)
        self.sun_outline_actor.SetPosition(0.0, 0.0, 0.01)
        self.sun_outline_actor.SetScale(sun_ang * scale, sun_ang * scale, 1.0)
        self.moon_actor.SetPosition(moon_x * scale, moon_y * scale, 0.02)
        self.moon_actor.SetScale(moon_ang * scale, moon_ang * scale, 1.0)


# =============================================================================
# HUD and Status Display
# =============================================================================


class StatusDisplay:
    """On-screen text overlay showing eclipse state and simulation info."""

    def __init__(self):
        self.actor = vtk.vtkTextActor()
        self.actor.SetInput("Initializing...")
        self.actor.SetPosition(15, 15)
        prop = self.actor.GetTextProperty()
        prop.SetFontSize(16)
        prop.SetColor(0.9, 0.9, 0.9)
        prop.SetFontFamilyToCourier()
        prop.SetShadow(True)

        self.title_actor = vtk.vtkTextActor()
        self.title_actor.SetInput("Solar Eclipse Simulation")
        prop2 = self.title_actor.GetTextProperty()
        prop2.SetFontSize(22)
        prop2.SetColor(1.0, 0.85, 0.2)
        prop2.SetBold(True)
        prop2.SetShadow(True)

        self.help_actor = vtk.vtkTextActor()
        self.help_actor.SetInput(
            "[1-5] Camera  [O] Orbits  [L] Labels  [R] Rays  [Space] Pause  [Q] Quit"
        )
        prop3 = self.help_actor.GetTextProperty()
        prop3.SetFontSize(13)
        prop3.SetColor(0.6, 0.6, 0.6)
        prop3.SetFontFamilyToCourier()
        prop3.SetShadow(True)

        self.inset_actor = vtk.vtkTextActor()
        self.inset_actor.SetInput("Earth View Inset")
        prop4 = self.inset_actor.GetTextProperty()
        prop4.SetFontSize(14)
        prop4.SetColor(0.8, 0.84, 0.96)
        prop4.SetBold(True)
        prop4.SetShadow(True)

    def position_for_window(self, width: int, height: int) -> None:
        self.title_actor.SetPosition(15, height - 40)
        self.help_actor.SetPosition(15, height - 65)
        self.actor.SetPosition(15, 15)
        self.inset_actor.SetPosition(width - 250, height - 335)

    def update(self, clock: SimulationClock, eclipse: EclipseState) -> None:
        day = clock.time % MOON_ORBITAL_PERIOD
        status = "PAUSED" if clock.paused else f"Speed: {clock.speed:.1f}x"
        if eclipse.penumbra_hits_earth:
            offset_pct = (
                100.0 * eclipse.centerline_offset
                / max(eclipse.penumbra_radius_at_earth, 1e-6)
            )
        else:
            offset_pct = 100.0

        if eclipse.eclipse_type == "Total Solar Eclipse":
            header = "TOTAL SOLAR ECLIPSE"
            color = (1.0, 0.88, 0.3)
        elif eclipse.eclipse_type == "Partial Solar Eclipse":
            header = "PARTIAL SOLAR ECLIPSE"
            color = (1.0, 0.82, 0.42)
        else:
            header = "ALIGNMENT MISS"
            color = (0.85, 0.88, 0.95)

        text = (
            f"{header}\n"
            f"Day {day:5.2f}/{MOON_ORBITAL_PERIOD:.2f}  |  {status}\n"
            f"Axis offset: {eclipse.centerline_offset:5.2f}  "
            f"({offset_pct:5.1f}% of penumbra radius)\n"
            f"Alignment: {eclipse.alignment_angle:5.2f} deg  |  "
            f"Moon-Earth-Sun: {eclipse.moon_earth_sun_angle:5.2f} deg"
        )
        self.actor.SetInput(text)
        self.actor.GetTextProperty().SetColor(*color)


# =============================================================================
# Slider Factory
# =============================================================================


class SliderFactory:
    """Creates VTK slider widgets with consistent styling."""

    @staticmethod
    def create(
        interactor: vtk.vtkRenderWindowInteractor,
        title: str,
        min_val: float,
        max_val: float,
        initial: float,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        callback,
    ) -> vtk.vtkSliderWidget:
        slider_rep = vtk.vtkSliderRepresentation2D()
        slider_rep.SetMinimumValue(min_val)
        slider_rep.SetMaximumValue(max_val)
        slider_rep.SetValue(initial)
        slider_rep.SetTitleText(title)

        # Positioning
        slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint1Coordinate().SetValue(x1, y1)
        slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint2Coordinate().SetValue(x2, y2)

        # Compact sizing to avoid overlap
        slider_rep.SetSliderLength(0.010)
        slider_rep.SetSliderWidth(0.012)
        slider_rep.SetTubeWidth(0.003)
        slider_rep.SetEndCapLength(0.005)
        slider_rep.SetEndCapWidth(0.010)
        slider_rep.SetTitleHeight(0.012)
        slider_rep.SetLabelHeight(0.010)

        # High-contrast colors on dark background
        slider_rep.GetSliderProperty().SetColor(0.3, 0.85, 0.95)
        slider_rep.GetTubeProperty().SetColor(0.22, 0.25, 0.35)
        slider_rep.GetCapProperty().SetColor(0.3, 0.35, 0.45)
        slider_rep.GetSelectedProperty().SetColor(0.5, 1.0, 1.0)
        slider_rep.GetTitleProperty().SetColor(0.75, 0.8, 0.9)
        slider_rep.GetTitleProperty().SetShadow(True)
        slider_rep.GetTitleProperty().SetShadowOffset(1, 1)
        slider_rep.GetTitleProperty().SetFontFamilyToArial()
        slider_rep.GetLabelProperty().SetColor(0.6, 0.7, 0.85)
        slider_rep.GetLabelProperty().SetShadow(True)
        slider_rep.GetLabelProperty().SetFontFamilyToArial()

        widget = vtk.vtkSliderWidget()
        widget.SetInteractor(interactor)
        widget.SetRepresentation(slider_rep)
        widget.SetAnimationModeToAnimate()
        widget.EnabledOn()

        widget.AddObserver("InteractionEvent", callback)

        return widget


# =============================================================================
# Scene Controller — Orchestrates Everything
# =============================================================================


class SceneController:
    """
    Master controller that owns all model objects and views, and orchestrates
    the update pipeline each frame.
    """

    def __init__(self):
        # Clock
        self.clock = SimulationClock()
        self.clock.speed = 0.5

        # Orbital parameters
        self.earth_orbit = OrbitParams(
            radius=EARTH_ORBIT_RADIUS,
            angular_speed=EARTH_ANGULAR_VELOCITY,
        )
        self.moon_orbit = OrbitParams(
            radius=MOON_ORBIT_RADIUS,
            angular_speed=MOON_ANGULAR_VELOCITY,
            inclination=MOON_INCLINATION_DEFAULT,
        )

        # Body states
        self.sun_state = CelestialBodyState("Sun", SUN_RADIUS)
        self.earth_state = CelestialBodyState("Earth", EARTH_RADIUS)
        self.moon_state = CelestialBodyState("Moon", MOON_RADIUS)

        # Eclipse state
        self.eclipse = EclipseState()

        # Scale exaggeration
        self.scale_factor = 1.0

        # Earth rotation speed (degrees per day of sim time)
        self.earth_spin_speed = 360.0

        # Views — per-body label colors for visibility on dark background
        self.sun_view = BodyView(
            "Sun", SUN_RADIUS, (1.0, 0.9, 0.3),
            label_color=(1.0, 0.95, 0.5), emissive=True,
        )
        self.earth_view = BodyView(
            "Earth", EARTH_RADIUS, (0.2, 0.4, 0.9),
            label_color=(0.5, 0.75, 1.0),
            glow_color=(0.28, 0.58, 1.0),
            glow_scale=1.22,
            glow_opacity=0.2,
        )
        self.moon_view = BodyView(
            "Moon", MOON_RADIUS, (0.7, 0.7, 0.7),
            label_color=(0.85, 0.85, 0.9),
            glow_color=(0.82, 0.86, 0.98),
            glow_scale=1.3,
            glow_opacity=0.16,
        )

        self.earth_orbit_view = OrbitView((0.15, 0.3, 0.7))
        self.moon_orbit_view = OrbitView((0.5, 0.5, 0.5))

        self.shadow_view = ShadowView()
        self.ray_view = SunRayView()
        self.footprint_view = EclipseFootprintView()
        self.inset_view = EclipseInsetView()
        self.status = StatusDisplay()

        # Visibility toggles
        self.show_orbits = True
        self.show_labels = True
        self.show_rays = False

    def get_all_actors(self) -> list:
        """Return all VTK actors that should be added to the renderer."""
        actors = [
            self.earth_view.glow_actor,
            self.moon_view.glow_actor,
            self.sun_view.actor,
            self.earth_view.actor,
            self.moon_view.actor,
            self.earth_orbit_view.actor,
            self.moon_orbit_view.actor,
            self.shadow_view.umbra_actor,
            self.shadow_view.penumbra_actor,
            self.ray_view.actor,
            self.footprint_view.disk_actor,
            self.footprint_view.ring_actor,
            self.sun_view.label,
            self.earth_view.label,
            self.moon_view.label,
            self.status.actor,
            self.status.title_actor,
            self.status.help_actor,
            self.status.inset_actor,
        ]
        return [actor for actor in actors if actor is not None]

    def update_scene(self) -> None:
        """
        Master update: recompute all positions, eclipse geometry, and views.

        Called every frame and whenever a slider changes.
        """
        t = self.clock.time

        # Compute body positions
        sun_pos = np.array([0.0, 0.0, 0.0])
        earth_pos = OrbitModel.compute_position(sun_pos, self.earth_orbit, t)
        moon_pos = OrbitModel.compute_position(earth_pos, self.moon_orbit, t)

        # Store states
        self.sun_state.position = sun_pos
        self.earth_state.position = earth_pos
        self.moon_state.position = moon_pos

        # Earth self-rotation
        self.earth_state.rotation_angle = (t * self.earth_spin_speed) % 360.0
        earth_tilt_axis = np.array(
            [math.sin(EARTH_AXIAL_TILT), 0.0, math.cos(EARTH_AXIAL_TILT)]
        )
        self.earth_state.rotation_axis = earth_tilt_axis

        # Compute scaled radii
        sun_r = SUN_RADIUS * self.scale_factor
        earth_r = EARTH_RADIUS * self.scale_factor
        moon_r = MOON_RADIUS * self.scale_factor

        # Eclipse computation
        self.eclipse = EclipseModel.compute(
            sun_pos, sun_r, earth_pos, earth_r, moon_pos, moon_r
        )

        # Update body views
        self.sun_view.set_position(sun_pos)
        self.sun_view.set_scale(self.scale_factor)

        self.earth_view.set_position(earth_pos)
        self.earth_view.set_scale(self.scale_factor)
        self.earth_view.set_rotation(
            self.earth_state.rotation_axis, self.earth_state.rotation_angle
        )

        self.moon_view.set_position(moon_pos, label_z_offset=1.5)
        self.moon_view.set_scale(self.scale_factor)

        # Update orbit paths
        self.earth_orbit_view.update_path(sun_pos, self.earth_orbit)
        self.moon_orbit_view.update_path(earth_pos, self.moon_orbit)

        # Update shadow cones
        self.shadow_view.update(
            moon_pos, moon_r, self.eclipse.shadow_axis, self.eclipse
        )

        # Update sun rays
        self.ray_view.update(sun_pos, sun_r, moon_pos, earth_pos)

        # Update footprint
        self.footprint_view.update(
            earth_pos, earth_r, self.eclipse.shadow_axis, self.eclipse
        )
        self.inset_view.update(sun_pos, sun_r, earth_pos, moon_pos, moon_r)

        # Update visibility
        self.earth_orbit_view.actor.SetVisibility(self.show_orbits)
        self.moon_orbit_view.actor.SetVisibility(self.show_orbits)
        self.sun_view.label.SetVisibility(self.show_labels)
        self.earth_view.label.SetVisibility(self.show_labels)
        self.moon_view.label.SetVisibility(self.show_labels)
        self.ray_view.actor.SetVisibility(self.show_rays)

        # Update HUD
        self.status.update(self.clock, self.eclipse)

    def project_all_labels(self, renderer) -> None:
        """Project all body labels from 3D world to 2D screen coordinates."""
        for view in (self.sun_view, self.earth_view, self.moon_view):
            view.project_label(renderer)


# =============================================================================
# Camera Presets
# =============================================================================


class CameraController:
    """Provides named camera presets for different viewing angles."""

    def __init__(self, renderer: vtk.vtkRenderer, scene: SceneController):
        self.renderer = renderer
        self.scene = scene

    def _apply(
        self,
        position: tuple,
        focal_point: tuple,
        view_up: tuple = (0, 0, 1),
    ) -> None:
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(*position)
        camera.SetFocalPoint(*focal_point)
        camera.SetViewUp(*view_up)
        self.renderer.ResetCameraClippingRange()

    def solar_system_view(self) -> None:
        """Wide view showing all bodies."""
        self._apply(
            position=(0, -120, 80),
            focal_point=(0, 0, 0),
        )

    def earth_view(self) -> None:
        """Close-up centered on Earth."""
        ep = self.scene.earth_state.position
        self._apply(
            position=(ep[0] - 15, ep[1] - 15, ep[2] + 12),
            focal_point=tuple(ep),
        )

    def moon_view(self) -> None:
        """Close-up centered on Moon."""
        mp = self.scene.moon_state.position
        self._apply(
            position=(mp[0] - 5, mp[1] - 5, mp[2] + 5),
            focal_point=tuple(mp),
        )

    def alignment_view(self) -> None:
        """View along the Sun-Earth line to see eclipse alignment."""
        ep = self.scene.earth_state.position
        sp = self.scene.sun_state.position
        direction = ep - sp
        dist = np.linalg.norm(direction)
        if dist > 1e-10:
            direction = direction / dist
        perp = np.cross(direction, np.array([0, 0, 1]))
        if np.linalg.norm(perp) < 1e-10:
            perp = np.array([0, 1, 0])
        perp = perp / np.linalg.norm(perp)
        cam_pos = sp - direction * 20 + perp * 5 + np.array([0, 0, 5])
        self._apply(
            position=tuple(cam_pos),
            focal_point=tuple(ep),
        )

    def shadow_view(self) -> None:
        """Close-up of the Earth-Moon system and shadow geometry."""
        ep = self.scene.earth_state.position
        mp = self.scene.moon_state.position
        midpoint = (ep + mp) / 2

        shadow_dir = self.scene.eclipse.shadow_axis
        if np.linalg.norm(shadow_dir) < 1e-10:
            shadow_dir = mp - ep
        if np.linalg.norm(shadow_dir) > 1e-10:
            shadow_dir = shadow_dir / np.linalg.norm(shadow_dir)
        else:
            shadow_dir = np.array([1.0, 0.0, 0.0])

        up = np.array([0.0, 0.0, 1.0])
        side = np.cross(shadow_dir, up)
        if np.linalg.norm(side) < 1e-10:
            side = np.array([0.0, 1.0, 0.0])
        side = side / np.linalg.norm(side)
        elevated = up + 0.45 * side
        elevated = elevated / np.linalg.norm(elevated)
        separation = max(np.linalg.norm(mp - ep), 1.0)

        self._apply(
            position=tuple(midpoint - shadow_dir * separation * 1.1 + elevated * separation * 1.45),
            focal_point=tuple(midpoint),
            view_up=tuple(up),
        )


# =============================================================================
# Main Application
# =============================================================================


def main():
    """
    Launch the interactive solar eclipse simulation.

    Sets up the VTK rendering pipeline, creates all simulation objects,
    wires slider callbacks, registers keyboard shortcuts, and starts the
    animation loop.
    """
    # -------------------------------------------------------------------------
    # Scene controller (model + views)
    # -------------------------------------------------------------------------
    scene = SceneController()

    # -------------------------------------------------------------------------
    # VTK rendering pipeline
    # -------------------------------------------------------------------------
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.02, 0.02, 0.06)
    renderer.SetBackground2(0.0, 0.0, 0.15)
    renderer.GradientBackgroundOn()
    renderer.SetUseFXAA(True)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(WINDOW_WIDTH, WINDOW_HEIGHT)
    render_window.SetWindowName("Interactive Solar Eclipse Simulation")

    inset_renderer = vtk.vtkRenderer()
    inset_renderer.SetViewport(0.72, 0.62, 0.97, 0.94)
    inset_renderer.SetBackground(0.04, 0.05, 0.1)
    inset_renderer.SetUseFXAA(True)
    render_window.AddRenderer(inset_renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)

    # -------------------------------------------------------------------------
    # Add actors to renderer
    # -------------------------------------------------------------------------
    for actor in scene.get_all_actors():
        renderer.AddActor(actor)
    scene.inset_view.add_to_renderer(inset_renderer)

    # Sun as VTK light source
    sun_light = vtk.vtkLight()
    sun_light.SetPosition(0, 0, 0)
    sun_light.SetFocalPoint(EARTH_ORBIT_RADIUS, 0, 0)
    sun_light.SetColor(1.0, 0.95, 0.8)
    sun_light.SetIntensity(1.0)
    sun_light.SetPositional(True)
    sun_light.SetConeAngle(180)
    renderer.AddLight(sun_light)

    # Dim ambient light so the unlit side of Earth is visible but dark
    ambient_light = vtk.vtkLight()
    ambient_light.SetColor(0.15, 0.15, 0.2)
    ambient_light.SetIntensity(0.3)
    ambient_light.SetPositional(False)
    renderer.AddLight(ambient_light)

    # Stars background — scatter small white points on a large sphere
    star_points = vtk.vtkPointSource()
    star_points.SetNumberOfPoints(2000)
    star_points.SetRadius(200)
    star_points.SetCenter(0, 0, 0)
    star_mapper = vtk.vtkPolyDataMapper()
    star_mapper.SetInputConnection(star_points.GetOutputPort())
    star_actor = vtk.vtkActor()
    star_actor.SetMapper(star_mapper)
    star_actor.GetProperty().SetColor(1.0, 1.0, 1.0)
    star_actor.GetProperty().SetPointSize(1.5)
    star_actor.GetProperty().SetOpacity(0.6)
    star_actor.GetProperty().SetAmbient(1.0)
    star_actor.GetProperty().SetDiffuse(0.0)
    renderer.AddActor(star_actor)

    # Ecliptic plane grid (subtle)
    grid_source = vtk.vtkPlaneSource()
    grid_source.SetOrigin(-100, -100, 0)
    grid_source.SetPoint1(100, -100, 0)
    grid_source.SetPoint2(-100, 100, 0)
    grid_source.SetXResolution(20)
    grid_source.SetYResolution(20)
    grid_mapper = vtk.vtkPolyDataMapper()
    grid_mapper.SetInputConnection(grid_source.GetOutputPort())
    grid_actor = vtk.vtkActor()
    grid_actor.SetMapper(grid_mapper)
    grid_actor.GetProperty().SetRepresentationToWireframe()
    grid_actor.GetProperty().SetColor(0.12, 0.12, 0.2)
    grid_actor.GetProperty().SetOpacity(0.2)
    grid_actor.GetProperty().SetAmbient(1.0)
    grid_actor.GetProperty().SetDiffuse(0.0)
    renderer.AddActor(grid_actor)

    # Axes widget for orientation
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.12, 0.12)
    axes_widget.EnabledOn()
    axes_widget.InteractiveOff()

    # -------------------------------------------------------------------------
    # Camera controller
    # -------------------------------------------------------------------------
    camera_ctrl = CameraController(renderer, scene)

    # -------------------------------------------------------------------------
    # Initial scene update and camera
    # -------------------------------------------------------------------------
    scene.update_scene()
    scene.status.position_for_window(WINDOW_WIDTH, WINDOW_HEIGHT)
    camera_ctrl.shadow_view()
    scene.project_all_labels(renderer)

    # -------------------------------------------------------------------------
    # Slider widgets
    # -------------------------------------------------------------------------
    # Layout: compact sliders in a left-side column to avoid overlap
    slider_x1 = 0.02
    slider_x2 = 0.18
    slider_gap = 0.072

    def make_y(index):
        base = 0.15
        return base + index * slider_gap, base + index * slider_gap

    # 1. Time slider
    def time_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.clock.set_time(val)
        scene.update_scene()
        scene.project_all_labels(renderer)
        render_window.Render()

    y1, y2 = make_y(0)
    time_slider = SliderFactory.create(
        interactor, "Time (days)", 0.0, MOON_ORBITAL_PERIOD * 2,
        0.0, slider_x1, y1, slider_x2, y2, time_cb,
    )

    # 2. Speed slider
    def speed_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.clock.set_speed(val)

    y1, y2 = make_y(1)
    speed_slider = SliderFactory.create(
        interactor, "Speed", 0.0, 5.0,
        0.5, slider_x1, y1, slider_x2, y2, speed_cb,
    )

    # 3. Scale slider
    def scale_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.scale_factor = val
        scene.update_scene()
        scene.project_all_labels(renderer)
        render_window.Render()

    y1, y2 = make_y(2)
    scale_slider = SliderFactory.create(
        interactor, "Body Scale", 0.5, 5.0,
        1.0, slider_x1, y1, slider_x2, y2, scale_cb,
    )

    # 4. Moon inclination slider
    def incl_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.moon_orbit.inclination = math.radians(val)
        scene.update_scene()
        scene.project_all_labels(renderer)
        render_window.Render()

    y1, y2 = make_y(3)
    incl_slider = SliderFactory.create(
        interactor, "Moon Incl (deg)", 0.0, 10.0,
        math.degrees(MOON_INCLINATION_DEFAULT),
        slider_x1, y1, slider_x2, y2, incl_cb,
    )

    # 5. Shadow opacity slider
    def shadow_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.shadow_view.set_opacity(val)
        render_window.Render()

    y1, y2 = make_y(4)
    shadow_slider = SliderFactory.create(
        interactor, "Shadow Alpha", 0.0, 1.0,
        0.35, slider_x1, y1, slider_x2, y2, shadow_cb,
    )

    # 6. Earth spin speed slider
    def spin_cb(obj, event):
        val = obj.GetRepresentation().GetValue()
        scene.earth_spin_speed = val

    y1, y2 = make_y(5)
    spin_slider = SliderFactory.create(
        interactor, "Earth Spin", 0.0, 1000.0,
        360.0, slider_x1, y1, slider_x2, y2, spin_cb,
    )

    # Keep references to prevent garbage collection
    _sliders = [time_slider, speed_slider, scale_slider, incl_slider,
                shadow_slider, spin_slider]

    # -------------------------------------------------------------------------
    # Keyboard shortcuts
    # -------------------------------------------------------------------------
    def on_key_press(obj, event):
        key = interactor.GetKeySym()

        if key == "1":
            camera_ctrl.solar_system_view()
        elif key == "2":
            camera_ctrl.earth_view()
        elif key == "3":
            camera_ctrl.moon_view()
        elif key == "4":
            camera_ctrl.alignment_view()
        elif key == "5":
            camera_ctrl.shadow_view()
        elif key in ("o", "O"):
            scene.show_orbits = not scene.show_orbits
            scene.update_scene()
        elif key in ("l", "L"):
            scene.show_labels = not scene.show_labels
            scene.update_scene()
        elif key in ("r", "R"):
            scene.show_rays = not scene.show_rays
            scene.update_scene()
        elif key == "space":
            scene.clock.toggle_pause()
        elif key in ("q", "Q", "Escape"):
            render_window.Finalize()
            interactor.TerminateApp()
            return

        scene.project_all_labels(renderer)
        render_window.Render()

    interactor.AddObserver("KeyPressEvent", on_key_press)

    def on_resize(obj, event):
        width, height = render_window.GetSize()
        scene.status.position_for_window(width, height)
        render_window.Render()

    interactor.AddObserver("ConfigureEvent", on_resize)

    # -------------------------------------------------------------------------
    # Animation timer
    # -------------------------------------------------------------------------
    dt_per_frame = 0.2  # base days per frame at speed 1.0

    def on_timer(obj, event):
        scene.clock.advance(dt_per_frame)

        # Sync time slider readout with current clock
        time_rep = time_slider.GetRepresentation()
        time_rep.SetValue(scene.clock.time % (MOON_ORBITAL_PERIOD * 2))

        scene.update_scene()
        scene.project_all_labels(renderer)

        # Keep sun light at Sun position
        sun_light.SetPosition(*scene.sun_state.position)
        ep = scene.earth_state.position
        sun_light.SetFocalPoint(*ep)

        render_window.Render()

    interactor.AddObserver("TimerEvent", on_timer)

    # -------------------------------------------------------------------------
    # Launch
    # -------------------------------------------------------------------------
    render_window.Render()
    interactor.Initialize()
    interactor.CreateRepeatingTimer(TIMER_INTERVAL_MS)
    interactor.Start()

    render_window.Finalize()
    interactor.TerminateApp()


if __name__ == "__main__":
    main()
