"""D2Q9 BGK lattice-Boltzmann flow past a rectangular obstacle.

The domain is open in x and periodic in y. It represents a periodic channel,
not a channel with solid top and bottom walls. Inlet velocity is imposed with
the Zou/He non-equilibrium bounce-back rule, the outlet uses population
extrapolation, and obstacle links use halfway bounce-back. The example is an
educational transient calculation; vortex shedding depends on the grid,
parameters, and simulated time and is not asserted by this solver.
"""

from __future__ import annotations

import numpy as np
import vtk

import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


# Domain and lattice-unit defaults.
MAX_ITER = 3000
# This lower-Re teaching case runs cleanly for the default duration. Larger
# configured Reynolds numbers remain useful experiments but can destabilize BGK.
REYNOLDS_NUMBER = 100.0
NX = 520
NY = 180
OBSTACLE_LENGTH = 50
OBSTACLE_HEIGHT = 30
ULB = 0.06
SNAPSHOT_INTERVAL = 100


def setup_lattice() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return D2Q9 velocities, weights, and opposite-direction indices."""
    c = np.array(
        [
            [0, 0],
            [1, 0],
            [0, 1],
            [-1, 0],
            [0, -1],
            [1, 1],
            [-1, 1],
            [-1, -1],
            [1, -1],
        ],
        dtype=np.int8,
    )
    weights = np.array(
        [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36],
        dtype=float,
    )
    opposite = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.intp)
    return c, weights, opposite


C, T, NOSLIP = setup_lattice()
Q = len(T)


def equilibrium(
    rho: np.ndarray | float,
    u: np.ndarray,
    c: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Compute second-order equilibrium for a weighted velocity lattice.

    ``u`` has shape ``(dimensions, *spatial_shape)`` and ``rho`` has the
    matching spatial shape. The sound speed squared is obtained from the
    weighted second moment of the supplied lattice, so D1Q3 and other
    isotropic lattices can use this routine as well as D2Q9.
    """
    lattice = C if c is None else np.asarray(c, dtype=float)
    lattice_weights = T if weights is None else np.asarray(weights, dtype=float)
    velocity = np.asarray(u, dtype=float)
    density = np.asarray(rho, dtype=float)

    if lattice.ndim != 2 or lattice_weights.shape != (lattice.shape[0],):
        raise ValueError("c must be (q, dimensions) and weights must have shape (q,)")
    if velocity.ndim < 1 or velocity.shape[0] != lattice.shape[1]:
        raise ValueError("u must have shape (dimensions, *spatial_shape)")
    if density.shape != velocity.shape[1:]:
        raise ValueError("rho shape must match the spatial dimensions of u")
    if not np.isfinite(lattice).all() or not np.isfinite(lattice_weights).all():
        raise ValueError("lattice velocities and weights must be finite")
    if not np.isfinite(velocity).all() or not np.isfinite(density).all():
        raise ValueError("rho and u must be finite")
    if np.any(density <= 0) or np.any(lattice_weights <= 0):
        raise ValueError("rho and lattice weights must be positive")

    weight_sum = float(np.sum(lattice_weights))
    first_moment = np.einsum("q,qd->d", lattice_weights, lattice)
    second_moment = np.einsum("q,qi,qj->ij", lattice_weights, lattice, lattice)
    cs2 = float(np.trace(second_moment) / lattice.shape[1])
    if not np.isfinite(cs2) or cs2 <= 0:
        raise ValueError("lattice must have a positive second moment")
    if not np.isclose(weight_sum, 1.0) or not np.allclose(first_moment, 0.0):
        raise ValueError("lattice weights must sum to one and have zero first moment")
    if not np.allclose(second_moment, cs2 * np.eye(lattice.shape[1]), rtol=1e-10, atol=1e-12):
        raise ValueError("lattice second moments must be isotropic")

    cu = np.einsum("qd,d...->q...", lattice, velocity)
    speed_squared = np.sum(velocity * velocity, axis=0)
    scale = lattice_weights.reshape((len(lattice_weights),) + (1,) * density.ndim)
    return scale * density * (
        1.0 + cu / cs2 + cu**2 / (2.0 * cs2**2) - speed_squared / (2.0 * cs2)
    )


def compute_macroscopic(
    fin: np.ndarray, c: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Recover density and velocity; raise on invalid populations or density."""
    lattice = C if c is None else np.asarray(c, dtype=float)
    distributions = np.asarray(fin, dtype=float)
    if lattice.ndim != 2 or distributions.ndim < 2:
        raise ValueError("c must be (q, dimensions) and fin must have spatial axes")
    if distributions.shape[0] != lattice.shape[0]:
        raise ValueError("fin direction count must match c")
    if not np.isfinite(distributions).all():
        raise FloatingPointError("distribution populations are non-finite")
    if np.any(distributions < 0):
        raise FloatingPointError("distribution populations became negative")

    rho = np.sum(distributions, axis=0)
    if not np.isfinite(rho).all() or np.any(rho <= 0):
        raise FloatingPointError("density must remain finite and positive")
    momentum = np.tensordot(lattice.T, distributions, axes=(1, 0))
    velocity = momentum / rho
    if not np.isfinite(velocity).all():
        raise FloatingPointError("velocity became non-finite")
    return rho, velocity


def _opposite_directions(c: np.ndarray) -> np.ndarray:
    """Find each lattice direction's unique opposite velocity."""
    opposite = []
    for direction in c:
        matches = np.flatnonzero(np.all(c == -direction, axis=1))
        if len(matches) != 1:
            raise ValueError("every lattice velocity must have one opposite direction")
        opposite.append(matches[0])
    return np.asarray(opposite, dtype=np.intp)


def apply_zou_he_velocity_inlet(
    fin: np.ndarray,
    inlet_velocity: np.ndarray,
    c: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Reconstruct the three incoming D2Q9 populations at the left inlet.

    The density follows from the mass constraint. Unknown populations use the
    Zou/He non-equilibrium bounce-back relations (Zou & He, 1997,
    https://arxiv.org/abs/comp-gas/9611001). ``fin`` is modified and returned.
    """
    lattice = C if c is None else np.asarray(c, dtype=float)
    lattice_weights = T if weights is None else np.asarray(weights, dtype=float)
    if lattice.shape != (9, 2) or lattice_weights.shape != (9,):
        raise ValueError("the Zou/He inlet currently requires a D2Q9 lattice")
    if not np.array_equal(lattice, C):
        raise ValueError("the Zou/He inlet expects the D2Q9 direction ordering in C")
    if not np.array_equal(lattice_weights, T):
        raise ValueError("the Zou/He inlet expects the standard D2Q9 weights")
    if fin.ndim != 3 or fin.shape[0] != 9:
        raise ValueError("fin must have shape (9, nx, ny)")
    velocity = np.asarray(inlet_velocity, dtype=float)
    if velocity.shape == (2,):
        velocity = np.repeat(velocity[:, None], fin.shape[2], axis=1)
    if velocity.shape != (2, fin.shape[2]) or not np.isfinite(velocity).all():
        raise ValueError("inlet_velocity must be finite with shape (2,) or (2, ny)")
    ux, uy = velocity
    if np.any(np.abs(ux) >= 1.0):
        raise ValueError("inlet normal velocity must have magnitude below one")

    normal_zero = lattice[:, 0] == 0
    known_outgoing = lattice[:, 0] < 0
    rho = (
        np.sum(fin[normal_zero, 0, :], axis=0)
        + 2.0 * np.sum(fin[known_outgoing, 0, :], axis=0)
    ) / (1.0 - ux)
    if not np.isfinite(rho).all() or np.any(rho <= 0):
        raise FloatingPointError("Zou/He inlet reconstruction produced invalid density")

    # C's direction order is E, N, W, S, NE, NW, SW, SE. The diagonal
    # relations include the transverse non-equilibrium correction (fS-fN),
    # which preserves prescribed uy when known N and S populations differ.
    east, north, west, south, northeast, northwest, southwest, southeast = range(1, 9)
    fin[east, 0, :] = fin[west, 0, :] + (2.0 / 3.0) * rho * ux
    fin[northeast, 0, :] = (
        fin[southwest, 0, :]
        + rho * ux / 6.0
        + rho * uy / 2.0
        + (fin[south, 0, :] - fin[north, 0, :]) / 2.0
    )
    fin[southeast, 0, :] = (
        fin[northwest, 0, :]
        + rho * ux / 6.0
        - rho * uy / 2.0
        + (fin[north, 0, :] - fin[south, 0, :]) / 2.0
    )
    return fin


def apply_open_boundaries(
    fin: np.ndarray,
    inlet_velocity: np.ndarray,
    c: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Apply a Zou/He left inlet and zero-gradient right outlet in place."""
    lattice = C if c is None else np.asarray(c, dtype=float)
    fin[lattice[:, 0] < 0, -1, :] = fin[lattice[:, 0] < 0, -2, :]
    return apply_zou_he_velocity_inlet(fin, inlet_velocity, lattice, weights)


def stream_distributions(
    fout: np.ndarray, c: np.ndarray | None = None
) -> np.ndarray:
    """Stream populations with periodic y and non-periodic x boundaries.

    The incoming populations at the two x edges are placeholders; apply the
    inlet and outlet rules after streaming before taking macroscopic moments.
    """
    lattice = C if c is None else np.asarray(c, dtype=int)
    if fout.ndim != 3 or fout.shape[0] != len(lattice):
        raise ValueError("fout must have shape (q, nx, ny) matching c")
    streamed = np.zeros_like(fout)
    for i, (cx, cy) in enumerate(lattice):
        shifted_y = np.roll(fout[i], shift=int(cy), axis=1)
        if cx > 0:
            streamed[i, 1:, :] = shifted_y[:-1, :]
        elif cx < 0:
            streamed[i, :-1, :] = shifted_y[1:, :]
        else:
            streamed[i] = shifted_y
    return streamed


def apply_link_bounce_back(
    fout: np.ndarray,
    streamed: np.ndarray,
    obstacle: np.ndarray,
    c: np.ndarray | None = None,
) -> np.ndarray:
    """Reflect post-collision populations halfway across fluid-solid links."""
    lattice = C if c is None else np.asarray(c, dtype=int)
    solid = np.asarray(obstacle, dtype=bool)
    if solid.shape != fout.shape[1:] or streamed.shape != fout.shape:
        raise ValueError("obstacle and distribution shapes do not match")
    opposite = _opposite_directions(lattice)
    for i, (cx, cy) in enumerate(lattice):
        if i == opposite[i]:
            continue
        shifted_solid = np.roll(solid, shift=-int(cy), axis=1)
        neighbor_solid = np.zeros_like(solid)
        if cx > 0:
            neighbor_solid[:-1, :] = shifted_solid[1:, :]
        elif cx < 0:
            neighbor_solid[1:, :] = shifted_solid[:-1, :]
        else:
            neighbor_solid = shifted_solid
        fluid_link = neighbor_solid & ~solid
        streamed[opposite[i], fluid_link] = fout[i, fluid_link]
    return streamed


def lbm_step(
    fin: np.ndarray,
    obstacle: np.ndarray,
    inlet_velocity: np.ndarray,
    omega: float,
    c: np.ndarray | None = None,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Advance one BGK step and return distributions, density, and velocity."""
    lattice = C if c is None else np.asarray(c, dtype=float)
    lattice_weights = T if weights is None else np.asarray(weights, dtype=float)
    if not np.isfinite(omega) or not 0.0 < omega < 2.0:
        raise ValueError("omega must be finite and between zero and two")
    if fin.shape[1] < 3:
        raise ValueError("the x domain must contain at least three nodes")

    current = apply_open_boundaries(fin.copy(), inlet_velocity, lattice, lattice_weights)
    rho, velocity = compute_macroscopic(current, lattice)
    feq = equilibrium(rho, velocity, lattice, lattice_weights)
    post_collision = current - omega * (current - feq)
    streamed = stream_distributions(post_collision, lattice)
    streamed = apply_link_bounce_back(post_collision, streamed, obstacle, lattice)
    streamed = apply_open_boundaries(streamed, inlet_velocity, lattice, lattice_weights)
    rho_next, velocity_next = compute_macroscopic(streamed, lattice)
    velocity_next[:, np.asarray(obstacle, dtype=bool)] = 0.0
    return streamed, rho_next, velocity_next


def create_obstacle_mask(
    nx: int,
    ny: int,
    obstacle_length: int,
    obstacle_height: int,
    x_center: int | None = None,
    y_center: int | None = None,
) -> np.ndarray:
    """Create an interior rectangular solid mask with shape ``(nx, ny)``."""
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) for value in (nx, ny, obstacle_length, obstacle_height)):
        raise ValueError("domain and obstacle dimensions must be integers")
    if nx < 5 or ny < 5:
        raise ValueError("the domain must be at least 5 by 5 lattice nodes")
    if obstacle_length < 1 or obstacle_height < 1:
        raise ValueError("obstacle dimensions must be positive integers")
    x_center = nx // 4 if x_center is None else x_center
    y_center = ny // 2 if y_center is None else y_center
    x0 = int(x_center) - int(obstacle_length) // 2
    y0 = int(y_center) - int(obstacle_height) // 2
    x1 = x0 + int(obstacle_length)
    y1 = y0 + int(obstacle_height)
    if x0 < 1 or x1 >= nx - 1 or y0 < 1 or y1 >= ny - 1:
        raise ValueError("obstacle must lie inside the domain and away from boundaries")
    mask = np.zeros((nx, ny), dtype=bool)
    mask[x0:x1, y0:y1] = True
    return mask


def simulate_flow(
    nx: int = NX,
    ny: int = NY,
    max_iter: int = MAX_ITER,
    reynolds: float = REYNOLDS_NUMBER,
    inlet_speed: float = ULB,
    obstacle_length: int = OBSTACLE_LENGTH,
    obstacle_height: int = OBSTACLE_HEIGHT,
    snapshot_interval: int = SNAPSHOT_INTERVAL,
    progress_callback=None,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    """Run LBM without opening a window; return post-step velocity snapshots.

    Snapshots and the final velocity have shape ``(2, ny, nx)`` and contain
    lattice-unit velocity components. The obstacle mask has shape ``(ny, nx)``
    for direct use with plotting routines.
    """
    if nx < 5 or ny < 5:
        raise ValueError("nx and ny must be at least five")
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) for value in (max_iter, snapshot_interval)):
        raise ValueError("max_iter and snapshot_interval must be integers")
    if max_iter < 1 or snapshot_interval < 1:
        raise ValueError("max_iter and snapshot_interval must be positive")
    if not np.isfinite(reynolds) or reynolds <= 0:
        raise ValueError("reynolds must be finite and positive")
    if not np.isfinite(inlet_speed) or not 0.0 < inlet_speed < 0.1:
        raise ValueError("inlet_speed must be finite, positive, and below 0.1")

    obstacle_xy = create_obstacle_mask(nx, ny, obstacle_length, obstacle_height)
    obstacle = obstacle_xy.T.copy()
    nu_lb = abs(inlet_speed) * obstacle_height / reynolds
    if nu_lb <= 0:
        raise ValueError("inlet_speed and obstacle height must give positive viscosity")
    omega = 1.0 / (3.0 * nu_lb + 0.5)
    if not 0.0 < omega < 2.0:
        raise ValueError("chosen Reynolds number gives an unstable BGK relaxation rate")

    velocity = np.zeros((2, nx, ny), dtype=float)
    velocity[0, :, :] = inlet_speed
    velocity[:, obstacle_xy] = 0.0
    fin = equilibrium(np.ones((nx, ny)), velocity)
    inlet = np.zeros((2, ny), dtype=float)
    inlet[0, :] = inlet_speed
    snapshots: list[np.ndarray] = []

    for iteration in range(1, max_iter + 1):
        fin, _, velocity = lbm_step(fin, obstacle_xy, inlet, omega)
        if iteration % snapshot_interval == 0 or iteration == max_iter:
            snapshots.append(velocity.transpose(0, 2, 1).copy())
            if progress_callback is not None:
                progress_callback(iteration, max_iter)

    return snapshots, obstacle, fin, velocity.transpose(0, 2, 1).copy()


def create_structured_grid(
    nx: int,
    ny: int,
    velocity: np.ndarray,
    obstacle_mask: np.ndarray,
) -> vtk.vtkStructuredGrid:
    """Map ``(2, ny, nx)`` velocities to point scalars, vectors, and quads."""
    if velocity.shape != (2, ny, nx) or obstacle_mask.shape != (ny, nx):
        raise ValueError("velocity or obstacle mask has the wrong shape")

    points = vtk.vtkPoints()
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(i, j, 0.0)

    speed = np.hypot(velocity[0], velocity[1]).copy()
    speed[obstacle_mask] = np.nan
    speed_array = vtk.vtkFloatArray()
    speed_array.SetName("Velocity Magnitude (lattice units)")
    vectors = vtk.vtkFloatArray()
    vectors.SetName("Velocity (lattice units)")
    vectors.SetNumberOfComponents(3)
    vectors.SetNumberOfTuples(nx * ny)
    for j in range(ny):
        for i in range(nx):
            speed_array.InsertNextValue(float(speed[j, i]))
            vx = 0.0 if obstacle_mask[j, i] else float(velocity[0, j, i])
            vy = 0.0 if obstacle_mask[j, i] else float(velocity[1, j, i])
            vectors.SetTuple3(j * nx + i, vx, vy, 0.0)

    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(nx, ny, 1)
    grid.SetPoints(points)
    grid.GetPointData().SetScalars(speed_array)
    grid.GetPointData().SetVectors(vectors)
    return grid


def vtk_velocity_plot(
    nx: int,
    ny: int,
    snapshots: list[np.ndarray],
    obstacle_mask: np.ndarray,
) -> None:
    """Animate post-step lattice velocity magnitude; solids are gray."""
    if not snapshots:
        raise ValueError("at least one velocity snapshot is required")
    finite_speeds = [
        np.hypot(snapshot[0], snapshot[1])[~obstacle_mask] for snapshot in snapshots
    ]
    global_max = max(float(np.max(speeds)) for speeds in finite_speeds)
    color_map = vtk.vtkLookupTable()
    color_map.SetRange(0.0, max(global_max, 1e-12))
    color_map.SetHueRange(0.667, 0.0)
    color_map.SetNanColor(0.25, 0.27, 0.3, 1.0)
    color_map.Build()

    renderer = vtk.vtkRenderer()
    renderer.SetUseFXAA(True)
    renderer.SetBackground(0.1, 0.1, 0.1)
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 500)
    render_window.SetWindowName("Flow Around an Obstacle - Lattice Boltzmann")
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    interactor.Initialize()

    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Velocity magnitude\n(lattice units)")
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.GetTitleTextProperty().SetFontSize(16)
    scalar_bar.GetLabelTextProperty().SetFontSize(14)
    scalar_bar.SetNumberOfLabels(5)
    renderer.AddActor2D(scalar_bar)

    surface_actor = vtk.vtkActor()
    mapper = vtk.vtkDataSetMapper()
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(0.0, max(global_max, 1e-12))
    mapper.UseLookupTableScalarRangeOn()
    mapper.InterpolateScalarsBeforeMappingOn()
    surface_actor.SetMapper(mapper)
    renderer.AddActor(surface_actor)

    status = vtk.vtkTextActor()
    status.GetTextProperty().SetFontSize(16)
    status.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    status.SetDisplayPosition(18, 18)
    renderer.AddActor2D(status)
    snapshot_index = 0

    def render_next_frame(_obj, _event):
        nonlocal snapshot_index
        snapshot = snapshots[snapshot_index]
        current_index = snapshot_index
        snapshot_index = (snapshot_index + 1) % len(snapshots)
        grid = create_structured_grid(nx, ny, snapshot, obstacle_mask)
        mapper.SetInputData(grid)
        status.SetInput(
            "Transient LBM; periodic y; open x; coordinates/velocity in lattice units\n"
            f"Snapshot {current_index + 1}/{len(snapshots)}"
        )
        if snapshot_index == 1:
            renderer.ResetCamera()
        render_window.Render()

    interactor.AddObserver("TimerEvent", render_next_frame)
    interactor.CreateRepeatingTimer(500)
    render_next_frame(None, None)
    interactor.Start()
    render_window.Finalize()
    interactor.TerminateApp()


def run_simulation() -> None:
    """Run the default transient and open its VTK animation."""
    nu_lb = abs(ULB) * OBSTACLE_HEIGHT / REYNOLDS_NUMBER
    omega = 1.0 / (3.0 * nu_lb + 0.5)
    print("D2Q9 BGK flow past a rectangular obstacle")
    print(f"Domain: {NX} x {NY} lattice nodes; y is periodic")
    print(f"Obstacle: {OBSTACLE_LENGTH} x {OBSTACLE_HEIGHT} lattice nodes")
    print(f"Reynolds number based on obstacle height: {REYNOLDS_NUMBER:g}")
    print(f"Inlet speed: {ULB:g}; nu: {nu_lb:.6g}; omega: {omega:.6f}")
    print(f"Steps: {MAX_ITER}; snapshots every {SNAPSHOT_INTERVAL} steps")

    snapshots, obstacle, _, _ = simulate_flow(
        progress_callback=lambda step, total: print(f"Step {step}/{total}")
    )
    print(f"Completed {MAX_ITER} steps; collected {len(snapshots)} post-step snapshots")
    vtk_velocity_plot(NX, NY, snapshots, obstacle)


if __name__ == "__main__":
    run_simulation()
