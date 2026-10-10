"""
Lid-Driven Cavity Flow Simulator using Artificial Compressibility Method

This module implements a 2D incompressible Navier-Stokes solver for the classic
lid-driven cavity problem, a standard benchmark in computational fluid dynamics (CFD).

Physical Problem:
-----------------
The lid-driven cavity is a square/rectangular cavity with:
- Three stationary walls (left, right, bottom) with no-slip boundary conditions
- A top wall (lid) moving at constant velocity U_lid
- The flow is driven by viscous shear from the moving lid

This creates a recirculating flow with primary and secondary vortices depending
on the Reynolds number: Re = (U_lid * L) / nu

Governing Equations:
--------------------
The incompressible Navier-Stokes equations in 2D:

1. Continuity (mass conservation):
   ∂u/∂x + ∂v/∂y = 0

2. Momentum (x-direction):
   ∂u/∂t + u(∂u/∂x) + v(∂u/∂y) = -∂p/∂x + ν(∂²u/∂x² + ∂²u/∂y²)

3. Momentum (y-direction):
   ∂v/∂t + u(∂v/∂x) + v(∂v/∂y) = -∂p/∂y + ν(∂²v/∂x² + ∂²v/∂y²)

where:
- u, v = velocity components in x and y directions (m/s)
- p = kinematic pressure p_physical/ρ (m²/s²)
- ρ = fluid density (kg/m³)
- ν = kinematic viscosity (m²/s)

Numerical Method:
-----------------
The Artificial Compressibility Method (ACM) by Chorin is used:

The continuity equation is modified to include a pressure time derivative:
   (1/β) ∂p/∂t + ∂u/∂x + ∂v/∂y = 0

where β is the artificial compressibility parameter (pseudo-speed of sound squared).

As the solution converges to steady state, ∂p/∂t → 0, recovering incompressibility.
This eliminates the need to solve a Poisson equation at each time step.

Discretization:
- Collocated node-centered velocity and pressure fields
- Second-order central differences for pressure gradient and diffusion
- First-order upwind convection
- Explicit Euler pseudo-time integration with CFL and diffusion limits

Limitations:
- Explicit Euler is conditionally stable; the time step is limited by
  convection, artificial-acoustic propagation, and viscous diffusion.
- Centered pressure and divergence operators on a collocated grid can support
  checkerboard pressure modes. This compact educational solver does not use a
  Rhie-Chow face-flux correction and its pressure field is not a benchmark
  quality pressure solution.
- Steady state is accepted only when both the discrete divergence and the
  momentum-equation residual are below tolerance.

Boundary Conditions:
- No-slip on all walls: u = v = 0
- Moving lid: u = U_lid, v = 0
- Zero pressure gradient normal to walls (Neumann)

References:
-----------
1. Chorin, A.J. (1967). "A numerical method for solving incompressible viscous
   flow problems" Journal of Computational Physics, 2(1), 12-26.
2. Ghia, U., Ghia, K.N., & Shin, C.T. (1982). "High-Re solutions for
   incompressible flow using the Navier-Stokes equations and a multigrid method"
   Journal of Computational Physics, 48(3), 387-411.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.numerics import (
    ConvergenceError,
    validate_count,
    validate_finite,
    validate_iteration_controls,
    validate_positive,
)
from src.common.field_visualization import field_points, quad_field


def initialize_grid(
    length: float, height: float, nx: int, ny: int
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Initialize the computational grid for the cavity domain.

    Creates a uniform Cartesian grid with specified number of points.
    Uses node-based (vertex-centered) values for all variables.

    Args:
        length: Length of the cavity in x-direction (m).
        height: Height of the cavity in y-direction (m).
        nx: Number of grid points in x-direction.
        ny: Number of grid points in y-direction.

    Returns:
        x: 1D array of x-coordinates.
        y: 1D array of y-coordinates.
        dx: Grid spacing in x-direction (m).
        dy: Grid spacing in y-direction (m).
    """
    validate_positive(length=length, height=height)
    validate_count("nx", nx, minimum=3)
    validate_count("ny", ny, minimum=3)
    dx = length / (nx - 1)
    dy = height / (ny - 1)
    x = np.linspace(0, length, nx)
    y = np.linspace(0, height, ny)
    return x, y, dx, dy


def initialize_fields(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Initialize velocity and pressure fields to zero.

    Args:
        nx: Number of grid points in x-direction.
        ny: Number of grid points in y-direction.

    Returns:
        u: x-velocity field, shape (ny, nx).
        v: y-velocity field, shape (ny, nx).
        p: Pressure field, shape (ny, nx).
    """
    u = np.zeros((ny, nx))
    v = np.zeros((ny, nx))
    p = np.zeros((ny, nx))
    return u, v, p


def apply_boundary_conditions(
    u: np.ndarray, v: np.ndarray, u_lid: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply boundary conditions for the lid-driven cavity.

    Boundary conditions:
    - Left wall (x=0): u=0, v=0 (no-slip)
    - Right wall (x=L): u=0, v=0 (no-slip)
    - Bottom wall (y=0): u=0, v=0 (no-slip)
    - Top wall (y=H): u=u_lid, v=0 (moving lid)

    Boundary nodes are assigned directly. The side walls own the corners, so
    the moving lid applies only between them.

    Args:
        u: x-velocity field.
        v: y-velocity field.
        u_lid: Velocity of the moving lid (m/s).

    Returns:
        u: Updated x-velocity field.
        v: Updated y-velocity field.
    """
    if u.shape != v.shape or u.ndim != 2 or min(u.shape) < 3:
        raise ValueError("u and v must be matching two-dimensional fields of at least 3x3")
    validate_finite(u_lid=u_lid)

    # The side walls own the corners; the moving lid acts between them.
    u[:, 0] = 0.0
    u[:, -1] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0

    # Bottom wall: no-slip
    u[0, :] = 0.0
    v[0, :] = 0.0

    # Top wall: moving lid (u = u_lid, v = 0), excluding wall corners.
    u[-1, :] = 0.0
    u[-1, 1:-1] = u_lid
    v[-1, :] = 0.0

    return u, v


def apply_pressure_bc(p: np.ndarray) -> np.ndarray:
    """
    Apply Neumann boundary conditions (zero pressure gradient) at walls.

    Args:
        p: Pressure field.

    Returns:
        p: Pressure field with boundary conditions applied.
    """
    p[:, 0] = p[:, 1]  # Left wall
    p[:, -1] = p[:, -2]  # Right wall
    p[0, :] = p[1, :]  # Bottom wall
    p[-1, :] = p[-2, :]  # Top wall
    return p


def step_acm(
    u: np.ndarray,
    v: np.ndarray,
    p: np.ndarray,
    dx: float,
    dy: float,
    dt: float,
    nu: float,
    beta: float,
    u_lid: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform one time step using the Artificial Compressibility Method (ACM).

    The ACM modifies the continuity equation to include a pressure time derivative:
       (1/β) ∂p/∂t + ∂u/∂x + ∂v/∂y = 0

    Combined with the momentum equations:
       ∂u/∂t + u·∇u = -(1/ρ)∂p/∂x + ν∇²u
       ∂v/∂t + u·∇v = -(1/ρ)∂p/∂y + ν∇²v

    As the solution marches to steady state, ∂p/∂t → 0, recovering incompressibility.

    Uses first-order upwind convection and centered pressure/diffusion terms.
    The explicit update is conditionally stable; choose dt with
    :func:`compute_stable_timestep`.

    Args:
        u: x-velocity field.
        v: y-velocity field.
        p: Pressure field.
        dx: Grid spacing in x (m).
        dy: Grid spacing in y (m).
        dt: Time step (s).
        nu: Kinematic viscosity (m²/s).
        beta: Artificial compressibility parameter.
        u_lid: Lid velocity (m/s).

    Returns:
        u_new: Updated x-velocity field.
        v_new: Updated y-velocity field.
        p_new: Updated pressure field.
    """
    if u.shape != v.shape or u.shape != p.shape or u.ndim != 2 or min(u.shape) < 3:
        raise ValueError("u, v, and p must be matching two-dimensional fields of at least 3x3")
    validate_positive(dx=dx, dy=dy, dt=dt, nu=nu, beta=beta)
    validate_finite(u_lid=u_lid)

    rhs_u, rhs_v = momentum_rhs(u, v, p, dx, dy, nu)
    u_new = u.copy()
    v_new = v.copy()
    p_new = p.copy()
    u_new[1:-1, 1:-1] += dt * rhs_u
    v_new[1:-1, 1:-1] += dt * rhs_v

    div = (
        (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)
        + (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy)
    )
    p_new[1:-1, 1:-1] -= dt * beta * div

    # Apply boundary conditions
    u_new, v_new = apply_boundary_conditions(u_new, v_new, u_lid)
    p_new = apply_pressure_bc(p_new)

    return u_new, v_new, p_new


def momentum_rhs(
    u: np.ndarray, v: np.ndarray, p: np.ndarray, dx: float, dy: float, nu: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the explicit momentum-equation right-hand sides on the interior."""
    ui = u[1:-1, 1:-1]
    vi = v[1:-1, 1:-1]
    up = u[1:-1, 2:]
    um = u[1:-1, :-2]
    un = u[2:, 1:-1]
    us = u[:-2, 1:-1]
    vp = v[1:-1, 2:]
    vm = v[1:-1, :-2]
    vn = v[2:, 1:-1]
    vs = v[:-2, 1:-1]

    dudx = np.where(ui >= 0, (ui - um) / dx, (up - ui) / dx)
    dudy = np.where(vi >= 0, (ui - us) / dy, (un - ui) / dy)
    dvdx = np.where(ui >= 0, (vi - vm) / dx, (vp - vi) / dx)
    dvdy = np.where(vi >= 0, (vi - vs) / dy, (vn - vi) / dy)

    diff_u = nu * ((up - 2 * ui + um) / dx**2 + (un - 2 * ui + us) / dy**2)
    diff_v = nu * ((vp - 2 * vi + vm) / dx**2 + (vn - 2 * vi + vs) / dy**2)
    dpdx = (p[1:-1, 2:] - p[1:-1, :-2]) / (2 * dx)
    dpdy = (p[2:, 1:-1] - p[:-2, 1:-1]) / (2 * dy)
    rhs_u = -ui * dudx - vi * dudy - dpdx + diff_u
    rhs_v = -ui * dvdx - vi * dvdy - dpdy + diff_v
    return rhs_u, rhs_v


def compute_momentum_residual(
    u: np.ndarray, v: np.ndarray, p: np.ndarray, dx: float, dy: float, nu: float
) -> float:
    """Return the maximum absolute steady momentum residual in the interior."""
    rhs_u, rhs_v = momentum_rhs(u, v, p, dx, dy, nu)
    return float(max(np.max(np.abs(rhs_u)), np.max(np.abs(rhs_v))))


def compute_divergence(u: np.ndarray, v: np.ndarray, dx: float, dy: float) -> float:
    """
    Compute the maximum divergence of the velocity field.

    For an incompressible flow, ∇·u = ∂u/∂x + ∂v/∂y should be zero.
    This is used to verify the pressure correction is working.

    Args:
        u: x-velocity field.
        v: y-velocity field.
        dx: Grid spacing in x-direction (m).
        dy: Grid spacing in y-direction (m).

    Returns:
        max_div: Maximum absolute divergence in the domain.
    """
    if u.shape != v.shape or u.ndim != 2 or min(u.shape) < 3:
        raise ValueError("u and v must be matching two-dimensional fields of at least 3x3")
    validate_positive(dx=dx, dy=dy)
    if not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError("velocity fields must contain only finite values")
    # Compute divergence using the same centered operator used by pressure update.
    dudx = (u[1:-1, 2:] - u[1:-1, :-2]) / (2 * dx)
    dvdy = (v[2:, 1:-1] - v[:-2, 1:-1]) / (2 * dy)
    div = np.abs(dudx + dvdy)
    return float(np.max(div))


def compute_stable_timestep(
    dx: float,
    dy: float,
    nu: float,
    beta: float,
    u: np.ndarray,
    v: np.ndarray,
    safety: float = 0.25,
) -> float:
    """Choose a conservative explicit step from current flow and grid scales.

    The convective/acoustic bound uses actual maximum component speeds and the
    diffusion bound uses both grid spacings. These estimates are safeguards,
    not a claim of unconditional stability.
    """
    validate_positive(dx=dx, dy=dy, nu=nu, beta=beta, safety=safety)
    if safety >= 1:
        raise ValueError("safety must be less than one")
    if u.shape != v.shape or not np.isfinite(u).all() or not np.isfinite(v).all():
        raise ValueError("u and v must be matching finite velocity fields")

    umax = float(np.max(np.abs(u)))
    vmax = float(np.max(np.abs(v)))
    acoustic_rate = np.sqrt(beta) * np.hypot(1.0 / dx, 1.0 / dy)
    convective_rate = umax / dx + vmax / dy
    dt_wave = safety / max(acoustic_rate + convective_rate, np.finfo(float).tiny)
    diffusion_rate = 2.0 * nu * (1.0 / dx**2 + 1.0 / dy**2)
    dt_diffusion = safety / diffusion_rate
    return float(min(dt_wave, dt_diffusion))


def create_points(
    x: np.ndarray, y: np.ndarray, velocity_mag: np.ndarray
) -> tuple[vtk.vtkPoints, vtk.vtkFloatArray]:
    """
    Create VTK points and velocity values for visualization.

    Args:
        x: 1D array of x-coordinates.
        y: 1D array of y-coordinates.
        velocity_mag: 2D array of velocity magnitude values.

    Returns:
        points: VTK points object.
        velocity_values: VTK float array with velocity magnitude.
    """
    return field_points(x, y, velocity_mag, name="Velocity Magnitude")


def create_poly_data(
    nx: int, ny: int, points: vtk.vtkPoints, velocity_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """
    Create VTK poly data for surface visualization.

    Args:
        nx: Number of grid points in x-direction.
        ny: Number of grid points in y-direction.
        points: VTK points object.
        velocity_values: VTK float array with scalar values.

    Returns:
        polyData: VTK poly data object for rendering.
    """
    if points.GetNumberOfPoints() != nx * ny or velocity_values.GetNumberOfTuples() != nx * ny:
        raise ValueError("Point and scalar counts must match the grid dimensions")
    return quad_field(nx, ny, points, velocity_values)


def create_color_map(
    data: np.ndarray, hue_range: tuple[float, float] = (0.667, 0)
) -> vtk.vtkLookupTable:
    """
    Create a color lookup table for VTK visualization.

    Args:
        data: Array of values to map colors to.
        hue_range: Range of hue values (default blue to red).

    Returns:
        color_map: VTK lookup table.
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(data), np.max(data)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(*hue_range)
    color_map.Build()
    return color_map


def vtk_velocity_plot(
    x: np.ndarray, y: np.ndarray, u: np.ndarray, v: np.ndarray, p: np.ndarray
) -> None:
    """
    Visualize the final velocity and pressure fields using VTK.

    Displays:
    - Velocity magnitude as colored surface
    - Scalar bar legend

    Args:
        x: 1D array of x-coordinates.
        y: 1D array of y-coordinates.
        u: x-velocity field.
        v: y-velocity field.
        p: Pressure field.
    """
    # Compute velocity magnitude
    velocity_mag = np.sqrt(u**2 + v**2)

    # Create VTK objects
    points, velocity_values = create_points(x, y, velocity_mag)
    polyData = create_poly_data(len(x), len(y), points, velocity_values)
    color_map = create_color_map(velocity_mag)

    # Mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(velocity_mag), np.max(velocity_mag))
    mapper.UseLookupTableScalarRangeOn()
    mapper.InterpolateScalarsBeforeMappingOn()

    # Surface actor
    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Velocity (m/s)")
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.GetTitleTextProperty().SetFontSize(16)
    scalar_bar.GetLabelTextProperty().SetFontSize(14)
    scalar_bar.SetNumberOfLabels(5)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(scalar_bar)
    renderer.SetBackground(0.1, 0.1, 0.1)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)
    render_window.SetWindowName("Lid-Driven Cavity Flow - Velocity Magnitude")

    # Interactor
    render_interactor = vtk.vtkRenderWindowInteractor()
    render_interactor.SetRenderWindow(render_window)
    render_interactor.Initialize()

    # Interaction style
    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    render_interactor.SetInteractorStyle(interact_style)

    # Axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(render_interactor)
    axes_widget.SetViewport(0, 0, 0.2, 0.2)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Reset camera and render
    renderer.ResetCamera()
    render_window.Render()
    render_interactor.Start()

    # Cleanup
    render_window.Finalize()
    render_interactor.TerminateApp()


def run_simulation(
    length: float = 1.0,
    height: float = 1.0,
    nx: int = 41,
    ny: int = 41,
    reynolds: float = 100.0,
    u_lid: float = 1.0,
    max_iter: int = 100000,
    print_interval: int = 500,
    nu: float | None = None,
    tolerance: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Run the lid-driven cavity flow simulation using Artificial Compressibility.

    Simulates incompressible viscous flow in a cavity driven by a moving lid
    using the Artificial Compressibility Method. It returns only when both
    divergence and steady momentum residual meet ``tolerance``; exhausted
    iteration limits raise :class:`ConvergenceError`.

    Args:
        length: Cavity length in x-direction (m).
        height: Cavity height in y-direction (m).
        nx: Number of grid points in x-direction.
        ny: Number of grid points in y-direction.
        reynolds: Reynolds number Re = abs(U_lid) * L / nu, used if ``nu`` is omitted.
        u_lid: Lid velocity (m/s).
        max_iter: Maximum number of iterations.
        print_interval: Iteration interval for printing progress.
        nu: Optional kinematic viscosity. Required when the lid is stationary.
        tolerance: Absolute thresholds for divergence and momentum residual.

    Returns:
        x: 1D array of x-coordinates.
        y: 1D array of y-coordinates.
        u: Final x-velocity field.
        v: Final y-velocity field.
        p: Final pressure field.
    """
    validate_positive(length=length, height=height, reynolds=reynolds)
    nx = validate_count("nx", nx, minimum=3)
    ny = validate_count("ny", ny, minimum=3)
    tolerance, max_iter = validate_iteration_controls(tolerance, max_iter)
    validate_count("print_interval", print_interval)
    validate_finite(u_lid=u_lid)
    if nu is None:
        if u_lid == 0:
            raise ValueError("nu must be supplied when u_lid is zero")
        nu = abs(u_lid) * length / reynolds
    else:
        validate_positive(nu=nu)

    # Retain a nonzero pseudo-acoustic speed for the stationary-lid case.
    beta = max(1.0, abs(u_lid) ** 2)

    print("=" * 60)
    print("Lid-Driven Cavity Flow Simulation")
    print("Artificial Compressibility Method")
    print("=" * 60)
    print(f"Domain: {length} x {height} m")
    print(f"Grid: {nx} x {ny} points")
    print(f"Reynolds number: {reynolds}")
    print(f"Lid velocity: {u_lid} m/s")
    print(f"Kinematic viscosity: {nu:.6f} m²/s")
    print(f"Artificial compressibility (β): {beta}")
    print("=" * 60)

    # Initialize grid and fields
    x, y, dx, dy = initialize_grid(length, height, nx, ny)
    u, v, p = initialize_fields(nx, ny)

    # Apply initial boundary conditions
    u, v = apply_boundary_conditions(u, v, u_lid)
    p = apply_pressure_bc(p)

    # Compute stable time step
    dt = compute_stable_timestep(dx, dy, nu, beta, u, v)
    print(f"Initial time step: {dt:.6g} s")

    for iteration in range(1, max_iter + 1):
        dt = compute_stable_timestep(dx, dy, nu, beta, u, v)
        # Perform one ACM time step
        u, v, p = step_acm(u, v, p, dx, dy, dt, nu, beta, u_lid)

        divergence = compute_divergence(u, v, dx, dy)
        momentum_residual = compute_momentum_residual(u, v, p, dx, dy, nu)
        residual = max(divergence, momentum_residual)

        if iteration % print_interval == 0 or iteration == 1:
            u_max = np.max(np.abs(u[1:-1, 1:-1]))
            v_max = np.max(np.abs(v[1:-1, 1:-1]))
            print(
                f"Iter {iteration:5d} | div: {divergence:.2e} | "
                f"momentum: {momentum_residual:.2e} | "
                f"u_max: {u_max:.4f} | v_max: {v_max:.4f}"
            )

        if divergence <= tolerance and momentum_residual <= tolerance:
            print(f"\nSteady state reached at iteration {iteration}")
            break
    else:
        raise ConvergenceError(
            "Cavity solver did not satisfy divergence and momentum tolerances",
            max_iter,
            residual,
        )

    # Final statistics
    div_max = compute_divergence(u, v, dx, dy)
    print("\n" + "=" * 60)
    print("Simulation Complete")
    print("=" * 60)
    print(f"Final maximum divergence: {div_max:.2e}")
    print(
        "Final maximum momentum residual: "
        f"{compute_momentum_residual(u, v, p, dx, dy, nu):.2e}"
    )
    print(f"Maximum u-velocity (interior): {np.max(u[1:-1, 1:-1]):.6f} m/s")
    print(f"Minimum u-velocity (interior): {np.min(u[1:-1, 1:-1]):.6f} m/s")
    print(f"Maximum v-velocity (interior): {np.max(v[1:-1, 1:-1]):.6f} m/s")
    print(f"Minimum v-velocity (interior): {np.min(v[1:-1, 1:-1]):.6f} m/s")
    print("=" * 60)

    return x, y, u, v, p


def plot_results(
    x: np.ndarray, y: np.ndarray, u: np.ndarray, v: np.ndarray, p: np.ndarray
) -> None:
    """
    Plot simulation results using matplotlib.

    Creates a 2x2 figure showing:
    - Velocity magnitude contour
    - Streamlines
    - Pressure contour
    - Velocity profiles at cavity center

    Args:
        x: 1D array of x-coordinates.
        y: 1D array of y-coordinates.
        u: x-velocity field.
        v: y-velocity field.
        p: Pressure field.
    """
    X, Y = np.meshgrid(x, y)
    velocity_mag = np.sqrt(u**2 + v**2)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("Lid-Driven Cavity Flow Results", fontsize=14)

    # Velocity magnitude
    ax1 = axes[0, 0]
    c1 = ax1.contourf(X, Y, velocity_mag, levels=20, cmap="jet")
    plt.colorbar(c1, ax=ax1, label="Velocity (m/s)")
    ax1.set_xlabel("x (m)")
    ax1.set_ylabel("y (m)")
    ax1.set_title("Velocity Magnitude")
    ax1.set_aspect("equal")

    # Streamlines
    ax2 = axes[0, 1]
    ax2.streamplot(X, Y, u, v, color=velocity_mag, cmap="jet", density=1.5)
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel("y (m)")
    ax2.set_title("Streamlines")
    ax2.set_aspect("equal")

    # Pressure
    ax3 = axes[1, 0]
    c3 = ax3.contourf(X, Y, p, levels=20, cmap="coolwarm")
    plt.colorbar(c3, ax=ax3, label="Kinematic pressure (m²/s²)")
    ax3.set_xlabel("x (m)")
    ax3.set_ylabel("y (m)")
    ax3.set_title("Kinematic Pressure Field")
    ax3.set_aspect("equal")

    # Velocity profiles at center
    ax4 = axes[1, 1]
    mid_x = len(x) // 2
    mid_y = len(y) // 2
    ax4.plot(u[:, mid_x], y, "b-", label=f"u at x={x[mid_x]:.2f}")
    ax4.plot(v[mid_y, :], x, "r--", label=f"v at y={y[mid_y]:.2f}")
    ax4.set_xlabel("Velocity (m/s)")
    ax4.set_ylabel("Position (m)")
    ax4.set_title("Velocity Profiles at Cavity Center")
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


def main():
    """
    Main function to run the lid-driven cavity simulation.

    Configures simulation parameters, runs the solver, and displays
    results using both matplotlib and VTK visualization.
    """
    # Simulation parameters
    # Adjust Reynolds number for different flow regimes:
    # Re = 100: Steady, single primary vortex
    # Re = 400: Steady, with secondary corner vortices
    # Re = 1000+: May require finer grid for accuracy

    x, y, u, v, p = run_simulation(
        length=1.0,  # Cavity length (m)
        height=1.0,  # Cavity height (m)
        nx=41,  # Grid points in x
        ny=41,  # Grid points in y
        reynolds=100.0,  # Reynolds number
        u_lid=1.0,  # Lid velocity (m/s)
        max_iter=100000,  # Maximum iterations
        print_interval=500,  # Print every N iterations
    )

    # Matplotlib visualization
    plot_results(x, y, u, v, p)

    # VTK visualization
    vtk_velocity_plot(x, y, u, v, p)


if __name__ == "__main__":
    main()
