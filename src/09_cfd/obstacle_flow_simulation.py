"""
Flow Around an Obstacle using the Lattice Boltzmann Method (LBM)

This module implements a 2D incompressible flow simulation around a rectangular
obstacle using the Lattice Boltzmann Method with the BGK collision operator.
This is a classic CFD benchmark problem demonstrating vortex shedding (von Kármán
vortex street) at moderate Reynolds numbers.

Physical Problem:
-----------------
Flow past a bluff body (rectangular obstacle) in a channel:
- Uniform inlet flow from the left boundary
- Outflow condition at the right boundary
- The obstacle creates flow separation and recirculation zones
- At Re > 40, periodic vortex shedding occurs (von Kármán vortex street)

The Reynolds number is defined as: Re = U * D / ν
where:
- U = characteristic velocity (inlet velocity)
- D = characteristic length (obstacle height)
- ν = kinematic viscosity

The Lattice Boltzmann Method:
-----------------------------
LBM is a mesoscopic approach that models fluid dynamics by tracking the evolution
of particle distribution functions f_i on a discrete lattice.

The fundamental equation is the discrete Boltzmann equation with BGK collision:
   f_i(x + c_i*Δt, t + Δt) = f_i(x, t) - (1/τ)[f_i(x, t) - f_i^eq(x, t)]

where:
- f_i = distribution function for velocity direction i
- c_i = discrete velocity vector
- τ = relaxation time (related to viscosity)
- f_i^eq = equilibrium distribution (Maxwell-Boltzmann)

D2Q9 Lattice:
-------------
The D2Q9 lattice (2D with 9 velocities) is used:

        6   2   5
          \\ | /
        3 - 0 - 1
          / | \\
        7   4   8

Velocity vectors c_i:
- c_0 = (0, 0)     : rest particles, weight w_0 = 4/9
- c_1,2,3,4        : axis-aligned, weight w_i = 1/9
- c_5,6,7,8        : diagonal, weight w_i = 1/36

Equilibrium Distribution:
-------------------------
f_i^eq = w_i * ρ * [1 + (c_i · u)/c_s² + (c_i · u)²/(2*c_s⁴) - u²/(2*c_s²)]

where c_s = 1/√3 is the lattice speed of sound.

Macroscopic Variables:
----------------------
The macroscopic density and velocity are recovered as moments:
- ρ = Σ f_i              (density)
- ρu = Σ f_i * c_i       (momentum)

The kinematic viscosity is related to relaxation time:
- ν = c_s² * (τ - 0.5) * Δt = (1/3) * (τ - 0.5)  [in lattice units]

Boundary Conditions:
--------------------
1. Inlet (left): Zou/He velocity boundary condition
   - Prescribed velocity profile
   - Non-equilibrium bounce-back for unknown distributions

2. Outlet (right): Zero-gradient (extrapolation)
   - f_i(-1) = f_i(-2) for outgoing populations

3. Obstacle: Half-way bounce-back (no-slip)
   - f_i = f_(-i) at obstacle nodes
   - Creates second-order accurate no-slip walls

Algorithm Steps:
----------------
1. Compute macroscopic variables (ρ, u) from distribution functions
2. Apply boundary conditions
3. Compute equilibrium distribution f^eq
4. Collision: f_out = f_in - (1/τ)(f_in - f^eq)
5. Apply bounce-back at obstacles
6. Streaming: f_i(x + c_i, t+1) = f_out_i(x, t)

References:
-----------
1. Chen, S., & Doolen, G.D. (1998). "Lattice Boltzmann Method for Fluid Flows"
   Annual Review of Fluid Mechanics, 30, 329-364.
2. Zou, Q., & He, X. (1997). "On pressure and velocity boundary conditions
   for the lattice Boltzmann BGK model"
   Physics of Fluids, 9(6), 1591-1598.
3. Succi, S. (2001). "The Lattice Boltzmann Equation for Fluid Dynamics
   and Beyond" Oxford University Press.
"""

import numpy as np
import vtk


# =============================================================================
# Simulation Configuration
# =============================================================================

# Domain and time parameters
MAX_ITER: int = 3000  # Total number of time iterations
REYNOLDS_NUMBER: float = 220.0  # Reynolds number (>40 for vortex shedding)

# Grid dimensions (in lattice units)
NX: int = 520  # Length of domain
NY: int = 180  # Height of domain
LY: float = NY - 1.0  # Reference length for velocity profile

# D2Q9 lattice parameters
Q: int = 9  # Number of discrete velocities

# Lattice velocity and derived quantities
ULB: float = 0.04  # Characteristic velocity in lattice units (keep < 0.1 for stability)
NU_LB: float = ULB * (NY / 9) / REYNOLDS_NUMBER  # Kinematic viscosity in lattice units
OMEGA: float = 1.0 / (3.0 * NU_LB + 0.5)  # BGK relaxation parameter (1/τ)


# =============================================================================
# D2Q9 Lattice Definition
# =============================================================================


def setup_lattice() -> (
    tuple[np.ndarray, np.ndarray, list[int], np.ndarray, np.ndarray, np.ndarray]
):
    """
    Set up the D2Q9 lattice velocities and weights.

    The D2Q9 lattice has 9 discrete velocities arranged as:
        6   2   5
          \\ | /
        3 - 0 - 1
          / | \\
        7   4   8

    Returns:
        c: Lattice velocity vectors, shape (9, 2).
        t: Lattice weights, shape (9,).
        noslip: Indices for bounce-back (opposite direction).
        i1: Indices for velocities pointing left (for right boundary).
        i2: Indices for velocities pointing vertically (for boundary calc).
        i3: Indices for velocities pointing right (for left boundary).
    """
    # Lattice velocities: c[i] = (cx, cy)
    # Ordering: 0=rest, 1=E, 2=N, 3=W, 4=S, 5=NE, 6=NW, 7=SW, 8=SE
    c = np.array(
        [
            [0, 0],  # 0: rest
            [1, 0],  # 1: East
            [0, 1],  # 2: North
            [-1, 0],  # 3: West
            [0, -1],  # 4: South
            [1, 1],  # 5: NE
            [-1, 1],  # 6: NW
            [-1, -1],  # 7: SW
            [1, -1],  # 8: SE
        ]
    )

    # Lattice weights for D2Q9
    t = np.array(
        [
            4.0 / 9.0,  # rest: w_0 = 4/9
            1.0 / 9.0,  # E
            1.0 / 9.0,  # N
            1.0 / 9.0,  # W
            1.0 / 9.0,  # S
            1.0 / 36.0,  # NE
            1.0 / 36.0,  # NW
            1.0 / 36.0,  # SW
            1.0 / 36.0,  # SE
        ]
    )

    # Bounce-back indices: noslip[i] gives the opposite direction
    noslip = [0, 3, 4, 1, 2, 7, 8, 5, 6]

    # Velocity direction groups for boundary conditions
    i1 = np.array([3, 6, 7])  # Velocities pointing left (cx < 0)
    i2 = np.array([0, 2, 4])  # Velocities with cx = 0
    i3 = np.array([1, 5, 8])  # Velocities pointing right (cx > 0)

    return c, t, noslip, i1, i2, i3


# Initialize lattice constants
C, T, NOSLIP, I1, I2, I3 = setup_lattice()


# =============================================================================
# LBM Core Functions
# =============================================================================


def equilibrium(rho: np.ndarray, u: np.ndarray) -> np.ndarray:
    """
    Compute the equilibrium distribution function.

    The equilibrium distribution for incompressible flow is:
    f_i^eq = w_i * ρ * [1 + (c_i · u)/c_s² + (c_i · u)²/(2*c_s⁴) - u²/(2*c_s²)]

    For D2Q9 with c_s² = 1/3:
    f_i^eq = w_i * ρ * [1 + 3*(c_i · u) + 4.5*(c_i · u)² - 1.5*u²]

    Args:
        rho: Density field, shape (NX, NY).
        u: Velocity field, shape (2, NX, NY) where u[0]=u_x, u[1]=u_y.

    Returns:
        feq: Equilibrium distribution, shape (Q, NX, NY).
    """
    # Compute c_i · u for all velocities: shape (Q, NX, NY)
    cu = 3.0 * np.dot(C, u.transpose(1, 0, 2))

    # Compute u² = u_x² + u_y²: shape (NX, NY)
    usqr = 1.5 * (u[0] ** 2 + u[1] ** 2)

    # Compute equilibrium for each velocity direction
    feq = np.zeros((Q, NX, NY))
    for i in range(Q):
        feq[i] = rho * T[i] * (1.0 + cu[i] + 0.5 * cu[i] ** 2 - usqr)

    return feq


def compute_macroscopic(fin: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute macroscopic density and velocity from distribution functions.

    The macroscopic variables are moments of the distribution:
    - ρ = Σ_i f_i
    - ρu = Σ_i f_i * c_i

    Args:
        fin: Distribution function, shape (Q, NX, NY).

    Returns:
        rho: Density field, shape (NX, NY).
        u: Velocity field, shape (2, NX, NY).
    """
    rho = np.sum(fin, axis=0)
    # Add small epsilon to prevent division by zero in rare cases
    rho_safe = np.maximum(rho, 1e-10)
    u = np.dot(C.T, fin.transpose(1, 0, 2)) / rho_safe
    return rho, u


# =============================================================================
# Visualization Functions
# =============================================================================


def create_points(
    nx: int, ny: int, u: np.ndarray
) -> tuple[vtk.vtkPoints, vtk.vtkFloatArray]:
    """
    Create VTK points and velocity values for visualization.

    Args:
        nx: Number of points in x-direction.
        ny: Number of points in y-direction.
        u: Velocity magnitude field, shape (ny, nx).

    Returns:
        points: VTK points object.
        velocity_values: VTK float array with velocity magnitudes.
    """
    points = vtk.vtkPoints()
    velocity_values = vtk.vtkFloatArray()
    velocity_values.SetName("Velocity Magnitude")

    assert (
        u.shape[0] == ny and u.shape[1] == nx
    ), f"Dimension mismatch: expected ({ny}, {nx}), got {u.shape}"

    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(i, j, 0)
            velocity_values.InsertNextValue(u[j, i])

    return points, velocity_values


def create_poly_data(
    nx: int, ny: int, points: vtk.vtkPoints, velocity_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """
    Create VTK poly data for surface visualization.

    Args:
        nx: Number of points in x-direction.
        ny: Number of points in y-direction.
        points: VTK points object.
        velocity_values: VTK float array with scalar values.

    Returns:
        polyData: VTK poly data object for rendering.
    """
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    cells = vtk.vtkCellArray()
    for j in range(ny - 1):
        for i in range(nx - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, j * nx + i)
            quad.GetPointIds().SetId(1, (j + 1) * nx + i)
            quad.GetPointIds().SetId(2, (j + 1) * nx + (i + 1))
            quad.GetPointIds().SetId(3, j * nx + (i + 1))
            cells.InsertNextCell(quad)

    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(velocity_values)

    return polyData


def create_color_map(
    u: np.ndarray, hue_range: tuple[float, float] = (0.667, 0.0)
) -> vtk.vtkLookupTable:
    """
    Create a color lookup table for VTK visualization.

    Args:
        u: Array of values to map colors to.
        hue_range: Range of hue values (default blue to red).

    Returns:
        color_map: VTK lookup table.
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(u), np.max(u)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(*hue_range)
    color_map.Build()
    return color_map


def vtk_velocity_plot(nx: int, ny: int, snapshots: list[np.ndarray]) -> None:
    """
    Animate the velocity field evolution using VTK.

    Creates an interactive visualization with:
    - Velocity magnitude as colored surface
    - Scalar bar legend
    - Timer-based animation through snapshots

    Args:
        nx: Number of points in x-direction.
        ny: Number of points in y-direction.
        snapshots: List of velocity magnitude arrays, each shape (ny, nx).
    """
    if not snapshots:
        print("No snapshots to visualize.")
        return

    # Calculate global min/max for consistent color scaling
    global_min_u = min(np.min(s) for s in snapshots)
    global_max_u = max(np.max(s) for s in snapshots)

    # Setup renderer
    renderer = vtk.vtkRenderer()
    renderer.SetUseFXAA(True)
    renderer.SetBackground(0.1, 0.1, 0.1)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 400)
    render_window.SetWindowName("Flow Around Obstacle - Lattice Boltzmann Method")

    # Interactor
    render_interactor = vtk.vtkRenderWindowInteractor()
    render_interactor.SetRenderWindow(render_window)
    render_interactor.Initialize()

    # State for animation
    surface_actor = None
    scalar_bar = None
    snapshot_idx = [0]  # Use list for mutable closure

    def update_visualization(obj, event):
        nonlocal surface_actor, scalar_bar

        # Remove current actors
        if surface_actor:
            renderer.RemoveActor(surface_actor)
        if scalar_bar:
            renderer.RemoveActor(scalar_bar)

        # Check if more snapshots available
        if snapshot_idx[0] >= len(snapshots):
            snapshot_idx[0] = 0  # Loop animation

        # Get current snapshot
        u_np = snapshots[snapshot_idx[0]]
        snapshot_idx[0] += 1

        # Create visualization
        points, velocity_values = create_points(nx, ny, u_np)
        polyData = create_poly_data(nx, ny, points, velocity_values)
        color_map = create_color_map(u_np)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)
        mapper.SetLookupTable(color_map)
        mapper.SetScalarModeToUsePointData()
        mapper.SetScalarRange(global_min_u, global_max_u)

        surface_actor = vtk.vtkActor()
        surface_actor.SetMapper(mapper)

        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(color_map)
        scalar_bar.SetTitle("Velocity (lattice units)")
        scalar_bar.SetNumberOfLabels(5)

        renderer.AddActor(surface_actor)
        renderer.AddActor(scalar_bar)
        renderer.ResetCamera()
        render_window.Render()

    # Add timer for animation
    render_interactor.AddObserver("TimerEvent", update_visualization)
    render_interactor.CreateRepeatingTimer(500)  # Update every 500ms

    # Initial render
    update_visualization(None, None)
    render_window.Render()
    render_interactor.Start()

    # Cleanup
    render_window.Finalize()
    render_interactor.TerminateApp()


# =============================================================================
# Main Simulation
# =============================================================================


def run_simulation() -> None:
    """
    Run the Lattice Boltzmann simulation for flow around an obstacle.

    This function:
    1. Initializes the rectangular obstacle and velocity field
    2. Runs the LBM time-stepping loop
    3. Collects snapshots of velocity magnitude
    4. Launches VTK visualization

    The simulation demonstrates:
    - Vortex shedding behind the obstacle (von Kármán vortex street)
    - Flow separation and recirculation zones
    - Time-dependent flow at Re = 220
    """
    print("=" * 70)
    print("Lattice Boltzmann Simulation: Flow Around an Obstacle")
    print("=" * 70)
    print(f"Domain size: {NX} x {NY} lattice units")
    print(f"Reynolds number: {REYNOLDS_NUMBER}")
    print(f"Relaxation parameter (omega): {OMEGA:.4f}")
    print(f"Kinematic viscosity (lattice units): {NU_LB:.6f}")
    print(f"Maximum iterations: {MAX_ITER}")
    print("=" * 70)

    # Define rectangular obstacle (bluff body)
    obstacle_length = 50  # Length in x-direction
    obstacle_height = 10  # Height in y-direction
    obstacle_x_center = NX // 4  # Position at 1/4 of domain
    obstacle_y_center = NY // 2  # Centered vertically

    # Create obstacle mask
    obstacle = np.fromfunction(
        lambda x, y: (
            (np.abs(x - obstacle_x_center) < obstacle_length / 2)
            & (np.abs(y - obstacle_y_center) < obstacle_height / 2)
        ),
        (NX, NY),
    )
    print(
        f"Obstacle: {obstacle_length}x{obstacle_height} at ({obstacle_x_center}, {obstacle_y_center})"
    )

    # Initialize velocity field with slight perturbation to trigger instability
    vel = np.fromfunction(
        lambda d, x, y: (1 - d) * ULB * (1.0 + 1e-4 * np.sin(y / LY * 2 * np.pi)),
        (2, NX, NY),
    )

    # Initialize distribution function to equilibrium
    fin = equilibrium(1.0, vel).copy()

    # Storage for visualization snapshots
    snapshots = []
    snapshot_interval = 100  # Save every N time steps

    print("\nRunning simulation...")

    for time_step in range(MAX_ITER):
        # -----------------------------------------------------------------
        # Step 1: Apply outlet boundary condition (right wall)
        # Zero-gradient extrapolation for outgoing populations
        # -----------------------------------------------------------------
        fin[I1, -1, :] = fin[I1, -2, :]

        # -----------------------------------------------------------------
        # Step 2: Compute macroscopic density and velocity
        # -----------------------------------------------------------------
        rho, u = compute_macroscopic(fin)

        # -----------------------------------------------------------------
        # Step 3: Apply inlet boundary condition (left wall)
        # Zou/He velocity boundary condition
        # -----------------------------------------------------------------
        # Prescribe inlet velocity
        u[:, 0, :] = vel[:, 0, :]

        # Compute density from known populations (Zou/He)
        # Note: u[0,0,:] should always be < 1 for stability (typically < 0.1)
        # Add small safety margin to prevent division by zero
        velocity_factor = np.maximum(1.0 - u[0, 0, :], 1e-10)
        rho[0, :] = (
            1.0
            / velocity_factor
            * (np.sum(fin[I2, 0, :], axis=0) + 2.0 * np.sum(fin[I1, 0, :], axis=0))
        )

        # Compute equilibrium and set unknown populations
        feq = equilibrium(rho, u)
        fin[I3, 0, :] = fin[I1, 0, :] + feq[I3, 0, :] - fin[I1, 0, :]

        # -----------------------------------------------------------------
        # Step 4: BGK Collision step
        # f_out = f_in - omega * (f_in - f_eq)
        # -----------------------------------------------------------------
        fout = fin - OMEGA * (fin - feq)

        # -----------------------------------------------------------------
        # Step 5: Apply no-slip boundary at obstacle (bounce-back)
        # -----------------------------------------------------------------
        for i in range(Q):
            fout[i, obstacle] = fin[NOSLIP[i], obstacle]

        # -----------------------------------------------------------------
        # Step 6: Streaming step
        # Move distributions to neighboring nodes
        # -----------------------------------------------------------------
        for i in range(Q):
            fin[i, :, :] = np.roll(
                np.roll(fout[i, :, :], C[i, 0], axis=0), C[i, 1], axis=1
            )

        # -----------------------------------------------------------------
        # Step 7: Data collection
        # -----------------------------------------------------------------
        if time_step % snapshot_interval == 0:
            u_magnitude = np.sqrt(u[0] ** 2 + u[1] ** 2)
            # Transpose for visualization (y, x) ordering
            snapshots.append(u_magnitude.T.copy())
            print(
                f"  Time step {time_step:5d}/{MAX_ITER}: "
                f"max|u| = {np.max(u_magnitude):.4f}"
            )

    print("\nSimulation complete!")
    print(f"Collected {len(snapshots)} snapshots for visualization")
    print("=" * 70)

    # Launch visualization
    vtk_velocity_plot(NX, NY, snapshots)


if __name__ == "__main__":
    run_simulation()
