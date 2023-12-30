import numpy as np
import vtk

# Configuration and Constants
MAX_ITER: int = 3000  # Total number of time iterations
REYNOLDS_NUMBER: float = 220.0
NX: int = 520
NY: int = 180
LY: float = NY - 1.0
Q: int = 9  # Lattice dimensions and populations
ULB: float = 0.04  # Velocity in lattice units
NU_LB: float = ULB * (NY / 9) / REYNOLDS_NUMBER
OMEGA: float = 1.0 / (3.0 * NU_LB + 0.5)  # Relaxation parameter

# Lattice Constants
C: np.ndarray = np.array(
    [(x, y) for x in [0, -1, 1] for y in [0, -1, 1]]
)  # Lattice velocities
T: np.ndarray = 1.0 / 36.0 * np.ones(Q)  # Lattice weights
T[np.linalg.norm(C, axis=1) < 1.1] = 1.0 / 9.0
T[0] = 4.0 / 9.0
NOSLIP: list[int] = [np.where(np.all(C == -C[i], axis=1))[0][0] for i in range(Q)]


I1: np.ndarray = np.arange(Q)[[ci[0] < 0 for ci in C]]  # Unknown on right wall
I2: np.ndarray = np.arange(Q)[[ci[0] == 0 for ci in C]]  # Vertical middle
I3: np.ndarray = np.arange(Q)[[ci[0] > 0 for ci in C]]  # Unknown on left wall


def equilibrium(rho: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Equilibrium distribution function."""
    cu = 3.0 * np.dot(C, u.transpose(1, 0, 2))
    usqr = 3.0 / 2.0 * (u[0] ** 2 + u[1] ** 2)
    feq = np.zeros((Q, NX, NY))
    for i in range(Q):
        feq[i, :, :] = rho * T[i] * (1.0 + cu[i] + 0.5 * cu[i] ** 2 - usqr)
    return feq


# Function to create vtkPoints and velocity values
def create_points(
    nx: int, ny: int, u: np.ndarray
) -> tuple[vtk.vtkPoints, vtk.vtkFloatArray]:
    """Create vtkPoints and velocity values."""
    points = vtk.vtkPoints()
    velocity_values = vtk.vtkFloatArray()
    velocity_values.SetName("Velocity")

    # Added checks for array dimensions
    assert (
        u.shape[0] == ny and u.shape[1] == nx
    ), f"Dimension mismatch in u array {u.shape}"

    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(i, j, 0)
            velocity_values.InsertNextValue(u[j, i])  # Ensuring correct indexing

    return points, velocity_values


# Function to create vtkPolyData for visualization
def create_poly_data(
    nx: int, ny: int, points: vtk.vtkPoints, velocity_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """Create vtkPolyData for visualization."""
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    cells = vtk.vtkCellArray()
    for i in range(ny - 1):
        for j in range(nx - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, i * nx + j)
            quad.GetPointIds().SetId(1, (i + 1) * nx + j)
            quad.GetPointIds().SetId(2, (i + 1) * nx + (j + 1))
            quad.GetPointIds().SetId(3, i * nx + (j + 1))
            cells.InsertNextCell(quad)

    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(velocity_values)

    return polyData


def create_color_map(u: np.ndarray) -> vtk.vtkLookupTable:
    """Create a color map for VTK visualization."""
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(u), np.max(u)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


def vtk_velocity_plot(nx: int, ny: int, snapshots: list[np.ndarray]):
    """Visualize the velocity field using VTK."""
    # Calculate the global minimum and maximum velocity values across all snapshots
    global_min_u = min(np.min(u) for u in snapshots)
    global_max_u = max(np.max(u) for u in snapshots)

    # Setup renderer and render window
    renderer = vtk.vtkRenderer()
    renderer.SetUseFXAA(True)  # Enable or disable FXAA anti-aliasing
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Initialize variables to store actors
    surface_actor = None
    scalar_bar = None

    # Define a callback function for the timer
    def update_visualization(obj, event):
        nonlocal surface_actor, scalar_bar
        nonlocal snapshots

        # Remove the current actors
        if surface_actor:
            renderer.RemoveActor(surface_actor)
            renderer.RemoveActor(scalar_bar)

        # Stop the timer if no more snapshots are available
        if not snapshots:
            return

        # Get the next snapshot
        u_np = snapshots.pop(0)

        points, velocity_values = create_points(ny, nx, u_np)
        polyData = create_poly_data(nx, ny, points, velocity_values)
        color_map = create_color_map(u_np)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)
        mapper.SetLookupTable(color_map)
        mapper.SetScalarModeToUsePointData()
        mapper.SetScalarRange(global_min_u, global_max_u)  # Use global min/max

        # Create and add new actors
        surface_actor = vtk.vtkActor()
        surface_actor.SetMapper(mapper)
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(color_map)
        scalar_bar.SetTitle("Velocity")
        scalar_bar.SetNumberOfLabels(5)

        renderer.AddActor(surface_actor)
        renderer.AddActor(scalar_bar)
        renderer.ResetCamera()
        renderWindow.Render()

    # Add a timer to update the visualization
    renderWindowInteractor.AddObserver("TimerEvent", update_visualization)
    timerId = renderWindowInteractor.CreateRepeatingTimer(1000)

    # Start the interaction
    renderWindow.Render()
    renderWindowInteractor.Start()

    # Clean up
    renderWindow.Finalize()
    renderWindowInteractor.TerminateApp()


def run_simulation():
    """Main simulation loop with data collection."""
    rod_length, rod_width = 50, 10
    rod_x_center, rod_y_center = NX / 4, NY / 2
    obstacle = np.fromfunction(
        lambda x, y: (np.abs(x - rod_x_center) < rod_length / 2)
        & (np.abs(y - rod_y_center) < rod_width / 2),
        (NX, NY),
    )
    vel = np.fromfunction(
        lambda d, x, y: (1 - d) * ULB * (1.0 + 1e-4 * np.sin(y / LY * 2 * np.pi)),
        (2, NX, NY),
    )
    fin = equilibrium(1.0, vel).copy()

    snapshots = []
    for time in range(MAX_ITER):
        # Right wall: outflow condition
        fin[I1, -1, :] = fin[I1, -2, :]

        # Calculate macroscopic density and velocity
        rho = np.sum(fin, axis=0)
        u = np.dot(C.transpose(), fin.transpose((1, 0, 2))) / rho

        # Left wall: compute density from known populations
        u[:, 0, :] = vel[:, 0, :]
        rho[0, :] = (
            1.0
            / (1.0 - u[0, 0, :])
            * (np.sum(fin[I2, 0, :], axis=0) + 2.0 * np.sum(fin[I1, 0, :], axis=0))
        )

        # Left wall: Zou/He boundary condition
        feq = equilibrium(rho, u)
        fin[I3, 0, :] = fin[I1, 0, :] + feq[I3, 0, :] - fin[I1, 0, :]

        # Collision step
        fout = fin - OMEGA * (fin - feq)
        for i in range(Q):
            fout[i, obstacle] = fin[NOSLIP[i], obstacle]

        # Streaming step
        for i in range(Q):
            fin[i, :, :] = np.roll(
                np.roll(fout[i, :, :], C[i, 0], axis=0), C[i, 1], axis=1
            )

        # Debugging information
        if time % 300 == 0:
            print(f"Time step: {time}")
            u_magnitude = np.sqrt(u[0] ** 2 + u[1] ** 2)
            snapshots.append(u_magnitude)

    vtk_velocity_plot(NX, NY, snapshots)


if __name__ == "__main__":
    run_simulation()
