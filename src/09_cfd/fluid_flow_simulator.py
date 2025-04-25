import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtkmodules.util import numpy_support


def initialize_grid(l, h, dx, dy, u0):
    nx, ny = int(l / dx), int(h / dy)
    x, y = np.arange(0, l + dx, dx), np.arange(0, h + dy, dy)
    u = u0 * np.ones((ny + 1, nx + 1))
    u[0, :] = u[-1, :] = 0  # No-slip boundary condition
    p = np.zeros((ny + 1, nx + 1))
    return x, y, u, p, nx, ny


def update_velocity(u, p, rho, mu, dx, dy, dt, ua, nx, ny):
    bmax = 0
    u_new = np.copy(u)

    # Compute velocity field
    for i in range(1, ny):
        for j in range(1, nx):
            # Compute pressure gradient
            dpdx = (p[i, j - 1] - p[i, j + 1]) / (2 * dx)
            dpdy = (p[i - 1, j] - p[i + 1, j]) / (2 * dy)

            # Compute viscous term
            viscous_term_x = mu * (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / dx**2
            viscous_term_y = mu * (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / dy**2

            # Update velocity
            u_new[i, j] = u[i, j] + dt * (
                -u[i, j] * dpdx + viscous_term_x - dpdy + viscous_term_y
            )

            # Update max change for convergence check
            bmax = max(bmax, abs(u_new[i, j] - u[i, j]))

    return u_new, bmax


def create_points(x, y, u):
    points = vtk.vtkPoints()
    velocity_values = vtk.vtkFloatArray()
    velocity_values.SetName("Velocity")

    nx, ny = len(x), len(y)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(x[i], y[j], 0)
            velocity_values.InsertNextValue(u[j, i])  # Ensure correct indexing

    return points, velocity_values


# Create vtkPolyData for visualization
def create_poly_data(nx, ny, points, velocity_values):
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


# Create a color map for VTK visualization
def create_color_map(u):
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(u), np.max(u)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


# Visualizes the velocity field using VTK
def vtk_velocity_plot(x, y, snapshots):
    # Setup renderer and render window
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Enable user interaction with the scene
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Set up interaction style for camera controls
    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interact_style)

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
        u_vtk = snapshots.pop(0)

        # Process the snapshot
        if isinstance(u_vtk, vtk.vtkTypeFloat32Array):
            u_np = numpy_support.vtk_to_numpy(u_vtk)
            u_np = u_np.reshape(len(y), len(x))
        else:
            u_np = u_vtk

        points, velocity_values = create_points(x, y, u_np)
        polyData = create_poly_data(len(x), len(y), points, velocity_values)
        color_map = create_color_map(u_np)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)
        mapper.SetLookupTable(color_map)
        mapper.SetScalarModeToUsePointData()
        mapper.SetScalarRange(np.min(u_np), np.max(u_np))

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
    timerId = renderWindowInteractor.CreateRepeatingTimer(1000)  # Time in milliseconds

    # Start the interaction
    renderWindow.Render()
    renderWindowInteractor.Start()

    # Clean up
    renderWindow.Finalize()
    renderWindowInteractor.TerminateApp()


def main():
    # Simulation parameters
    l, h = 20.0, 2.0
    dx, dy = 0.1, 0.1
    dt = 0.01  # Time step
    rho, mu, ua = 1.0, 0.01, 0.5
    u0 = 1.0  # Inlet velocity
    btol = 1e-4  # Convergence tolerance

    # Initialize the grid and variables
    x, y, u, p, nx, ny = initialize_grid(l, h, dx, dy, u0)

    snapshots = []
    snapshot_interval = 100  # Define how often to take snapshots

    # Iterative solver
    iter_count = 0
    time_elapsed = 0

    while True:
        u, bmax = update_velocity(u, p, rho, mu, dx, dy, dt, ua, nx, ny)
        time_elapsed += dt
        print(
            f"Iteration {iter_count}, Time: {time_elapsed:.2f} s, Residual: {bmax:.6f}"
        )

        if iter_count % snapshot_interval == 0:
            # Convert u to a VTK array and store it
            vtk_u = numpy_support.numpy_to_vtk(
                num_array=u.ravel(), deep=True, array_type=vtk.VTK_FLOAT
            )
            vtk_u.SetName("Velocity")
            snapshots.append(vtk_u)

        if bmax < btol:
            print(f"Convergence achieved in {time_elapsed:.2f} seconds.")
            break

        iter_count += 1
        if iter_count > 10000:  # Prevent infinite loop
            print("Convergence not achieved")
            break

    # Plotting the velocity field
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, u, levels=np.linspace(u.min(), u.max(), 20))
    plt.colorbar()
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Velocity Field")
    plt.show()

    # Post-simulation visualization
    vtk_velocity_plot(x, y, snapshots)


if __name__ == "__main__":
    main()
