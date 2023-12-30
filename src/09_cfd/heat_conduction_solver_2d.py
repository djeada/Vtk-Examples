"""
This Python module provides an advanced numerical solver for two-dimensional heat conduction problems. It utilizes a finite difference method to solve the heat equation in a 2D domain, taking into account both conductive and convective heat transfer mechanisms.

Workflow:
- The domain is first initialized with a grid of specified dimensions.
- Boundary conditions and initial conditions are set up based on the physical problem.
- The heat equation is then solved iteratively using a Gauss-Seidel method with over-relaxation for better convergence.
- During each iteration, the temperature distribution is updated based on the neighboring node values, considering both conduction and convection effects.
- The iterative process continues until the solution converges to a specified residual tolerance.

Physics and Mathematics:
- The module models heat transfer in a 2D domain, considering both conduction and convection.
- Conduction is modeled using Fourier's law, and convective effects are incorporated through a velocity field.
- The discretization of the heat equation is done using a finite difference method, specifically focusing on a staggered grid approach.
- The solver handles various boundary conditions, including fixed temperature and heat flux boundaries.
- It also includes the calculation of heat transfer coefficients and Nusselt numbers for further analysis of the heat transfer process.
"""
import matplotlib.pyplot as plt
import numpy as np
import vtk


def initialize_grid(l: float, h: float, dx: float, dy: float) -> tuple:
    """
    Initialize the grid for the simulation.

    Args:
        l (float): Length of the domain (x-direction).
        h (float): Height of the domain (y-direction).
        dx (float): Grid spacing in the x-direction.
        dy (float): Grid spacing in the y-direction.

    Returns:
        tuple: A tuple containing (nx, ny, x, y).
    """
    nx, ny = int(l / dx), int(h / dy)
    x, y = np.arange(0, l + dx, dx), np.arange(0, h + dy, dy)
    return nx, ny, x, y


def initialize_temperature_array(ny: int, nx: int, Ti: float) -> np.ndarray:
    """
    Initialize the temperature array.

    Args:
        ny (int): Number of grid points in the y-direction.
        nx (int): Number of grid points in the x-direction.
        Ti (float): Initial temperature.

    Returns:
        np.ndarray: Initialized temperature array.
    """
    return Ti * np.ones((ny + 1, nx + 1))


def velocity(yp: float, h: float) -> float:
    """
    Calculate velocity.

    Args:
        yp (float): Y-coordinate for velocity calculation.
        h (float): Height of the domain (y-direction).

    Returns:
        float: Velocity value.
    """
    return 1.5 * (1 - (4 * yp / h) ** 2)


def update_temperature(
    Aw: float,
    Ae: float,
    An: float,
    As: float,
    b: float,
    T: np.ndarray,
    i: int,
    j: int,
    nx: int,
    ny: int,
) -> float:
    """
    Update temperature at a grid point.

    Args:
        Aw (float): West diffusion coefficient.
        Ae (float): East diffusion coefficient.
        An (float): North diffusion coefficient.
        As (float): South diffusion coefficient.
        b (float): Source term.
        T (np.ndarray): Temperature array.
        i (int): Row index of the grid point.
        j (int): Column index of the grid point.
        nx (int): Number of grid points in the x-direction.
        ny (int): Number of grid points in the y-direction.

    Returns:
        float: Absolute value of the temperature update.
    """
    Ap = Aw + Ae + An + As
    tmp = b - Ap * T[i, j]
    if j > 0:
        tmp += Aw * T[i, j - 1]
    if j < nx:
        tmp += Ae * T[i, j + 1]
    if i < ny:
        tmp += An * T[i + 1, j]
    if i > 0:
        tmp += As * T[i - 1, j]
    T[i, j] += tmp / Ap
    return abs(tmp)


def apply_boundary_conditions(
    T: np.ndarray,
    Ti: float,
    nx: int,
    ny: int,
    k: float,
    dx: float,
    dy: float,
    hf: float,
):
    """
    Apply boundary conditions to the temperature array.

    Args:
        T (np.ndarray): Temperature array.
        Ti (float): Inlet temperature.
        nx (int): Number of grid points in the x-direction.
        ny (int): Number of grid points in the y-direction.
        k (float): Thermal conductivity.
        dx (float): Grid spacing in the x-direction.
        dy (float): Grid spacing in the y-direction.
        hf (float): Heat flux.
    """
    T[:, 0], T[:, -1] = Ti, T[:, -2]  # Left (Inlet) and Right (Developed Flow)
    for i in range(1, nx):
        Aw, Ae, An, As = k * dy / (2 * dx), k * dy / (2 * dx), k * dx / dy, 0
        b = hf * dx
        update_temperature(Aw, Ae, An, As, b, T, 0, i, nx, ny)
        update_temperature(Aw, Ae, An, As, b, T, ny, i, nx, ny)


def run_simulation(
    T: np.ndarray,
    nx: int,
    ny: int,
    Rtol: float,
    dx: float,
    dy: float,
    y: np.ndarray,
    Ti: float,
    hf: float,
    rho: float,
    c: float,
    k: float,
    h: float,
) -> np.ndarray:
    """
    Run the heat transfer simulation.

    Args:
        T (np.ndarray): Temperature array.
        nx (int): Number of grid points in the x-direction.
        ny (int): Number of grid points in the y-direction.
        Rtol (float): Residual tolerance.
        dx (float): Grid spacing in the x-direction (m).
        dy (float): Grid spacing in the y-direction (m).
        y (np.ndarray): Y-coordinates (m).
        Ti (float): Inlet temperature (°C).
        hf (float): Heat transfer coefficient (W/(m^2*K)).
        rho (float): Density (kg/m^3).
        c (float): Specific heat capacity (J/(kg*K)).
        k (float): Thermal conductivity (W/(m*K)).
        h (float): Height of the domain in the y-direction (m).

    Returns:
        np.ndarray: Updated temperature array.
    """
    iteration, Rmax = 0, 1

    while Rmax > Rtol:
        Rmax = 0

        for i in range(1, ny):
            for j in range(1, nx):
                F = rho * c * velocity(y[i] - dy / 2, h)
                Aw = max(F, k / dx + F / 2, 0) * dy
                Ae = max(-F, k / dx - F / 2, 0) * dy
                An, As = k * dx / dy, k * dx / dy
                Rmax = max(Rmax, update_temperature(Aw, Ae, An, As, 0, T, i, j, nx, ny))

        apply_boundary_conditions(T, Ti, nx, ny, k, dx, dy, hf)
        iteration += 1

    return T


def calculate_bulk_temperature(
    T: np.ndarray, y: np.ndarray, dy: float, h: float, nx: int
) -> np.ndarray:
    """
    Calculate bulk temperature.

    Args:
        T (np.ndarray): Temperature array.
        y (np.ndarray): Y-coordinates.
        dy (float): Grid spacing in the y-direction.
        h (float): Height of the domain (y-direction).
        nx (int): Number of grid points in the x-direction.

    Returns:
        np.ndarray: Bulk temperature values.
    """
    Tb = np.zeros(nx + 1)
    for i in range(nx + 1):
        velocity_profile = np.array([velocity(yp, h) for yp in y[:-1]])
        Tb[i] = np.sum(T[:-1, i] * velocity_profile * dy) / np.sum(
            velocity_profile * dy
        )
    return Tb


def calculate_heat_flux(T: np.ndarray, k: float, dy: float) -> np.ndarray:
    """
    Calculate heat flux at the wall.

    Args:
        T (np.ndarray): Temperature array.
        k (float): Thermal conductivity.
        dy (float): Grid spacing in the y-direction.

    Returns:
        np.ndarray: Heat flux values.
    """
    q_w = -k * (T[1, :] - T[0, :]) / dy
    return q_w


def calculate_ht_and_nusselt(
    Tb: np.ndarray, q_w: np.ndarray, k: float, h: float, T: np.ndarray
) -> tuple:
    """
    Calculate heat transfer coefficient and Nusselt number.

    Args:
        Tb (np.ndarray): Bulk temperature values.
        q_w (np.ndarray): Heat flux at the wall.
        k (float): Thermal conductivity.
        h (float): Height of the domain (y-direction).
        T (np.ndarray): Temperature array.

    Returns:
        tuple: Tuple containing heat transfer coefficient and Nusselt number.
    """
    ht = np.zeros_like(Tb)
    Nu = np.zeros_like(Tb)
    for i in range(len(Tb)):
        ht[i] = q_w[i] / (T[0, i] - Tb[i])
        Nu[i] = ht[i] * (2 * h) / k
    return ht, Nu


def plot_results(
    x: np.ndarray,
    y: np.ndarray,
    T: np.ndarray,
    Ti: float,
    Tb: np.ndarray,
    ht: np.ndarray,
    Nu: np.ndarray,
):
    """
    Plot simulation results.

    Args:
        x (np.ndarray): X-coordinates.
        y (np.ndarray): Y-coordinates.
        T (np.ndarray): Temperature array.
        Ti (float): Inlet temperature.
        Tb (np.ndarray): Bulk temperature values.
        ht (np.ndarray): Heat transfer coefficients.
        Nu (np.ndarray): Nusselt numbers.
    """
    plt.figure(figsize=(12, 8))
    plot_titles = [
        "Temperature Distribution",
        "Bulk Temperature",
        "Heat Transfer Coefficient",
        "Nusselt Number",
    ]

    for i in range(4):
        plt.subplot(2, 2, i + 1)
        if i == 0:
            plt.contourf(x, y, T, levels=np.linspace(Ti, np.max(T), 50))
            plt.colorbar()
        elif i == 1:
            plt.plot(x, Tb)
            plt.ylabel("Bulk Temperature (°C)")
        elif i == 2:
            plt.plot(x, ht)
            plt.ylabel("Heat Transfer Coefficient")
        elif i == 3:
            plt.plot(x, Nu)
            plt.ylabel("Nusselt Number")

        plt.xlabel("x (m)")
        plt.title(plot_titles[i])

    plt.tight_layout()
    plt.show()


def create_points(x: np.ndarray, y: np.ndarray, T: np.ndarray) -> tuple:
    """
    Create points and temperature values for VTK visualization.

    Args:
        x (np.ndarray): X-coordinates.
        y (np.ndarray): Y-coordinates.
        T (np.ndarray): Temperature array.

    Returns:
        tuple: Tuple containing vtkPoints and vtkFloatArray for temperature values.
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    # Adjusted loop to iterate over both x and y
    for i in range(len(y)):  # Iterate over y
        for j in range(len(x)):  # Iterate over x
            points.InsertNextPoint(x[j], y[i], 0)
            temperature_values.InsertNextValue(T[i, j])

    return points, temperature_values


def create_poly_data(
    n_x: int, n_y: int, points: vtk.vtkPoints, temperature_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """
    Create vtkPolyData for VTK visualization.

    Args:
        n_x (int): Number of grid points in the x-direction.
        n_y (int): Number of grid points in the y-direction.
        points (vtk.vtkPoints): vtkPoints containing grid points.
        temperature_values (vtk.vtkFloatArray): Temperature values.

    Returns:
        vtk.vtkPolyData: vtkPolyData for visualization.
    """
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    cells = vtk.vtkCellArray()
    for i in range(n_y - 1):
        for j in range(n_x - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, i * n_x + j)
            quad.GetPointIds().SetId(1, (i + 1) * n_x + j)
            quad.GetPointIds().SetId(2, (i + 1) * n_x + (j + 1))
            quad.GetPointIds().SetId(3, i * n_x + (j + 1))
            cells.InsertNextCell(quad)

    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(temperature_values)

    return polyData


def create_color_map(
    T: np.ndarray, color_range=None, hue_range=(0.667, 0)
) -> vtk.vtkLookupTable:
    """
    Create a color map for VTK visualization.

    Args:
        T (np.ndarray): Temperature array.
        color_range (tuple): Range for color mapping (min, max).
        hue_range (tuple): Range for hue variation.

    Returns:
        vtk.vtkLookupTable: VTK lookup table for color mapping.
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(T), np.max(T)
    color_map.SetRange(color_range if color_range else (min_val, max_val))
    color_map.SetHueRange(*hue_range)
    color_map.Build()
    return color_map


def vtk_surface_plot(x: np.ndarray, y: np.ndarray, T: np.ndarray):
    """
    Visualize the temperature distribution using VTK.

    Args:
        x (np.ndarray): X-coordinates.
        y (np.ndarray): Y-coordinates.
        T (np.ndarray): Temperature array.
        color_range (tuple): Range for color mapping (min, max).
        hue_range (tuple): Range for hue variation.
    """
    points, temperature_values = create_points(x, y, T)
    polyData = create_poly_data(len(x), len(y), points, temperature_values)
    color_map = create_color_map(T)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(T), np.max(T))

    # Surface actor setup
    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)
    surface_actor.GetProperty().SetLineWidth(5)

    # Scalar bar (legend)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature")
    scalar_bar.SetNumberOfLabels(5)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(scalar_bar)

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Enable user interaction with the scene
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    # Set up interaction style for camera controls
    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interact_style)

    # Create an axes actor
    axes = vtk.vtkAxesActor()

    # Create the widget for the axes
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0, 0, 0.3, 0.3)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Reset the camera to include all visible actors
    renderer.ResetCamera()

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    # Simulation parameters
    l = 20.0
    h = 1.0
    dx = 0.1
    dy = 0.1
    rho = 1.0
    c = 100.0
    k = 1.0
    hf = 1000.0
    Ti = 50.0
    Rtol = 1e-3

    # Initialize grid and temperature arrays
    nx, ny, x, y = initialize_grid(l, h, dx, dy)
    T = initialize_temperature_array(ny, nx, Ti)

    # Run the simulation
    T = run_simulation(
        T=T,
        nx=nx,
        ny=ny,
        Rtol=Rtol,
        dx=dx,
        dy=dy,
        y=y,
        Ti=Ti,
        hf=hf,
        rho=rho,
        c=c,
        k=k,
        h=h,
    )

    # Calculate bulk temperature and heat flux
    Tb = calculate_bulk_temperature(T, y, dy, h, nx)
    q_w = calculate_heat_flux(T, k, dy)

    # Calculate heat transfer and Nusselt number
    ht, Nu = calculate_ht_and_nusselt(Tb, q_w, k, h, T)

    # Plot the results
    plot_results(x, y, T, Ti, Tb, ht, Nu)
    vtk_surface_plot(x, y, T)


if __name__ == "__main__":
    main()
