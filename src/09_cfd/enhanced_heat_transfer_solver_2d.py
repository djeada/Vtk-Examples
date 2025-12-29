"""
2D Steady-State Heat Transfer Solver

Part 4 of Heat Transfer Series: Two-Dimensional Extension

This module extends the heat transfer analysis to two dimensions, solving the 2D
steady-state heat conduction equation in a square domain using iterative methods.

Series Overview:
- Part 1: Fixed temperature BCs - simplest case with Dirichlet conditions
- Part 2: Convective BCs - more realistic with convection at boundaries
- Part 3: Temperature-dependent conductivity - non-linear problem
- Part 4 (this file): 2D extension - full two-dimensional heat transfer

Workflow:
1. Initialize a 2D grid with specified dimensions and properties.
2. Set boundary conditions (fixed temperature or symmetry).
3. Iteratively solve using Gauss-Seidel method (point-by-point).
4. Apply boundary conditions at each iteration.
5. Check residual for convergence.
6. Visualize the 2D temperature field.

Physics:
This models 2D steady-state heat conduction in a square domain with:
- Uniform thermal conductivity k (W/m·K)
- Internal volumetric heat generation S (W/m³)
- Mixed boundary conditions:
  * Fixed temperature (Dirichlet) on top and right edges
  * Zero flux/symmetry (Neumann) on bottom and left edges

Governing Equation:
    ∂/∂x(k ∂T/∂x) + ∂/∂y(k ∂T/∂y) + S = 0

For constant k:
    k(∂²T/∂x² + ∂²T/∂y²) + S = 0

Discretization (central difference):
For interior node (i,j):
    Aw*T[i,j-1] + Ae*T[i,j+1] + An*T[i+1,j] + As*T[i-1,j] - Ap*T[i,j] = -S*dx*dy

Where:
    Aw = Ae = k*dy/dx
    An = As = k*dx/dy
    Ap = Aw + Ae + An + As
"""

import matplotlib.pyplot as plt
import numpy as np
import vtk


def tdma_solve(
    a: np.ndarray, ad: np.ndarray, au: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """
    Solves a tridiagonal matrix equation using the Thomas algorithm (TDMA).

    :param a: Main diagonal elements (numpy.ndarray)
    :param ad: Sub-diagonal (lower) elements (numpy.ndarray)
    :param au: Super-diagonal (upper) elements (numpy.ndarray)
    :param b: Right-hand side vector (numpy.ndarray)
    :return: Solution vector (numpy.ndarray)
    """
    n = len(b)
    # Make copies to avoid modifying original arrays
    a = a.copy()
    b = b.copy()
    T = np.zeros(n)

    # Forward elimination
    for i in range(1, n):
        f = ad[i] / a[i - 1]
        a[i] -= f * au[i - 1]
        b[i] -= f * b[i - 1]

    # Back substitution
    T[-1] = b[-1] / a[-1]
    for i in range(n - 2, -1, -1):
        T[i] = (b[i] - au[i] * T[i + 1]) / a[i]

    return T


def two_dimensional_heat_transfer_solver(
    n: int, L: float, k: float, S: float, T_boundary: float, Rtol: float
) -> tuple:
    """
    Solves the 2D steady-state heat conduction problem using Gauss-Seidel iteration.

    Boundary Conditions:
    - Top (y=L) and Right (x=L): Fixed temperature T_boundary
    - Bottom (y=0) and Left (x=0): Zero flux (symmetry)

    This represents a quarter of a symmetric domain with heat generation.

    :param n: Number of divisions in each direction (int)
    :param L: Length of the domain (float, meters)
    :param k: Thermal conductivity (float, W/m·K)
    :param S: Volumetric heat generation rate (float, W/m³)
    :param T_boundary: Fixed temperature on Dirichlet boundaries (float, °C)
    :param Rtol: Residual tolerance for convergence (float)
    :return: Tuple of (x, y, T, iterations) - coordinate vectors, temperature field, iteration count
    """
    dx = dy = L / n
    x = np.linspace(0, L, n + 1)
    y = np.linspace(0, L, n + 1)

    # Initial guess
    T = T_boundary * np.ones((n + 1, n + 1))
    Rmax = 1.0
    iterations = 0
    max_iterations = 50000

    while Rmax > Rtol and iterations < max_iterations:
        Rmax = 0.0
        iterations += 1

        # Point-by-Point Gauss-Seidel iteration for interior nodes
        for i in range(1, n):
            for j in range(1, n):
                # Coefficients for the 5-point stencil
                Aw = Ae = k * dy / dx  # West and East
                An = As = k * dx / dy  # North and South

                # Central coefficient
                Ap = Aw + Ae + An + As

                # Source term
                b_ij = S * dx * dy

                # Calculate new temperature at (i, j)
                T_new = (
                    Aw * T[i, j - 1]
                    + Ae * T[i, j + 1]
                    + An * T[i + 1, j]
                    + As * T[i - 1, j]
                    + b_ij
                ) / Ap

                # Update residual and temperature
                R = abs(T_new - T[i, j])
                Rmax = max(R, Rmax)
                T[i, j] = T_new

        # Apply Boundary Conditions
        # Dirichlet BCs: Top and Right edges have fixed temperature
        T[:, n] = T_boundary  # Right edge (x = L)
        T[n, :] = T_boundary  # Top edge (y = L)

        # Neumann BCs (zero flux): Left and Bottom edges (symmetry)
        T[:, 0] = T[:, 1]  # Left edge (x = 0): dT/dx = 0
        T[0, :] = T[1, :]  # Bottom edge (y = 0): dT/dy = 0

        if iterations % 1000 == 0:
            print(f"Iteration: {iterations}, Residual: {Rmax:.2e}")

    return x, y, T, iterations


def plot_temperature_profile(X: np.ndarray, Y: np.ndarray, T: np.ndarray, S: float):
    """
    Plots the 2D temperature distribution with contours.

    :param X: X-coordinate grid (numpy.ndarray)
    :param Y: Y-coordinate grid (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param S: Heat generation rate (float)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Filled contour plot
    levels = 50
    cf = ax1.contourf(X, Y, T, levels, cmap="coolwarm")
    plt.colorbar(cf, ax=ax1, label="Temperature (°C)")
    ax1.set_xlabel("x (m)", fontsize=12)
    ax1.set_ylabel("y (m)", fontsize=12)
    ax1.set_title("2D Temperature Distribution", fontsize=14)
    ax1.set_aspect("equal")

    # Contour lines
    cs = ax2.contour(X, Y, T, 20, colors="k", linewidths=0.5)
    ax2.clabel(cs, inline=True, fontsize=8)
    cf2 = ax2.contourf(X, Y, T, levels, cmap="coolwarm", alpha=0.7)
    plt.colorbar(cf2, ax=ax2, label="Temperature (°C)")
    ax2.set_xlabel("x (m)", fontsize=12)
    ax2.set_ylabel("y (m)", fontsize=12)
    ax2.set_title(f"Temperature Contours (S = {S:.0e} W/m³)", fontsize=14)
    ax2.set_aspect("equal")

    plt.tight_layout()
    plt.show()


def create_points(x: np.ndarray, y: np.ndarray, T: np.ndarray) -> tuple:
    """
    Creates points and temperature values for VTK visualization.

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :return: Tuple of (points, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    nx = len(x)
    ny = len(y)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(x[i], y[j], 0)
            temperature_values.InsertNextValue(T[j, i])

    return points, temperature_values


def create_poly_data(
    nx: int, ny: int, points: vtk.vtkPoints, temperature_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """
    Creates vtkPolyData for 2D visualization.

    :param nx: Number of points in x direction (int)
    :param ny: Number of points in y direction (int)
    :param points: vtkPoints
    :param temperature_values: vtkFloatArray
    :return: vtkPolyData
    """
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    cells = vtk.vtkCellArray()
    for j in range(ny - 1):
        for i in range(nx - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, j * nx + i)
            quad.GetPointIds().SetId(1, j * nx + (i + 1))
            quad.GetPointIds().SetId(2, (j + 1) * nx + (i + 1))
            quad.GetPointIds().SetId(3, (j + 1) * nx + i)
            cells.InsertNextCell(quad)

    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(temperature_values)

    return polyData


def create_color_map(T: np.ndarray) -> vtk.vtkLookupTable:
    """
    Creates a color map for vtk visualization.

    :param T: Temperature field (numpy.ndarray)
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(T), np.max(T)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


def create_boundary_labels(L: float) -> list:
    """
    Create text labels for boundary conditions.

    :param L: Domain size (float)
    :return: List of vtkBillboardTextActor3D actors
    """
    labels = []

    # Create boundary condition labels
    bc_info = [
        ("T = T_boundary\n(Dirichlet)", (L + 0.15, L / 2, 0), (1.0, 0.4, 0.4)),  # Right
        ("T = T_boundary\n(Dirichlet)", (L / 2, L + 0.15, 0), (1.0, 0.4, 0.4)),  # Top
        ("dT/dx = 0\n(Symmetry)", (-0.15, L / 2, 0), (0.4, 0.8, 1.0)),  # Left
        ("dT/dy = 0\n(Symmetry)", (L / 2, -0.15, 0), (0.4, 0.8, 1.0)),  # Bottom
    ]

    for text, pos, color in bc_info:
        text_actor = vtk.vtkBillboardTextActor3D()
        text_actor.SetInput(text)
        text_actor.SetPosition(pos)
        text_actor.GetTextProperty().SetFontSize(12)
        text_actor.GetTextProperty().SetColor(color)
        text_actor.GetTextProperty().SetJustificationToCentered()
        text_actor.GetTextProperty().SetBold(True)
        labels.append(text_actor)

    return labels


def vtk_3d_surface_plot(
    x: np.ndarray,
    y: np.ndarray,
    T: np.ndarray,
    L: float,
    S: float,
    T_boundary: float,
    title: str = "2D Heat Transfer with Internal Generation",
):
    """
    Creates an advanced 3D visualization of the 2D temperature field.

    Features:
    - 3D surface with height representing temperature
    - Color mapping for temperature
    - Boundary condition annotations
    - Improved lighting and camera

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param L: Domain size (float)
    :param S: Heat generation rate (float)
    :param T_boundary: Boundary temperature (float)
    :param title: Window title (str)
    """
    # Create points with Z = (T - T_min) / scale for 3D effect
    T_min, T_max = np.min(T), np.max(T)
    height_scale = L / (T_max - T_min) if T_max > T_min else 1.0

    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    nx, ny = len(x), len(y)
    for j in range(ny):
        for i in range(nx):
            z = (T[j, i] - T_min) * height_scale * 0.5
            points.InsertNextPoint(x[i], y[j], z)
            temperature_values.InsertNextValue(T[j, i])

    polyData = create_poly_data(nx, ny, points, temperature_values)
    color_map = create_color_map(T)

    # Mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(T_min, T_max)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)
    surface_actor.GetProperty().SetInterpolationToPhong()
    surface_actor.GetProperty().SetSpecular(0.3)
    surface_actor.GetProperty().SetSpecularPower(20)

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.08)
    scalar_bar.SetHeight(0.4)
    scalar_bar.SetPosition(0.9, 0.3)

    # Title
    title_actor = vtk.vtkTextActor()
    title_actor.SetInput(title)
    title_actor.GetTextProperty().SetFontSize(18)
    title_actor.GetTextProperty().SetColor(1, 1, 1)
    title_actor.GetTextProperty().SetBold(True)
    title_actor.SetPosition(10, 10)

    # Info text
    info_text = f"S = {S:.0e} W/m³\nT_boundary = {T_boundary}°C\nT_max = {T_max:.1f}°C"
    info_actor = vtk.vtkTextActor()
    info_actor.SetInput(info_text)
    info_actor.GetTextProperty().SetFontSize(14)
    info_actor.GetTextProperty().SetColor(0.9, 0.9, 0.9)
    info_actor.SetPosition(10, 60)

    # Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)
    renderer.AddActor2D(info_actor)

    # Add boundary labels
    for label in create_boundary_labels(L):
        renderer.AddActor(label)

    # Background gradient
    renderer.SetBackground(0.1, 0.1, 0.15)
    renderer.SetBackground2(0.02, 0.02, 0.05)
    renderer.GradientBackgroundOn()

    # Lighting
    light = vtk.vtkLight()
    light.SetPosition(L / 2, L / 2, L * 2)
    light.SetFocalPoint(L / 2, L / 2, 0)
    light.SetIntensity(1.0)
    renderer.AddLight(light)

    # Render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1200, 800)
    renderWindow.SetWindowName(title)

    # Interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interact_style)

    # Axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0, 0, 0.2, 0.2)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Camera position for good 3D view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(L * 2.5, -L * 1.5, L * 2)
    camera.SetFocalPoint(L / 2, L / 2, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    """
    Main function demonstrating 2D steady-state heat transfer.

    Part 4 of the Heat Transfer Series: Two-dimensional extension
    """
    # Physical parameters
    L = 1.0  # Domain size (m)
    n = 40  # Number of divisions
    k = 50.0  # Thermal conductivity (W/m·K)
    S = 5e4  # Heat generation rate (W/m³)
    T_boundary = 100.0  # Fixed boundary temperature (°C)
    Rtol = 1e-5  # Convergence tolerance

    print("=" * 70)
    print("Part 4: 2D Steady-State Heat Transfer Solver")
    print("=" * 70)
    print(f"Domain: {L}m × {L}m, Grid: {n+1} × {n+1}")
    print(f"Thermal conductivity: k = {k} W/m·K")
    print(f"Heat generation: S = {S:.2e} W/m³")
    print(f"Boundary temperature: T = {T_boundary}°C (Dirichlet on top/right)")
    print("Symmetry boundaries on left/bottom (zero flux)")
    print("-" * 70)
    print("Solving...")

    # Solve the 2D heat transfer problem
    x, y, T, iterations = two_dimensional_heat_transfer_solver(
        n, L, k, S, T_boundary, Rtol
    )

    # Create meshgrid for plotting
    X, Y = np.meshgrid(x, y)

    # Print results
    print("-" * 70)
    print(f"Converged in {iterations} iterations")
    print(f"Maximum temperature: {np.max(T):.2f}°C at corner (0,0)")
    print(f"Minimum temperature: {np.min(T):.2f}°C (at boundaries)")
    print(f"Temperature at center: {T[n//2, n//2]:.2f}°C")
    print("=" * 70)

    # Plot temperature profile
    plot_temperature_profile(X, Y, T, S)

    # 3D VTK visualization
    vtk_3d_surface_plot(x, y, T, L, S, T_boundary)


if __name__ == "__main__":
    main()
