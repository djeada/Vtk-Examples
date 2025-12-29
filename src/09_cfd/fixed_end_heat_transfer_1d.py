"""
1D Temperature Distribution in a Rod with Fixed Temperature Boundary Conditions

Part 1 of Heat Transfer Series: Fixed-End Boundary Conditions

This module simulates the steady-state temperature distribution along a one-dimensional
rod with fixed temperatures at both ends and internal heat generation. This is the
simplest boundary condition case in the series.

Series Overview:
- Part 1 (this file): Fixed temperature BCs - simplest case with Dirichlet conditions
- Part 2: Convective BCs - more realistic with convection at boundaries
- Part 3: Temperature-dependent conductivity - non-linear material properties
- Part 4: 2D extension - full two-dimensional heat transfer

Workflow:
1. The rod is discretized into 'n' equal segments, defining a grid along its length.
2. For each grid point, the discretized heat equation is written, forming a tri-diagonal
   matrix system (A, Ad, Au) and a right-hand side vector 'b'.
3. Fixed temperature boundary conditions are applied at both ends.
4. The TDMA solver efficiently solves the tri-diagonal system.
5. The temperature distribution is visualized in 3D with VTK.
6. Energy balance is verified to validate the solution.

Physics:
The problem models steady-state heat conduction in a rod of length L with:
- Constant thermal conductivity k (W/m·K)
- Uniform volumetric heat generation S (W/m³)
- Fixed temperatures at both ends: T(0) = T_left, T(L) = T_right

Governing Equation:
    d/dx(k * dT/dx) + S = 0
    
For constant k:
    k * d²T/dx² + S = 0

Discretization (central difference):
For interior node i:
    -k/dx * T[i-1] + 2k/dx * T[i] - k/dx * T[i+1] = S * dx

Analytical Solution:
    T(x) = T_left + (T_right - T_left) * x/L + S/(2k) * x * (L - x)

The parabolic profile results from heat generation, with maximum temperature
occurring in the interior of the rod.
"""

import numpy as np
import vtk
from matplotlib import pyplot as plt


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


def setup_problem(
    n: int, L: float, k: float, S: float, T_left: float, T_right: float
) -> tuple:
    """
    Sets up the problem matrix and vectors for fixed-end heat transfer.

    The discretized equation for interior nodes is:
        -(k/dx)*T[i-1] + (2k/dx)*T[i] - (k/dx)*T[i+1] = S*dx

    Boundary conditions (Dirichlet):
        T[0] = T_left, T[n] = T_right

    :param n: Number of divisions in the domain (int)
    :param L: Length of the domain (float, meters)
    :param k: Thermal conductivity (float, W/m·K)
    :param S: Source term - volumetric heat generation (float, W/m³)
    :param T_left: Fixed temperature at left boundary (float, °C)
    :param T_right: Fixed temperature at right boundary (float, °C)
    :return: Tuple of ((A, Ad, Au), b) - matrix diagonals and RHS vector
    """
    dx = L / n
    cond = k / dx

    # Interior nodes: standard discretization
    A = 2 * cond * np.ones(n + 1)
    Au = -cond * np.ones(n + 1)
    Ad = -cond * np.ones(n + 1)

    # Apply boundary conditions (Dirichlet type)
    A[0] = 1.0
    Au[0] = 0.0
    Ad[0] = 0.0

    A[-1] = 1.0
    Au[-1] = 0.0
    Ad[-1] = 0.0

    # Right-hand side: source term for interior, boundary values at ends
    b = (S * dx) * np.ones(n + 1)
    b[0] = T_left
    b[-1] = T_right

    return (A, Ad, Au), b


def analytical_solution(
    x: np.ndarray, L: float, k: float, S: float, T_left: float, T_right: float
) -> np.ndarray:
    """
    Computes the analytical solution for 1D heat conduction with fixed ends
    and uniform heat generation.

    The governing equation is: k * d²T/dx² + S = 0
    With boundary conditions: T(0) = T_left, T(L) = T_right

    Solution: T(x) = T_left + (T_right - T_left) * x/L + S/(2k) * x * (L - x)

    The solution is a parabola superimposed on a linear profile.
    The linear part accounts for the boundary condition difference.
    The parabolic part accounts for heat generation.

    :param x: Position array (numpy.ndarray)
    :param L: Length of the domain (float)
    :param k: Thermal conductivity (float)
    :param S: Heat generation rate (float)
    :param T_left: Temperature at left boundary (float)
    :param T_right: Temperature at right boundary (float)
    :return: Temperature array (numpy.ndarray)
    """
    # Linear part from boundary conditions
    T_linear = T_left + (T_right - T_left) * x / L
    # Parabolic part from heat generation
    T_generation = (S / (2 * k)) * x * (L - x)
    return T_linear + T_generation


def plot_temperature_distribution(
    x: np.ndarray, T_numerical: np.ndarray, T_analytical: np.ndarray
):
    """
    Plots the numerical and analytical temperature distributions.

    :param x: Array of positions (numpy.ndarray)
    :param T_numerical: Numerical solution (numpy.ndarray)
    :param T_analytical: Analytical solution (numpy.ndarray)
    """
    plt.figure(figsize=(10, 6))
    plt.plot(x, T_numerical, "b-", linewidth=2, label="Numerical (TDMA)")
    plt.plot(x, T_analytical, "r--", linewidth=2, label="Analytical")
    plt.xlabel("Position x (m)", fontsize=12)
    plt.ylabel("Temperature (°C)", fontsize=12)
    plt.title("1D Heat Conduction with Fixed Ends and Heat Generation", fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


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


def create_cylinder_rod(
    x: np.ndarray, T: np.ndarray, radius: float = 0.05, resolution: int = 20
) -> tuple:
    """
    Create a 3D cylindrical rod representation with temperature coloring.

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param radius: Radius of the cylinder (float)
    :param resolution: Number of segments around the cylinder circumference (int)
    :return: Tuple of (polydata, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")
    cells = vtk.vtkCellArray()

    # Create cylinder surface points for each x position
    for i, (xi, Ti) in enumerate(zip(x, T)):
        for j in range(resolution):
            theta = 2.0 * np.pi * j / resolution
            y = radius * np.cos(theta)
            z = radius * np.sin(theta)
            points.InsertNextPoint(xi, y, z)
            temperature_values.InsertNextValue(Ti)

    # Create quad cells connecting adjacent rings
    for i in range(len(x) - 1):
        for j in range(resolution):
            quad = vtk.vtkQuad()
            p0 = i * resolution + j
            p1 = i * resolution + (j + 1) % resolution
            p2 = (i + 1) * resolution + (j + 1) % resolution
            p3 = (i + 1) * resolution + j

            quad.GetPointIds().SetId(0, p0)
            quad.GetPointIds().SetId(1, p1)
            quad.GetPointIds().SetId(2, p2)
            quad.GetPointIds().SetId(3, p3)
            cells.InsertNextCell(quad)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(temperature_values)

    return polyData, temperature_values


def create_end_cap(
    x_pos: float, T_val: float, radius: float = 0.05, resolution: int = 20
) -> vtk.vtkPolyData:
    """
    Create an end cap (disk) for the cylinder.

    :param x_pos: X position of the cap (float)
    :param T_val: Temperature value for the cap (float)
    :param radius: Radius of the cap (float)
    :param resolution: Number of segments (int)
    :return: vtkPolyData for the end cap
    """
    disk = vtk.vtkDiskSource()
    disk.SetInnerRadius(0)
    disk.SetOuterRadius(radius)
    disk.SetRadialResolution(1)
    disk.SetCircumferentialResolution(resolution)
    disk.Update()

    transform = vtk.vtkTransform()
    transform.Translate(x_pos, 0, 0)
    transform.RotateY(90)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(disk.GetOutput())
    transformFilter.Update()

    output = transformFilter.GetOutput()
    temp_array = vtk.vtkFloatArray()
    temp_array.SetName("Temperature")
    for _ in range(output.GetNumberOfPoints()):
        temp_array.InsertNextValue(T_val)
    output.GetPointData().SetScalars(temp_array)

    return output


def create_boundary_annotation(
    text: str, position: tuple, color: tuple = (1, 1, 1)
) -> vtk.vtkBillboardTextActor3D:
    """
    Create a 3D text annotation for boundary conditions.

    :param text: Text to display (str)
    :param position: Position (x, y, z)
    :param color: Text color (r, g, b)
    :return: vtkBillboardTextActor3D
    """
    text_actor = vtk.vtkBillboardTextActor3D()
    text_actor.SetInput(text)
    text_actor.SetPosition(position)
    text_actor.GetTextProperty().SetFontSize(14)
    text_actor.GetTextProperty().SetColor(color)
    text_actor.GetTextProperty().SetJustificationToCentered()
    text_actor.GetTextProperty().SetBold(True)

    return text_actor


def create_heat_source_indicator(L: float, radius: float) -> vtk.vtkActor:
    """
    Create visual indicators for internal heat generation.

    :param L: Length of the rod (float)
    :param radius: Radius of the rod (float)
    :return: vtkActor showing heat generation symbols
    """
    # Create small spheres along the rod to indicate heat generation
    append = vtk.vtkAppendPolyData()

    n_indicators = 5
    for i in range(n_indicators):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(radius * 0.15)
        sphere.SetCenter(L * (i + 0.5) / n_indicators, 0, radius * 1.5)
        sphere.SetThetaResolution(10)
        sphere.SetPhiResolution(10)
        sphere.Update()
        append.AddInputData(sphere.GetOutput())

    append.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(append.GetOutput())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1.0, 0.8, 0.0)  # Yellow/gold for heat source

    return actor


def vtk_3d_rod_visualization(
    x: np.ndarray,
    T: np.ndarray,
    T_left: float,
    T_right: float,
    S: float,
    title: str = "1D Heat Conduction with Fixed Ends",
):
    """
    Create an advanced 3D visualization of the rod with temperature distribution.

    Features:
    - 3D cylindrical rod with temperature color mapping
    - End caps showing fixed boundary temperatures
    - Heat source indicators
    - Boundary condition annotations

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param T_left: Left boundary temperature (float)
    :param T_right: Right boundary temperature (float)
    :param S: Heat generation rate (float)
    :param title: Window title (str)
    """
    L = x[-1] - x[0]
    radius = 0.08 * L

    # Create 3D cylinder rod
    rod_polydata, _ = create_cylinder_rod(x, T, radius=radius)

    # Create end caps
    left_cap = create_end_cap(x[0], T[0], radius=radius)
    right_cap = create_end_cap(x[-1], T[-1], radius=radius)

    # Combine rod and caps
    append_filter = vtk.vtkAppendPolyData()
    append_filter.AddInputData(rod_polydata)
    append_filter.AddInputData(left_cap)
    append_filter.AddInputData(right_cap)
    append_filter.Update()

    # Color map
    color_map = create_color_map(T)

    # Rod mapper and actor
    rod_mapper = vtk.vtkPolyDataMapper()
    rod_mapper.SetInputData(append_filter.GetOutput())
    rod_mapper.SetLookupTable(color_map)
    rod_mapper.SetScalarModeToUsePointData()
    rod_mapper.SetScalarRange(np.min(T), np.max(T))

    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(rod_mapper)
    rod_actor.GetProperty().SetInterpolationToPhong()
    rod_actor.GetProperty().SetSpecular(0.3)
    rod_actor.GetProperty().SetSpecularPower(20)

    # Create heat source indicators
    heat_source_actor = create_heat_source_indicator(L, radius)

    # Create boundary annotations
    annotation_offset = radius * 3
    left_text = f"Fixed BC\nT = {T_left}°C"
    right_text = f"Fixed BC\nT = {T_right}°C"
    source_text = f"Heat Generation\nS = {S:.0e} W/m³"

    left_annotation = create_boundary_annotation(
        left_text, (x[0] - 0.12 * L, 0, annotation_offset), color=(0.3, 0.6, 1.0)
    )
    right_annotation = create_boundary_annotation(
        right_text, (x[-1] + 0.12 * L, 0, annotation_offset), color=(1.0, 0.4, 0.4)
    )
    source_annotation = create_boundary_annotation(
        source_text, (L / 2, 0, radius * 2.5), color=(1.0, 0.8, 0.0)
    )

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

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(rod_actor)
    renderer.AddActor(heat_source_actor)
    renderer.AddActor(left_annotation)
    renderer.AddActor(right_annotation)
    renderer.AddActor(source_annotation)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)

    # Gradient background
    renderer.SetBackground(0.1, 0.1, 0.2)
    renderer.SetBackground2(0.02, 0.02, 0.05)
    renderer.GradientBackgroundOn()

    # Add lighting
    light = vtk.vtkLight()
    light.SetPosition(L / 2, L, L)
    light.SetFocalPoint(L / 2, 0, 0)
    light.SetIntensity(1.0)
    renderer.AddLight(light)

    # Render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1200, 700)
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

    # Camera position
    camera = renderer.GetActiveCamera()
    camera.SetPosition(L / 2, -L * 0.8, L * 0.6)
    camera.SetFocalPoint(L / 2, 0, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    camera.Zoom(1.2)

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    """
    Main function demonstrating 1D heat conduction with fixed-end temperatures.

    Part 1 of the Heat Transfer Series: Dirichlet (fixed temperature) BCs
    """
    # Physical parameters
    L = 1.0  # Rod length (m)
    n = 100  # Number of divisions
    k = 50.0  # Thermal conductivity (W/m·K)
    S = 1e5  # Heat generation rate (W/m³)
    T_left = 100.0  # Left boundary temperature (°C)
    T_right = 200.0  # Right boundary temperature (°C)

    # Create position array
    x = np.linspace(0, L, n + 1)

    # Setup and solve the problem
    diagonals, b = setup_problem(n, L, k, S, T_left, T_right)
    T_numerical = tdma_solve(*diagonals, b)

    # Compute analytical solution
    T_analytical = analytical_solution(x, L, k, S, T_left, T_right)

    # Plot the results
    plot_temperature_distribution(x, T_numerical, T_analytical)

    # 3D VTK visualization
    vtk_3d_rod_visualization(x, T_numerical, T_left, T_right, S)

    # Print results and energy balance
    print("=" * 65)
    print("Part 1: 1D Heat Conduction with Fixed-End Temperatures")
    print("=" * 65)
    print(f"Rod length: {L} m, Grid points: {n + 1}")
    print(f"Thermal conductivity: {k} W/m·K")
    print(f"Heat generation: S = {S:.2e} W/m³")
    print(f"Left BC: T = {T_left}°C (fixed)")
    print(f"Right BC: T = {T_right}°C (fixed)")
    print("-" * 65)
    print(f"Maximum temperature: {np.max(T_numerical):.2f}°C")
    print(f"Location of max temp: x = {x[np.argmax(T_numerical)]:.3f} m")
    print("-" * 65)

    # Energy balance verification
    dx = L / n
    left_flux = -k * (T_numerical[1] - T_numerical[0]) / dx
    right_flux = -k * (T_numerical[-1] - T_numerical[-2]) / dx
    total_flux_out = -left_flux + right_flux  # Heat leaving the rod
    total_generation = S * L  # Total heat generated

    print("Energy Balance Verification:")
    print(f"  Heat flux at x=0: {left_flux:.2f} W/m² (into rod)")
    print(f"  Heat flux at x=L: {right_flux:.2f} W/m² (out of rod)")
    print(f"  Total heat out: {total_flux_out:.2f} W/m²")
    print(f"  Total generated: {total_generation:.2f} W/m²")
    print(f"  Balance error: {abs(total_flux_out - total_generation):.4f} W/m²")
    print("-" * 65)
    max_error = np.max(np.abs(T_numerical - T_analytical))
    print(f"Maximum numerical error: {max_error:.6f}°C")
    print("=" * 65)


if __name__ == "__main__":
    main()
