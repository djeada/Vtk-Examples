"""
Heat Transfer Modeling and TDMA Solver

This module provides functions to model one-dimensional steady-state heat transfer
in a rod and solve the resulting linear system using the Tri-Diagonal Matrix Algorithm (TDMA).
The rod is assumed to be homogeneous with constant thermal conductivity.

Workflow:
1. Define the problem's physical parameters (rod length, thermal conductivity, heat transfer coefficients, etc.).
2. Discretize the rod into a finite number of nodes and calculate the grid spacing.
3. Apply boundary conditions to the ends of the rod.
4. Formulate the linear system of equations based on the finite difference method.
5. Solve the linear system using TDMA.
6. Compare the numerical solution with an analytical solution (if available).
7. Visualize the temperature distribution along the rod using VTK.

Physics:
A rod of length L with thermal conductivity k has convective heat transfer at both ends.
The left end (x=0) is exposed to fluid at temperature T_l with heat transfer coefficient h_l.
The right end (x=L) is exposed to fluid at temperature T_r with heat transfer coefficient h_r.

Mathematics:

The heat transfer in the rod is modeled using the steady-state heat conduction equation:
    d^2T/dx^2 = 0
where T is the temperature and x is the position along the rod.

For interior nodes, the discretized equation using central difference is:
    -T[i-1] + 2*T[i] - T[i+1] = 0

Boundary Conditions (convective at both ends):
At x = 0 (left end), heat balance:
    -k * (T[1] - T[0])/dx = h_l * (T[0] - T_l)
    Rearranged: (k/dx + h_l)*T[0] - (k/dx)*T[1] = h_l*T_l

At x = L (right end), heat balance:
    k * (T[n] - T[n-1])/dx = h_r * (T[n] - T_r)
    Rearranged: -(k/dx)*T[n-1] + (k/dx + h_r)*T[n] = h_r*T_r

Analytical Solution:
Since d^2T/dx^2 = 0, T(x) = a*x + c (linear profile).
Using thermal resistance approach:
    R_total = 1/h_l + L/k + 1/h_r  (total thermal resistance)
    q = (T_r - T_l) / R_total       (heat flux)
    T(0) = T_l + q/h_l              (temperature at x=0)
    T(L) = T_r - q/h_r              (temperature at x=L)
    a = (T(L) - T(0)) / L           (slope)

TDMA Solver:
The TDMA solver efficiently solves tridiagonal systems of equations.
The system has the form: ad[i]*T[i-1] + a[i]*T[i] + au[i]*T[i+1] = b[i]
where 'a' is the main diagonal, 'ad' is the sub-diagonal (lower), 
'au' is the super-diagonal (upper), and 'b' is the right-hand side.

The algorithm involves two steps:
1. Forward elimination:
   For i = 1 to n-1:
       f = ad[i] / a[i-1]
       a[i] = a[i] - f * au[i-1]
       b[i] = b[i] - f * b[i-1]
2. Back substitution:
   T[n-1] = b[n-1] / a[n-1]
   For i = n-2 to 0:
       T[i] = (b[i] - au[i] * T[i+1]) / a[i]
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


def initialize_matrices(
    n: int, dx: float, k: float, h_l: float, h_r: float, T_l: float, T_r: float
) -> tuple:
    """
    Initialize the matrices and vectors for the TDMA solver.

    All equations are normalized with k/dx to ensure consistent scaling:
    - Interior nodes (i = 1 to n-1): -(k/dx)*T[i-1] + 2*(k/dx)*T[i] - (k/dx)*T[i+1] = 0
    - Left boundary (i = 0): (k/dx + h_l)*T[0] - (k/dx)*T[1] = h_l*T_l
    - Right boundary (i = n): -(k/dx)*T[n-1] + (k/dx + h_r)*T[n] = h_r*T_r

    :param n: Number of divisions in the domain (int)
    :param dx: Grid spacing (float)
    :param k: Thermal conductivity (float)
    :param h_l: Heat transfer coefficient at the left boundary (float)
    :param h_r: Heat transfer coefficient at the right boundary (float)
    :param T_l: Surrounding temperature at the left boundary (float)
    :param T_r: Surrounding temperature at the right boundary (float)
    :return: Tuple of (A, Ad, Au, b) representing the main, lower, and upper diagonals, and the right-hand side vector.
    """
    # Number of nodes is n + 1 (from 0 to n)
    num_nodes = n + 1

    # Conduction coefficient
    cond = k / dx

    # Initialize diagonals for interior nodes (using k/dx scaling)
    A = 2.0 * cond * np.ones(num_nodes)
    Au = -cond * np.ones(num_nodes)
    Ad = -cond * np.ones(num_nodes)

    # Apply left boundary condition: (k/dx + h_l)*T[0] - (k/dx)*T[1] = h_l*T_l
    A[0] = cond + h_l
    Au[0] = -cond
    Ad[0] = 0.0  # No lower diagonal at first node

    # Apply right boundary condition: -(k/dx)*T[n-1] + (k/dx + h_r)*T[n] = h_r*T_r
    A[-1] = cond + h_r
    Ad[-1] = -cond
    Au[-1] = 0.0  # No upper diagonal at last node

    # Initialize right-hand side vector
    b = np.zeros(num_nodes)
    b[0] = h_l * T_l  # Left boundary
    b[-1] = h_r * T_r  # Right boundary

    return A, Ad, Au, b


def analytical_solution(
    x: np.ndarray, L: float, k: float, h_l: float, h_r: float, T_l: float, T_r: float
) -> np.ndarray:
    """
    Compute the analytical solution for 1D steady-state heat conduction with
    convective boundary conditions at both ends.

    Since d^2T/dx^2 = 0, the temperature profile is linear: T(x) = a*x + c

    Using thermal resistance approach:
    - Total resistance: R_total = 1/h_l + L/k + 1/h_r
    - Heat flux: q = (T_r - T_l) / R_total
    - Temperature at x=0: T(0) = T_l + q/h_l
    - Temperature at x=L: T(L) = T_r - q/h_r

    :param x: Position array (numpy.ndarray)
    :param L: Length of the rod (float)
    :param k: Thermal conductivity (float)
    :param h_l: Heat transfer coefficient at left end (float)
    :param h_r: Heat transfer coefficient at right end (float)
    :param T_l: Surrounding temperature at left end (float)
    :param T_r: Surrounding temperature at right end (float)
    :return: Temperature array (numpy.ndarray)
    """
    # Total thermal resistance
    R_total = 1 / h_l + L / k + 1 / h_r

    # Heat flux (positive when flowing from left to right, i.e., when T_r > T_l)
    q = (T_r - T_l) / R_total

    # Temperature at boundaries
    T_0 = T_l + q / h_l
    T_L = T_r - q / h_r

    # Linear temperature profile
    a = (T_L - T_0) / L
    return a * x + T_0


def plot_results(x: np.ndarray, T_numerical: np.ndarray, T_analytical: np.ndarray):
    """
    Plot the numerical and analytical solutions.

    :param x: Array of positions (numpy.ndarray)
    :param T_numerical: Numerical solution (numpy.ndarray)
    :param T_analytical: Analytical solution (numpy.ndarray)
    """
    plt.plot(x, T_numerical, label="Numerical Solution")
    plt.plot(x, T_analytical, label="Analytical Solution", linestyle="--")
    plt.xlabel("x (m)")
    plt.ylabel("Temperature (°C)")
    plt.legend()
    plt.show()


def create_color_map(T: np.ndarray) -> vtk.vtkLookupTable:
    """
    Create a color map for vtk visualization.

    :param T: Temperature field (numpy.ndarray)
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(T), np.max(T)
    color_map.SetRange(min_val, max_val)

    # Blue to red color range (cold to hot)
    color_map.SetHueRange(0.667, 0)

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
            # Current ring, current point
            p0 = i * resolution + j
            # Current ring, next point
            p1 = i * resolution + (j + 1) % resolution
            # Next ring, next point
            p2 = (i + 1) * resolution + (j + 1) % resolution
            # Next ring, current point
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

    # Transform to position at x_pos
    transform = vtk.vtkTransform()
    transform.Translate(x_pos, 0, 0)
    transform.RotateY(90)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(disk.GetOutput())
    transformFilter.Update()

    # Add temperature scalar
    output = transformFilter.GetOutput()
    temp_array = vtk.vtkFloatArray()
    temp_array.SetName("Temperature")
    for _ in range(output.GetNumberOfPoints()):
        temp_array.InsertNextValue(T_val)
    output.GetPointData().SetScalars(temp_array)

    return output


def create_heat_flow_arrow(
    start_pos: tuple, direction: tuple, scale: float = 0.15
) -> vtk.vtkActor:
    """
    Create an arrow to represent heat flow direction.

    :param start_pos: Starting position (x, y, z)
    :param direction: Direction vector (dx, dy, dz)
    :param scale: Scale of the arrow (float)
    :return: vtkActor for the arrow
    """
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(20)
    arrow.SetShaftResolution(20)

    # Normalize direction
    length = np.sqrt(sum(d**2 for d in direction))
    if length == 0:
        length = 1

    transform = vtk.vtkTransform()
    transform.Translate(start_pos)

    # Calculate rotation to align arrow with direction
    # Default arrow points in +X direction
    if direction[0] < 0:
        transform.RotateZ(180)

    transform.Scale(scale, scale * 0.3, scale * 0.3)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputConnection(arrow.GetOutputPort())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transformFilter.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1.0, 0.5, 0.0)  # Orange color for heat flow

    return actor


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


def vtk_3d_rod_visualization(
    x: np.ndarray,
    T: np.ndarray,
    T_l: float,
    T_r: float,
    h_l: float,
    h_r: float,
    title: str = "1D Heat Conduction with Convective BCs",
):
    """
    Create an advanced 3D visualization of the rod with temperature distribution.

    Features:
    - 3D cylindrical rod with temperature color mapping
    - End caps for the rod
    - Heat flow arrows showing direction
    - Boundary condition annotations
    - Improved lighting and background

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param T_l: Left surrounding temperature (float)
    :param T_r: Right surrounding temperature (float)
    :param h_l: Left heat transfer coefficient (float)
    :param h_r: Right heat transfer coefficient (float)
    :param title: Window title (str)
    """
    L = x[-1] - x[0]
    radius = 0.08 * L  # Rod radius proportional to length

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

    # Create heat flow arrow (pointing in direction of heat flow)
    heat_flow_direction = 1 if T_r > T_l else -1
    arrow_y_offset = radius * 2.5
    arrow = create_heat_flow_arrow(
        start_pos=(L / 2 - 0.1 * L * heat_flow_direction, arrow_y_offset, 0),
        direction=(heat_flow_direction, 0, 0),
        scale=0.2 * L,
    )

    # Create boundary annotations
    annotation_offset = radius * 3
    left_text = f"Fluid\nT∞ = {T_l}°C\nh = {h_l} W/m²K"
    right_text = f"Fluid\nT∞ = {T_r}°C\nh = {h_r} W/m²K"

    left_annotation = create_boundary_annotation(
        left_text, (x[0] - 0.15 * L, 0, annotation_offset), color=(0.3, 0.6, 1.0)
    )
    right_annotation = create_boundary_annotation(
        right_text, (x[-1] + 0.15 * L, 0, annotation_offset), color=(1.0, 0.4, 0.4)
    )

    # Heat flow label
    flow_label = create_boundary_annotation(
        "Heat Flow →" if heat_flow_direction > 0 else "← Heat Flow",
        (L / 2, arrow_y_offset + radius, 0),
        color=(1.0, 0.6, 0.2),
    )

    # Scalar bar (legend)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.08)
    scalar_bar.SetHeight(0.4)
    scalar_bar.SetPosition(0.9, 0.3)
    scalar_bar.GetTitleTextProperty().SetFontSize(12)
    scalar_bar.GetLabelTextProperty().SetFontSize(10)

    # Title annotation
    title_actor = vtk.vtkTextActor()
    title_actor.SetInput(title)
    title_actor.GetTextProperty().SetFontSize(20)
    title_actor.GetTextProperty().SetColor(1, 1, 1)
    title_actor.GetTextProperty().SetBold(True)
    title_actor.SetPosition(10, 10)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(rod_actor)
    renderer.AddActor(arrow)
    renderer.AddActor(left_annotation)
    renderer.AddActor(right_annotation)
    renderer.AddActor(flow_label)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)

    # Set gradient background (dark blue to black)
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

    # Interaction style
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

    # Set camera position for better initial view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(L / 2, -L * 0.8, L * 0.6)
    camera.SetFocalPoint(L / 2, 0, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    camera.Zoom(1.2)

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    # Problem setup: Rod of length L with convective BCs at both ends
    L = 1.0  # Rod length (m)
    n = 50  # Number of divisions
    dx = L / n  # Grid spacing
    k = 50.0  # Thermal conductivity (W/m·K)
    h_l = 100.0  # Heat transfer coefficient at left end (W/m²·K)
    h_r = 50.0  # Heat transfer coefficient at right end (W/m²·K)
    T_l = 100.0  # Surrounding temperature at left end (°C)
    T_r = 300.0  # Surrounding temperature at right end (°C)

    # Create position array
    x = np.linspace(0, L, n + 1)

    # Initialize matrices for TDMA solver
    A, Ad, Au, b = initialize_matrices(n, dx, k, h_l, h_r, T_l, T_r)

    # Solve the linear system
    T_numerical = tdma_solve(A, Ad, Au, b)

    # Compute analytical solution
    T_analytical = analytical_solution(x, L, k, h_l, h_r, T_l, T_r)

    # Plot results
    plot_results(x, T_numerical, T_analytical)

    # 3D VTK visualization with boundary condition annotations
    vtk_3d_rod_visualization(x, T_numerical, T_l, T_r, h_l, h_r)

    # Print boundary values for verification
    print("=" * 60)
    print("1D Steady-State Heat Conduction with Convective BCs")
    print("=" * 60)
    print(f"Rod length: {L} m, Grid points: {n + 1}")
    print(f"Thermal conductivity: {k} W/m·K")
    print(f"Left BC:  h_l = {h_l} W/m²·K, T_l = {T_l} °C")
    print(f"Right BC: h_r = {h_r} W/m²·K, T_r = {T_r} °C")
    print("-" * 60)
    print(f"Numerical solution at x=0:  {T_numerical[0]:.4f} °C")
    print(f"Numerical solution at x=L:  {T_numerical[-1]:.4f} °C")
    print(f"Analytical solution at x=0: {T_analytical[0]:.4f} °C")
    print(f"Analytical solution at x=L: {T_analytical[-1]:.4f} °C")
    print("-" * 60)
    max_error = np.max(np.abs(T_numerical - T_analytical))
    print(f"Maximum error: {max_error:.6f} °C")
    print("=" * 60)


if __name__ == "__main__":
    main()
