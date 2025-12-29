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


def create_points(x: np.ndarray, T: np.ndarray) -> tuple:
    """
    Create points and temperature values for vtk visualization.

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :return: Tuple of (points, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    for i in range(len(x)):
        # Points along the x-axis at y=0, z=0
        points.InsertNextPoint(x[i], 0, 0)
        # Insert temperature values in correct order
        temperature_values.InsertNextValue(T[i])

    return points, temperature_values


def create_polyline(x: np.ndarray) -> vtk.vtkPolyLine:
    """
    Create a vtkPolyLine for visualization.

    :param x: X-coordinate vector (numpy.ndarray)
    :return: vtkPolyLine
    """
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(len(x))
    for i in range(len(x)):
        polyline.GetPointIds().SetId(i, i)
    return polyline


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


def vtk_line_plot(
    x: np.ndarray,
    T: np.ndarray,
    title: str = "Temperature Distribution",
):
    """
    Visualize the temperature field using VTK.

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param title: Title for the visualization (str)
    """
    points, temperature_values = create_points(x, T)
    polyline = create_polyline(x)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)
    polyData.GetPointData().AddArray(temperature_values)
    polyData.GetPointData().SetActiveScalars("Temperature")

    color_map = create_color_map(T)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(T), np.max(T))

    # Rod actor setup
    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(mapper)

    # Set the width of the rod
    rod_actor.GetProperty().SetLineWidth(5)
    # Scalar bar (legend)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature")
    scalar_bar.SetNumberOfLabels(5)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(rod_actor)
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

    # Alternatively, adjust the camera's clipping range manually if needed
    # camera = renderer.GetActiveCamera()
    # camera.SetClippingRange(minRange, maxRange)

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
    vtk_line_plot(x, T_numerical)

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
