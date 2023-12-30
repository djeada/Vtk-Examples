"""
This module provides an advanced solver for 1D heat transfer problems, incorporating both conduction and convection phenomena. It is designed to handle various boundary conditions and supports different discretization schemes, including Upwind Difference Scheme (UDS) and Central Difference Scheme (CDS). The solver is particularly useful for studying the effects of convective heat transfer in one-dimensional domains.

Workflow:
1. Initialization: Set up the computational grid, material properties (like diffusivity), and flow characteristics.
2. Matrix Assembly: Depending on the chosen scheme (UDS or CDS), the module constructs the tridiagonal matrix representing the discretized form of the heat convection-diffusion equation.
3. Boundary Conditions: Apply specified boundary conditions at the ends of the domain.
4. Solver: Utilize the Tridiagonal Matrix Algorithm (TDMA) to efficiently solve the linear system.
5. Post-Processing: Generate numerical results and compare them with the analytical solution for validation.

Physics:
The module addresses the one-dimensional convection-diffusion equation, which combines the effects of heat conduction (diffusion) and convection. It models how a scalar quantity (like temperature) is transported through a medium under the influence of these two processes.

Mathematics:
- Discretization: Spatial discretization of the convection-diffusion equation using finite difference methods.
- TDMA: The Tridiagonal Matrix Algorithm, an efficient method for solving tridiagonal systems of linear equations, crucial for computational efficiency in 1D problems.
- Analytical Solution Comparison: Provides a method to compare numerical results with an analytical solution, offering a means for validation and accuracy assessment.
"""

import matplotlib.pyplot as plt
import numpy as np
import vtk


def tdma_solve(
    a: np.ndarray, ad: np.ndarray, au: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """
    Solves a tridiagonal matrix equation.

    :param a: Main diagonal
    :param ad: Sub-diagonal (lower)
    :param au: Super-diagonal (upper)
    :param b: Right-hand side vector
    :return: Solution vector
    """
    n = len(b)
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


def setup_matrix(n, dx, T, F, uds):
    """
    Sets up the matrix for the TDMA solver based on the chosen scheme (UDS or CDS).

    Parameters:
    n (int): Number of grid points.
    dx (float): Grid spacing.
    T (float): Diffusivity.
    F (float): Convection term.
    uds (bool): Flag for Upwind Difference Scheme (True) or Central Difference Scheme (False).

    Returns:
    tuple: Matrices for TDMA solver.
    """
    A = np.zeros(n + 1)
    Au = np.zeros(n + 1)
    Ad = np.zeros(n + 1)

    if uds:
        A.fill(2 * (T / dx) + F)
        Au.fill(-T / dx)
        Ad.fill(-T / dx - F)
    else:
        A.fill(2 * (T / dx))
        Au.fill(-T / dx + F / 2)
        Ad.fill(-T / dx - F / 2)

    return A, Ad, Au


def plot_solution(x, phi, xt, y):
    """
    Plots the numerical and analytical solutions.

    Parameters:
    x (np.ndarray): Grid points for numerical solution.
    phi (np.ndarray): Numerical solution.
    xt (np.ndarray): Grid points for analytical solution.
    y (np.ndarray): Analytical solution.
    """
    plt.plot(xt, y, color="red", label="Analytical Solution")
    plt.plot(x, phi, color="blue", linestyle="--", label="Numerical Solution")
    plt.xlabel("x (m)")
    plt.ylabel("phi")
    plt.axis([0, 1, 0, 1])
    plt.legend()
    plt.show()


def create_points(x: np.ndarray, phi: np.ndarray) -> tuple:
    """
    Creates points and phi values for vtk visualization.

    :param x: X-coordinate vector
    :param phi: Phi field (can represent temperature, concentration, etc.)
    :return: Tuple of (points, phi_values) for vtk
    """
    points = vtk.vtkPoints()
    phi_values = vtk.vtkFloatArray()
    phi_values.SetName("Phi")

    for i in range(len(x)):
        points.InsertNextPoint(x[i], 0, 0)
        phi_values.InsertNextValue(phi[i])

    return points, phi_values


def create_polyline(n: int) -> vtk.vtkPolyLine:
    """
    Creates a vtkPolyLine for visualization.

    :param n: Number of points in the polyline
    :return: vtkPolyLine
    """
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(n)
    for i in range(n):
        polyline.GetPointIds().SetId(i, i)
    return polyline


def create_color_map(phi: np.ndarray) -> vtk.vtkLookupTable:
    """
    Creates a color map for vtk visualization.

    :param phi: Phi field
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(phi), np.max(phi)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


def vtk_line_plot(x: np.ndarray, phi: np.ndarray):
    """
    Visualizes the phi field using VTK.

    :param x: X-coordinate vector
    :param phi: Phi field
    """
    points, phi_values = create_points(x, phi)
    polyline = create_polyline(len(x))

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)
    polyData.GetPointData().AddArray(phi_values)
    polyData.GetPointData().SetActiveScalars("Phi")

    color_map = create_color_map(phi)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(phi), np.max(phi))

    # Rod actor setup
    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(mapper)

    # Set the width of the rod
    rod_actor.GetProperty().SetLineWidth(5)

    # Scalar bar (legend)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Phi")
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

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    n = 400
    l = 1.0
    dx = l / n
    rho = 1.0
    T = 0.05
    u = 2
    F = rho * u
    phi_0, phi_l = 0.0, 1.0
    x = np.arange(0, l + dx, dx)
    uds = True

    A, Ad, Au = setup_matrix(n, dx, T, F, uds)

    b = np.zeros(n + 1)
    b[0], b[-1] = phi_0, phi_l
    A[0], A[-1] = 1, 1
    Au[0], Ad[-1] = 0, 0

    phi = tdma_solve(A, Ad, Au, b)

    P = rho * u * l / T
    xt = np.arange(0, l + l / 200, 0.1 / 200)
    y = [
        phi_0 + (phi_l - phi_0) * (np.exp(P * i / l) - 1) / (np.exp(P) - 1) for i in xt
    ]

    plot_solution(x, phi, xt, y)
    vtk_line_plot(x, phi)


if __name__ == "__main__":
    main()
