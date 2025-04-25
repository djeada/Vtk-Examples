"""
This script simulates the steady-state temperature distribution along a one-dimensional rod using the Tri-Diagonal Matrix Algorithm (TDMA), also known as the Thomas Algorithm. It is designed to solve problems involving heat conduction in materials with temperature-dependent conductivity.

Workflow:
1. Initialization of material properties, grid size, boundary conditions, and initial temperature distribution.
2. Iterative solution of the discretized heat conduction equation using TDMA for a specified convergence criterion.
3. Visualization of the temperature distribution along the rod.

Physics:
The script models heat conduction in a rod, considering the linear heat conduction equation with the Fourier's law. It incorporates boundary conditions like constant temperature or convective heat transfer at the ends of the rod.

Mathematics:
The heat conduction equation is discretized using the finite difference method, resulting in a linear system of equations with a tridiagonal matrix. The TDMA efficiently solves this system iteratively until the solution converges below a predefined residual tolerance.
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


def k(T: np.ndarray, kc: float, beta: float) -> np.ndarray:
    """
    Calculates temperature-dependent conductivity.

    :param T: Temperature
    :param kc: Conductivity coefficient
    :param beta: Temperature coefficient
    :return: Conductivity
    """
    return kc * (1 + beta * T)


def solve_heat_transfer_equation(
    n: int,
    l: float,
    kc: float,
    beta: float,
    h_r: float,
    h_l: float,
    T_l: float,
    T_r: float,
    Rtol: float,
) -> tuple:
    """
    Solves the enhanced heat transfer problem.

    :param n: Number of divisions in the rod
    :param l: Length of the rod
    :param kc: Conductivity coefficient
    :param beta: Temperature coefficient
    :param h_r: Right-side heat transfer coefficient
    :param h_l: Left-side heat transfer coefficient
    :param T_l: Left-side boundary temperature
    :param T_r: Right-side boundary temperature
    :param Rtol: Residual tolerance for convergence
    :return: Tuple of (x, T) where x is the coordinate vector, and T is the temperature field
    """
    dx = l / n
    x = np.arange(0, l + dx, dx)
    T = 300 * np.ones(n + 1)
    A = np.zeros(n + 1)
    Au = np.ones(n + 1)
    Ad = np.ones(n + 1)
    Rmax = 1

    while Rmax > Rtol:
        Rmax = 0
        ks = k(T, kc, beta)

        for i in range(1, n):
            kw = 2 * ks[i] * ks[i - 1] / (ks[i] + ks[i - 1])
            ke = 2 * ks[i] * ks[i + 1] / (ks[i] + ks[i + 1])
            A[i] = (kw + ke) / dx
            Ad[i] = -kw / dx
            Au[i] = -ke / dx
            R = abs(A[i] * T[i] + Ad[i] * T[i - 1] + Au[i] * T[i + 1])
            Rmax = max(R, Rmax)

        A[0] = ks[0] / dx + h_l
        A[n] = ks[n] / dx + h_r
        Au[0] = -ks[0] / dx
        Ad[n] = -ks[n] / dx

        b = np.zeros(n + 1)
        b[0] = h_l * T_l
        b[n] = h_r * T_r

        R = abs(A[n] * T[n] + Ad[n] * T[n - 1] - b[n])
        Rmax = max(R, Rmax)
        R = abs(A[0] * T[0] + Au[0] * T[1] - b[0])
        Rmax = max(R, Rmax)

        T = tdma_solve(A, Ad, Au, b)

    return x, T


def plot_temperature_profile(x: np.ndarray, T: np.ndarray):
    """
    Plots the temperature profile.

    :param x: X-coordinate vector
    :param T: Temperature field
    """
    plt.plot(x, T)
    plt.xlabel("x (m)")
    plt.ylabel("T (degrees Celsius)")
    plt.title("Temperature Profile")
    plt.show()


def create_points(x: np.ndarray, T: np.ndarray) -> tuple:
    """
    Creates points and temperature values for vtk visualization.

    :param x: X-coordinate vector
    :param T: Temperature field
    :return: Tuple of (points, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    for i in range(len(x)):
        points.InsertNextPoint(x[i], 0, 0)
        temperature_values.InsertNextValue(T[i])

    return points, temperature_values


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


def create_color_map(T: np.ndarray) -> vtk.vtkLookupTable:
    """
    Creates a color map for vtk visualization.

    :param T: Temperature field
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(T), np.max(T)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


def vtk_line_plot(x: np.ndarray, T: np.ndarray):
    """
    Visualizes the temperature field using VTK.

    :param x: X-coordinate vector
    :param T: Temperature field
    """
    points, temperature_values = create_points(x, T)
    polyline = create_polyline(len(x))

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

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    # Parameters (can be modified as needed)
    n = 10
    l = 0.1
    kc = 0.05
    beta = 0.001
    h_r = 5.0
    h_l = 30.0
    T_l = 300.0
    T_r = 30.0
    Rtol = 1e-6

    # Solve heat transfer equation
    x, T = solve_heat_transfer_equation(n, l, kc, beta, h_r, h_l, T_l, T_r, Rtol)

    # Plot temperature profile
    plot_temperature_profile(x, T)
    vtk_line_plot(x, T)


if __name__ == "__main__":
    main()
