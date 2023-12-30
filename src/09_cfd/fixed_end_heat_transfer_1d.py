"""
1D Temperature Distribution in a Rod

This module simulates the steady-state temperature distribution along a one-dimensional rod using the finite difference method and solves the resulting linear system with the Tri-Diagonal Matrix Algorithm (TDMA), also known as the Thomas algorithm.


Workflow:
1. The rod is discretized into 'n' equal segments, defining a grid along its length.
2. For each grid point, the discretized heat equation is written, forming a tri-diagonal matrix system (A, Ad, Au) and a right-hand side vector 'b'.
3. Boundary conditions are applied at the ends of the rod.
4. The TDMA solver, implemented in the function 'tdma_solve', efficiently solves this tri-diagonal system to find the temperature at each grid point.
5. The temperature distribution is then plotted, and an energy balance calculation is performed to verify the solution's accuracy.

Physics and Mathematics:
- The problem models heat conduction along a rod of length 'l' with constant thermal conductivity 'k'.
- The rod generates heat uniformly with a source term 'S'.
- Boundary conditions are defined as constant temperatures at both ends of the rod.
- The heat conduction in the rod is described by the steady-state heat equation:
  d/dx(k * dT/dx) = -S, where T is the temperature.
- This differential equation is discretized using the finite difference method, leading to a linear system of equations.
"""
import numpy as np
import vtk
from matplotlib import pyplot as plt


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


def setup_problem(n: int, l: float, k: float, S: float) -> tuple:
    """
    Sets up the problem matrix and vectors based on physical parameters.

    :param n: Number of divisions in the domain.
    :param l: Length of the domain.
    :param k: Thermal conductivity.
    :param S: Source term (heat generation rate).
    :return: Diagonals of the matrix (A, Ad, Au) and vector b.
    """
    dx = l / n
    A = 2 * (k / dx) * np.ones(n)
    A[0], A[-1] = 1, 1  # Boundary conditions
    Au = -(k / dx) * np.ones(n)
    Au[0] = 0
    Ad = -(k / dx) * np.ones(n)
    Ad[-1] = 0
    b = (S * dx) * np.ones(n)
    b[0], b[-1] = 100, 200  # Boundary values
    return (A, Ad, Au), b


def plot_temperature_distribution(x: np.ndarray, T: np.ndarray):
    """
    Plots the temperature distribution.

    :param x: Array of positions.
    :param T: Array of temperature values.
    """
    plt.plot(x, T, color="blue")
    plt.xlabel("Position")
    plt.ylabel("Temperature")
    plt.title("Temperature Distribution")
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

    # Alternatively, adjust the camera's clipping range manually if needed
    # camera = renderer.GetActiveCamera()
    # camera.SetClippingRange(minRange, maxRange)

    renderWindow.Render()
    renderWindowInteractor.Start()


def main():
    # Problem parameters
    n = 100
    l = 0.02
    k = 0.5
    S = 1e6

    # Setup and solve the problem
    diagonals, b = setup_problem(n, l, k, S)
    x = np.linspace(0, l, n)
    T = tdma_solve(*diagonals, b)

    # Plot the results
    plot_temperature_distribution(x, T)
    vtk_line_plot(x, T)

    # Energy balance calculation
    dx = l / n
    left_flux = -(k / dx) * (T[0] - T[1])
    right_flux = -(k / dx) * (T[-1] - T[-2])
    total_flux = left_flux + right_flux
    source = S * l

    print(f"Left Flux: {left_flux}, Right Flux: {right_flux}")
    print(f"Total Flux: {total_flux}, Source: {source}")


if __name__ == "__main__":
    main()
