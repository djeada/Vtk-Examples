"""
This script provides a numerical solution to two-dimensional heat transfer problems in a rectangular domain using the Gauss-Seidel method with the option of Line by Line or Point by Point approaches. It is particularly useful for simulations where heat source distribution and boundary conditions are known.

Workflow:
1. Initialization of grid size, material properties, heat source term, and boundary temperatures.
2. Iterative solution using the Gauss-Seidel method (either Line by Line or Point by Point) to update the temperature distribution based on the discretized heat equation.
3. Applying boundary conditions at each iteration and computing the residual to check for convergence.
4. Visual representation of the temperature field using a contour plot.

Physics:
The script models steady-state heat transfer in a two-dimensional domain. It accounts for internal heat generation and heat conduction according to Fourier's law, along with specified boundary conditions.

Mathematics:
The heat conduction equation is discretized using the finite difference method, creating a grid of equations that can be solved iteratively. The Gauss-Seidel method is used for this iterative process, with the flexibility of using different approaches for computational efficiency.
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


def two_dimensional_heat_transfer_solver(
    n: int, l: float, k: float, S: float, Ts: float, Rtol: float
) -> tuple:
    """
    Solves the two-dimensional heat transfer problem.

    :param n: Number of divisions in each direction
    :param l: Length of the domain
    :param k: Thermal conductivity
    :param S: Source term
    :param Ts: Boundary temperature
    :param Rtol: Residual tolerance for convergence
    :return: Tuple of (x, y, T) where x and y are coordinate vectors, and T is the temperature field
    """
    dx = dy = l / n
    x = np.arange(0, l + dx, dx)
    y = np.arange(0, l + dy, dy)
    T = 5 * Ts * np.ones((n + 1, n + 1))
    Rmax = 1
    iter_num = 0

    while Rmax > Rtol:
        Rmax = 0

        # Point by Point Gauss Seidel Method
        for i in range(1, n):
            for j in range(1, n):
                # North, West, South, East Points
                Aw = Ae = k * dy / dx
                An = As = k * dx / dy

                # Ap is the coefficient of T[i, j]
                Ap = Aw + Ae + An + As

                # Generation term
                b_ij = S * dx * dy

                # Calculate new temperature at (i, j)
                T_new = (
                    Aw * T[i, j - 1]
                    + Ae * T[i, j + 1]
                    + An * T[i + 1, j]
                    + As * T[i - 1, j]
                    + b_ij
                ) / Ap

                # Calculate residual and update temperature
                R = abs(T_new - T[i, j])
                Rmax = max(R, Rmax)
                T[i, j] = T_new

        # Boundary Conditions
        # Top and Right, temperature given
        T[:, n] = Ts
        T[n, :] = Ts

        # Left and Bottom, zero flux (symmetry)
        T[:, 0] = T[:, 1]
        T[0, :] = T[1, :]

        iter_num += 1
        print(f"Iteration: {iter_num}, Residual: {Rmax}")

    return x, y, T


def plot_temperature_profile(X: np.ndarray, Y: np.ndarray, T: np.ndarray, Ts: float):
    """
    Plots the temperature profile.

    :param X: X-coordinate grid
    :param Y: Y-coordinate grid
    :param T: Temperature field
    :param Ts: Boundary temperature
    """
    plt.contourf(X, Y, T, np.arange(Ts, 700, 5))
    plt.colorbar(ticks=np.arange(Ts, 700, 100))
    plt.axis([0, 1, 0, 1])
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Two-Dimensional Heat Transfer Solver")
    plt.show()


def create_points(x: np.ndarray, y: np.ndarray, T: np.ndarray) -> tuple:
    """
    Creates points and temperature values for vtk visualization.

    :param x: X-coordinate vector
    :param y: Y-coordinate vector
    :param T: Temperature field
    :return: Tuple of (points, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    n = len(x)
    for i in range(n):
        for j in range(n):
            points.InsertNextPoint(x[i], y[j], 0)
            temperature_values.InsertNextValue(T[i, j])

    return points, temperature_values


def create_poly_data(
    n: int, points: vtk.vtkPoints, temperature_values: vtk.vtkFloatArray
) -> vtk.vtkPolyData:
    """
    Creates vtkPolyData for visualization.

    :param n: Number of divisions in each direction
    :param points: vtkPoints
    :param temperature_values: vtkFloatArray
    :return: vtkPolyData
    """
    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    cells = vtk.vtkCellArray()
    for i in range(n - 1):
        for j in range(n - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, i * n + j)
            quad.GetPointIds().SetId(1, (i + 1) * n + j)
            quad.GetPointIds().SetId(2, (i + 1) * n + (j + 1))
            quad.GetPointIds().SetId(3, i * n + (j + 1))
            cells.InsertNextCell(quad)

    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(temperature_values)

    return polyData


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


def vtk_surface_plot(x: np.ndarray, y: np.ndarray, T: np.ndarray):
    """
    Visualizes the temperature field using VTK.

    :param x: X-coordinate vector
    :param y: Y-coordinate vector
    :param T: Temperature field
    """
    points, temperature_values = create_points(x, y, T)
    polyData = create_poly_data(len(x), points, temperature_values)
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
    # Parameters (can be modified as needed)
    n = 20
    l = 1.0
    k = 0.5
    S = 1000
    Ts = 100
    Rtol = 1e-3

    # Solve and plot
    x, y, T = two_dimensional_heat_transfer_solver(n, l, k, S, Ts, Rtol)
    X, Y = np.meshgrid(x, y)

    # Plot temperature profile
    plot_temperature_profile(X, Y, T, Ts)
    vtk_surface_plot(x, y, T)


if __name__ == "__main__":
    main()
