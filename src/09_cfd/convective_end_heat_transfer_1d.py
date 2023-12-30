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

Mathematics:

The heat transfer in the rod is modeled using the steady-state heat conduction equation:
    -k * d^2T/dx^2 = 0
where k is the thermal conductivity, T is the temperature, and x is the position along the rod.

Discretization using the finite difference method leads to a system of linear equations:
    A[i-1]*T[i-1] + A[i]*T[i] + A[i+1]*T[i+1] = b[i]
for each interior node i, where A represents the coefficients matrix, T is the temperature vector, and b is the source term vector.

Boundary Conditions:
At x = 0 (left end):
    h_l*(T[0] - T_l) = -k*(T[1] - T[0])/dx
At x = L (right end):
    h_r*(T[n] - T_r) = k*(T[n] - T[n-1])/dx
where h_l, h_r are the convective heat transfer coefficients, T_l, T_r are the surrounding temperatures, and dx is the grid spacing.

TDMA Solver:
The TDMA solver efficiently solves tridiagonal systems of equations. The algorithm involves two steps:
1. Forward elimination:
   For i = 1 to n-1:
       w = c[i-1] / b[i-1]
       b[i] = b[i] - w * a[i]
       d[i] = d[i] - w * d[i-1]
2. Back substitution:
   For i = n-1 to 0:
       T[i] = (d[i] - a[i] * T[i+1]) / b[i]
"""
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import vtk


def tdma_solve(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """
    Solve a linear system with a tri-diagonal matrix using the Thomas algorithm.

    :param a: Lower diagonal elements (numpy.ndarray)
    :param b: Main diagonal elements (numpy.ndarray)
    :param c: Upper diagonal elements (numpy.ndarray)
    :param d: Right-hand side vector (numpy.ndarray)
    :return: Solution vector (numpy.ndarray)
    """
    n = len(d)
    w = np.zeros(n - 1, float)
    g = np.zeros(n, float)
    p = np.zeros(n, float)

    w[0] = c[0] / b[0]
    g[0] = d[0] / b[0]

    for i in range(1, n - 1):
        w[i] = c[i] / (b[i] - a[i - 1] * w[i - 1])
    for i in range(1, n):
        g[i] = (d[i] - a[i - 1] * g[i - 1]) / (b[i] - a[i - 1] * w[i - 1])

    p[n - 1] = g[n - 1]
    for i in range(n - 2, -1, -1):
        p[i] = g[i] - w[i] * p[i + 1]

    return p


def initialize_matrices(
    n: int, dx: float, k: float, h_l: float, h_r: float, T_l: float, T_r: float
) -> tuple:
    """
    Initialize the matrices and vectors for the TDMA solver.

    :param n: Number of divisions in the domain (int)
    :param dx: Grid spacing (float)
    :param k: Thermal conductivity (float)
    :param h_l: Heat transfer coefficient at the left boundary (float)
    :param h_r: Heat transfer coefficient at the right boundary (float)
    :param T_l: Temperature at the left boundary (float)
    :param T_r: Temperature at the right boundary (float)
    :return: Tuple of (A, Ad, Au, b) representing the main, lower, and upper diagonals, and the right-hand side vector.
    """
    # Initialize diagonals
    A = 2 * np.ones(n + 1)
    Au = -1 * np.ones(n + 1)
    Ad = -1 * np.ones(n + 1)

    # Apply boundary conditions
    A[0], A[-1] = k / dx + h_l, k / dx + h_r
    Au[0], Ad[-1] = -k / dx, -k / dx

    # Initialize right-hand side vector
    b = np.zeros(n + 1)
    b[0], b[-1] = h_l * T_l, h_r * T_r

    return A, Ad, Au, b


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


def create_points(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Create points and temperature values for vtk visualization.

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Temperature field (numpy.ndarray)
    :return: Tuple of (points, temperature_values) for vtk
    """
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")

    for i in range(len(x)):
        points.InsertNextPoint(x[i], y[i], 0)
        # Insert values in reverse order
        temperature_values.InsertNextValue(y[-(i + 1)])

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


def create_color_map(y: np.ndarray) -> vtk.vtkLookupTable:
    """
    Create a color map for vtk visualization.

    :param y: Temperature field (numpy.ndarray)
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(y), np.max(y)
    color_map.SetRange(min_val, max_val)

    # Reverse the color order
    color_map.SetHueRange(0.667, 0)  # Blue to red

    color_map.Build()
    return color_map


def apply_transformation(
    polyData: vtk.vtkPolyData, x_offset: float = 50, y_offset: float = 10
) -> vtk.vtkTransformPolyDataFilter:
    transform = vtk.vtkTransform()
    transform.Translate(x_offset, y_offset, 0)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(polyData)
    transformFilter.Update()

    return transformFilter


def vtk_line_plot(
    x: np.ndarray,
    y: np.ndarray,
    title: str = "Temperature Distribution",
    x_label: str = "Position",
    y_label: str = "Temperature",
):
    points, temperature_values = create_points(x, y)
    polyline = create_polyline(x)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)
    polyData.GetPointData().AddArray(temperature_values)
    polyData.GetPointData().SetActiveScalars("Temperature")

    transformFilter = apply_transformation(polyData)

    color_map = create_color_map(y)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transformFilter.GetOutputPort())
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(y), np.max(y))

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
    # Problem setup
    n = 10
    dx = 0.1 / n
    k = 0.1
    h_r, h_l = 5.0, 30.0
    T_r, T_l = 30.0, 300.0

    x = np.arange(0, 0.1 + dx, dx)
    A, Ad, Au, b = initialize_matrices(n, dx, k, h_l, h_r, T_l, T_r)

    # Solve the linear system
    T_numerical = tdma_solve(Ad, A, Au, b)

    # Analytical solution for comparison
    c2 = 292.70270270270271
    c1 = -300 * (300 - c2)
    T_analytical = c1 * x + c2

    # Plot results
    plot_results(x, T_numerical, T_analytical)
    vtk_line_plot(x, T_numerical)

    # Print boundary values
    print("Numerical solution at boundaries:", T_numerical[0], T_numerical[-1])
    print("Analytical solution at boundaries:", T_analytical[0], T_analytical[-1])


if __name__ == "__main__":
    main()
