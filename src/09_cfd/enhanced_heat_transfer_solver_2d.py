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

import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.field_visualization import field_points, quad_field
from src.common.heat_transfer import solve_quarter_conduction, validate_field


def two_dimensional_heat_transfer_solver(
    n: int,
    L: float,
    k: float,
    S: float,
    T_boundary: float,
    Rtol: float,
    max_iterations: int = 50000,
) -> tuple:
    """Solve a quarter of a symmetric heated square, per unit depth.

    Left/bottom faces are insulated; top/right nodes have T_boundary. Boundary
    control volumes include source generation and conserve their normal fluxes.
    The residual is measured in kelvin after each complete sweep. Failure to
    reach Rtol raises ConvergenceError. Returns x, y, T[y,x], iterations.
    """
    return solve_quarter_conduction(n, L, k, S, T_boundary, Rtol, max_iterations)


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
    """Attach T[y,x] to its physical coordinates without displacement."""
    return field_points(x, y, T)


def create_poly_data(nx: int, ny: int, points, temperature_values):
    """Connect the scalar grid with consistently oriented quads."""
    return quad_field(nx, ny, points, temperature_values)


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
    x, y, T = validate_field(x, y, T)
    T_min, T_max = np.min(T), np.max(T)
    points, temperature_values = field_points(x, y, T, height=L * 0.5)
    polyData = create_poly_data(len(x), len(y), points, temperature_values)
    color_map = create_color_map(T)

    # Mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.InterpolateScalarsBeforeMappingOn()
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetScalarRange(T_min, T_max)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)
    surface_actor.GetProperty().LightingOff()
    surface_actor.GetProperty().SetInterpolationToPhong()
    surface_actor.GetProperty().SetSpecular(0.3)
    surface_actor.GetProperty().SetSpecularPower(20)

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.GetTitleTextProperty().SetFontSize(16)
    scalar_bar.GetLabelTextProperty().SetFontSize(14)
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
    info_text = f"S = {S:.0e} W/m³\nT_boundary = {T_boundary}°C\nT_max = {T_max:.1f}°C\nHeight: normalized temperature (display only)"
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
        label.AddPosition(0, 0, L * 0.5)
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
    axes.SetZAxisLabelText("scalar (scaled)")
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0, 0.82, 0.15, 1)
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
    print(f"Domain: {L}m × {L}m, Grid: {n + 1} × {n + 1}")
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
    print(f"Temperature at center: {T[n // 2, n // 2]:.2f}°C")
    print("=" * 70)

    # Plot temperature profile
    plot_temperature_profile(X, Y, T, S)

    # 3D VTK visualization
    vtk_3d_surface_plot(x, y, T, L, S, T_boundary)


if __name__ == "__main__":
    main()
