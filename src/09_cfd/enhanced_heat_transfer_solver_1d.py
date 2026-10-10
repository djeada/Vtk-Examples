"""
1D Heat Transfer with Temperature-Dependent Conductivity

Part 3 of Heat Transfer Series: Non-Linear Material Properties

This module extends the basic 1D heat conduction to include temperature-dependent
thermal conductivity, making it a non-linear problem that requires iterative solution.

Series Overview:
- Part 1: Fixed temperature BCs - simplest case with Dirichlet conditions
- Part 2: Convective BCs - more realistic with convection at boundaries
- Part 3 (this file): Temperature-dependent conductivity - non-linear problem
- Part 4: 2D extension - full two-dimensional heat transfer

Workflow:
1. Initialize material properties including temperature coefficient for conductivity.
2. Set up convective boundary conditions at both ends.
3. Iteratively solve using TDMA with updated conductivity values.
4. Use harmonic mean for interface conductivity (for better accuracy).
5. Check residual for convergence.
6. Visualize showing the non-linear temperature profile.

Physics:
The thermal conductivity varies linearly with temperature:
    k(T) = k_c * (1 + β * T)

where:
- k_c: base conductivity at reference temperature (W/m·K)
- β: temperature coefficient (1/°C)
- T: local temperature (°C)

This is an assumed linear constitutive law with user-selected parameters. The
sign and magnitude of β are not intended to represent all materials; the law
must remain positive over the computed temperature range.

Governing Equation:
    d/dx[k(T) * dT/dx] = 0

This is a non-linear equation because k depends on T.

Interface Conductivity:
For accuracy at cell interfaces, we use the harmonic mean:
    k_{i+1/2} = 2 * k_i * k_{i+1} / (k_i + k_{i+1})

This properly handles discontinuities in conductivity.

Solution Method:
Since the equation is non-linear, we iterate:
1. Assume initial temperature distribution
2. Calculate k(T) at each node
3. Compute interface conductivities using harmonic mean
4. Set up and solve the linear system
5. Check residual, repeat if not converged
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.numerics import (
    ConvergenceError,
    tdma_solve,
    validate_finite,
    validate_iteration_controls,
    validate_count,
    validate_positive,
)



def conductivity(T: np.ndarray, k_c: float, beta: float) -> np.ndarray:
    """
    Calculates temperature-dependent thermal conductivity.

    k(T) = k_c * (1 + β * T)

    :param T: Temperature array (numpy.ndarray or float)
    :param k_c: Base conductivity coefficient (float, W/m·K)
    :param beta: Temperature coefficient (float, 1/°C)
    :return: Conductivity array (numpy.ndarray)
    """
    validate_positive(k_c=k_c)
    validate_finite(beta=beta)
    T = np.asarray(T, dtype=float)
    if not np.all(np.isfinite(T)):
        raise ValueError("temperature values must be finite")
    values = k_c * (1 + beta * T)
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("thermal conductivity must be finite and greater than zero")
    return values


def solve_heat_transfer_equation(
    n: int,
    L: float,
    k_c: float,
    beta: float,
    h_l: float,
    h_r: float,
    T_l: float,
    T_r: float,
    Rtol: float,
    max_iterations: int = 10000,
) -> tuple:
    """
    Solves the non-linear heat transfer problem with temperature-dependent conductivity.

    Uses iterative approach with TDMA solver until residual converges below tolerance.
    Interface conductivities are computed using harmonic mean for better accuracy.

    :param n: Number of divisions in the rod (int)
    :param L: Length of the rod (float, meters)
    :param k_c: Base conductivity coefficient (float, W/m·K)
    :param beta: Temperature coefficient for conductivity (float, 1/°C)
    :param h_l: Left-side heat transfer coefficient (float, W/m²·K)
    :param h_r: Right-side heat transfer coefficient (float, W/m²·K)
    :param T_l: Left-side fluid temperature (float, °C)
    :param T_r: Right-side fluid temperature (float, °C)
    :param Rtol: Residual tolerance for convergence (float)
    :return: Tuple of (x, T, iterations, k_profile) - positions, temperatures, iteration count, conductivity
    """
    n = validate_count("n", n)
    validate_positive(L=L, k_c=k_c, h_l=h_l, h_r=h_r)
    validate_finite(beta=beta, T_l=T_l, T_r=T_r)
    Rtol, max_iterations = validate_iteration_controls(Rtol, max_iterations)
    dx = L / n
    x = np.linspace(0, L, n + 1)

    # Start from the straight line between ambient temperatures.
    T = np.linspace(T_l, T_r, n + 1)

    def faces_for(temperatures):
        nodal_k = conductivity(temperatures, k_c, beta)
        lower_k = np.minimum(nodal_k[:-1], nodal_k[1:])
        upper_k = np.maximum(nodal_k[:-1], nodal_k[1:])
        return lower_k * (2.0 / (1.0 + lower_k / upper_k))

    def assemble(face_k):
        main = np.zeros(n + 1)
        lower = np.zeros(n + 1)
        upper = np.zeros(n + 1)
        rhs = np.zeros(n + 1)
        conductance = face_k / dx

        main[0] = conductance[0] + h_l
        upper[0] = -conductance[0]
        rhs[0] = h_l * T_l

        for index in range(1, n):
            west, east = conductance[index - 1], conductance[index]
            main[index] = west + east
            lower[index] = -west
            upper[index] = -east

        main[-1] = conductance[-1] + h_r
        lower[-1] = -conductance[-1]
        rhs[-1] = h_r * T_r
        return main, lower, upper, rhs

    def residual_norm(temperatures, face_k, system):
        main, lower, upper, rhs = system
        residual = main * temperatures - rhs
        residual[1:] += lower[1:] * temperatures[:-1]
        residual[:-1] += upper[:-1] * temperatures[1:]
        scale = max(
            1.0,
            float(np.max(np.abs(main * temperatures))),
            float(np.max(np.abs(rhs))),
            float(np.max(np.abs(face_k / dx * np.diff(temperatures)))),
        )
        return float(np.max(np.abs(residual)) / scale)

    residual = float("inf")
    for iterations in range(1, max_iterations + 1):
        system = assemble(faces_for(T))
        updated_T = tdma_solve(*system[:3], system[3])

        # Evaluate the original nonlinear equations using k(T) at the new
        # iterate, including the two convective boundary equations.
        final_faces = faces_for(updated_T)
        final_system = assemble(final_faces)
        residual = residual_norm(updated_T, final_faces, final_system)
        T = updated_T
        if residual <= Rtol:
            return x, T, iterations, conductivity(T, k_c, beta)

    raise ConvergenceError(
        "Temperature-dependent heat-transfer solve did not converge",
        max_iterations,
        residual,
    )


def plot_temperature_profile(
    x: np.ndarray, T: np.ndarray, k_profile: np.ndarray, k_c: float, beta: float
):
    """
    Plots the temperature and conductivity profiles showing the non-linear effects.

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param k_profile: Conductivity values at each point (numpy.ndarray)
    :param k_c: Base conductivity (float)
    :param beta: Temperature coefficient (float)
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Temperature profile
    ax1.plot(x, T, "b-", linewidth=2, label="Temperature")
    ax1.set_xlabel("Position x (m)", fontsize=12)
    ax1.set_ylabel("Temperature (°C)", fontsize=12)
    ax1.set_title("Temperature Profile with Variable Conductivity", fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Conductivity profile
    ax2.plot(x, k_profile, "g-", linewidth=2, label=f"k(T) = {k_c}(1 + {beta}T)")
    ax2.axhline(y=k_c, color="r", linestyle="--", label=f"k₀ = {k_c} W/m·K")
    ax2.set_xlabel("Position x (m)", fontsize=12)
    ax2.set_ylabel("Conductivity (W/m·K)", fontsize=12)
    ax2.set_title("Temperature-Dependent Conductivity Profile", fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

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
    x = np.asarray(x, dtype=float)
    T = np.asarray(T, dtype=float)
    if x.ndim != 1 or T.ndim != 1 or x.size != T.size:
        raise ValueError("x and T must be one-dimensional arrays of equal length")
    validate_count("profile points", x.size, minimum=2)
    validate_count("resolution", resolution, minimum=3)
    validate_positive(radius=radius)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(T)):
        raise ValueError("x and T must contain only finite values")
    points = vtk.vtkPoints()
    temperature_values = vtk.vtkFloatArray()
    temperature_values.SetName("Temperature")
    cells = vtk.vtkCellArray()

    for i, (xi, Ti) in enumerate(zip(x, T)):
        for j in range(resolution):
            theta = 2.0 * np.pi * j / resolution
            y = radius * np.cos(theta)
            z = radius * np.sin(theta)
            points.InsertNextPoint(xi, y, z)
            temperature_values.InsertNextValue(Ti)

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
    Create a 3D text annotation.
    """
    text_actor = vtk.vtkBillboardTextActor3D()
    text_actor.SetInput(text)
    text_actor.SetPosition(position)
    text_actor.GetTextProperty().SetFontSize(18)
    text_actor.GetTextProperty().SetColor(color)
    text_actor.GetTextProperty().SetJustificationToCentered()
    text_actor.GetTextProperty().SetBold(True)

    return text_actor


def vtk_3d_rod_visualization(
    x: np.ndarray,
    T: np.ndarray,
    k_profile: np.ndarray,
    T_l: float,
    T_r: float,
    h_l: float,
    h_r: float,
    k_c: float,
    beta: float,
    title: str = "1D Heat Transfer with Variable Conductivity",
):
    """
    Create an advanced 3D visualization showing temperature-dependent conductivity.

    :param x: X-coordinate vector (numpy.ndarray)
    :param T: Temperature field (numpy.ndarray)
    :param k_profile: Conductivity profile (numpy.ndarray)
    :param T_l: Left fluid temperature (float)
    :param T_r: Right fluid temperature (float)
    :param h_l: Left heat transfer coefficient (float)
    :param h_r: Right heat transfer coefficient (float)
    :param k_c: Base conductivity (float)
    :param beta: Temperature coefficient (float)
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
    rod_mapper.UseLookupTableScalarRangeOn()
    rod_mapper.InterpolateScalarsBeforeMappingOn()

    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(rod_mapper)
    # Preserve the exact lookup-table colors shown by the scalar bar.
    rod_actor.GetProperty().LightingOff()
    rod_actor.GetProperty().SetInterpolationToPhong()
    rod_actor.GetProperty().SetSpecular(0.3)
    rod_actor.GetProperty().SetSpecularPower(20)

    # Boundary annotations
    annotation_offset = radius * 3
    left_text = f"Convective BC\nT∞ = {T_l}°C\nh = {h_l} W/m²K"
    right_text = f"Convective BC\nT∞ = {T_r}°C\nh = {h_r} W/m²K"
    material_text = f"Variable k(T)\nk = {k_c}(1 + {beta}T)"

    left_annotation = create_boundary_annotation(
        left_text, (x[0] - 0.12 * L, 0, annotation_offset), color=(0.9, 0.95, 1.0)
    )
    right_annotation = create_boundary_annotation(
        right_text,
        (x[-1] + 0.12 * L, 0, -radius * 3.5),
        color=(1.0, 0.9, 0.9),
    )
    material_annotation = create_boundary_annotation(
        material_text, (L / 2, 0, radius * 2.5), color=(0.8, 1.0, 0.8)
    )

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature\n(°C)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.06)
    scalar_bar.SetHeight(0.4)
    scalar_bar.SetPosition(0.91, 0.3)
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.GetTitleTextProperty().SetFontSize(16)
    scalar_bar.GetLabelTextProperty().SetFontSize(14)

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
    renderer.AddActor(left_annotation)
    renderer.AddActor(right_annotation)
    renderer.AddActor(material_annotation)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)

    # Gradient background
    renderer.SetBackground(0.1, 0.15, 0.1)
    renderer.SetBackground2(0.02, 0.05, 0.02)
    renderer.GradientBackgroundOn()

    # Lighting
    light = vtk.vtkLight()
    light.SetLightTypeToHeadlight()
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
    Main function demonstrating 1D heat transfer with temperature-dependent conductivity.

    Part 3 of the Heat Transfer Series: Non-linear material properties
    """
    # Physical parameters
    L = 1.0  # Rod length (m)
    n = 50  # Number of divisions
    k_c = 50.0  # Base conductivity (W/m·K)
    beta = 0.001  # Temperature coefficient (1/°C) - positive means k increases with T
    h_l = 100.0  # Left heat transfer coefficient (W/m²·K)
    h_r = 50.0  # Right heat transfer coefficient (W/m²·K)
    T_l = 100.0  # Left fluid temperature (°C)
    T_r = 300.0  # Right fluid temperature (°C)
    Rtol = 1e-8  # Convergence tolerance

    # Solve the non-linear heat transfer problem
    x, T, iterations, k_profile = solve_heat_transfer_equation(
        n, L, k_c, beta, h_l, h_r, T_l, T_r, Rtol
    )

    # Plot results
    plot_temperature_profile(x, T, k_profile, k_c, beta)

    # 3D VTK visualization
    vtk_3d_rod_visualization(x, T, k_profile, T_l, T_r, h_l, h_r, k_c, beta)

    # Print results
    print("=" * 70)
    print("Part 3: 1D Heat Transfer with Temperature-Dependent Conductivity")
    print("=" * 70)
    print(f"Rod length: {L} m, Grid points: {n + 1}")
    print(f"Base conductivity: k₀ = {k_c} W/m·K")
    print(f"Temperature coefficient: β = {beta} 1/°C")
    print("Conductivity law: k(T) = k₀(1 + βT)")
    print("-" * 70)
    print(f"Left BC:  h = {h_l} W/m²·K, T∞ = {T_l}°C")
    print(f"Right BC: h = {h_r} W/m²·K, T∞ = {T_r}°C")
    print("-" * 70)
    print(f"Converged in {iterations} iterations")
    print(f"Temperature at x=0: {T[0]:.2f}°C")
    print(f"Temperature at x=L: {T[-1]:.2f}°C")
    print(f"Conductivity range: {k_profile.min():.3f} - {k_profile.max():.3f} W/m·K")
    print(
        f"Conductivity variation: {100 * (k_profile.max() / k_profile.min() - 1):.1f}%"
    )
    print("=" * 70)


if __name__ == "__main__":
    main()
