"""Steady thermal development in a laminar parallel-plate channel.

The prescribed parabolic velocity has unit mean speed and zero speed at y=0,h.
The energy equation is rho*c*u(y)*T_x = k*(T_xx + T_yy). Both walls receive
constant heat flux hf (W/m², positive into the fluid), the inlet is Ti, and the
outlet has zero conductive boundary flux. This solves a temperature field in
an assumed flow, not the velocity or pressure of a fluid. A nodal finite-volume
scheme with upwind convection conserves fluxes, including boundary volumes.
"""

import matplotlib.pyplot as plt
import numpy as np
import vtk

import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.heat_transfer import (
    channel_velocity,
    channel_weights,
    solve_channel,
    validate_field,
)
from src.common.field_visualization import field_points, quad_field
from src.common.numerics import validate_count, validate_finite, validate_positive


def initialize_grid(length: float, h: float, dx: float, dy: float) -> tuple:
    """Build an endpoint-inclusive nodal grid with the requested spacings.

    Both dimensions must contain a whole number of intervals. Roundoff in
    ratios such as 0.3/0.1 is tolerated; a genuinely nondivisible grid is rejected
    so coordinates, temperature dimensions, and the operator cannot disagree.
    """
    validate_positive(length=length, h=h, dx=dx, dy=dy)
    counts = []
    for extent, spacing in ((length, dx), (h, dy)):
        ratio = extent / spacing
        count = round(ratio)
        validate_count("grid intervals", count, minimum=2)
        if not np.isclose(ratio, count, rtol=1e-12, atol=1e-12):
            raise ValueError(
                "Domain dimensions must be whole multiples of grid spacing."
            )
        counts.append(count)
    nx, ny = counts
    return nx, ny, np.linspace(0, length, nx + 1), np.linspace(0, h, ny + 1)


def initialize_temperature_array(ny: int, nx: int, Ti: float) -> np.ndarray:
    """
    Initialize the temperature array.

    Args:
        ny (int): Number of grid points in the y-direction.
        nx (int): Number of grid points in the x-direction.
        Ti (float): Initial temperature.

    Returns:
        np.ndarray: Initialized temperature array.
    """
    validate_count("ny", ny, minimum=2)
    validate_count("nx", nx, minimum=2)
    validate_finite(Ti=Ti)
    return np.full((ny + 1, nx + 1), float(Ti))


def velocity(yp: float, h: float):
    """Prescribed parallel-plate velocity for y in [0,h], mean speed 1 m/s."""
    return channel_velocity(yp, h)


def update_temperature(
    Aw: float,
    Ae: float,
    An: float,
    As: float,
    b: float,
    T: np.ndarray,
    i: int,
    j: int,
    nx: int,
    ny: int,
) -> float:
    """
    Update temperature at a grid point.

    Args:
        Aw (float): West diffusion coefficient.
        Ae (float): East diffusion coefficient.
        An (float): North diffusion coefficient.
        As (float): South diffusion coefficient.
        b (float): Source term.
        T (np.ndarray): Temperature array.
        i (int): Row index of the grid point.
        j (int): Column index of the grid point.
        nx (int): Number of grid points in the x-direction.
        ny (int): Number of grid points in the y-direction.

    Returns:
        float: Absolute value of the temperature update.
    """
    Ap = Aw + Ae + An + As
    tmp = b - Ap * T[i, j]
    if j > 0:
        tmp += Aw * T[i, j - 1]
    if j < nx:
        tmp += Ae * T[i, j + 1]
    if i < ny:
        tmp += An * T[i + 1, j]
    if i > 0:
        tmp += As * T[i - 1, j]
    T[i, j] += tmp / Ap
    return abs(tmp)


def apply_boundary_conditions(
    T: np.ndarray,
    Ti: float,
    nx: int,
    ny: int,
    k: float,
    dx: float,
    dy: float,
    hf: float,
):
    """Relax the two wall control volumes; positive hf heats both walls.

    This standalone boundary demonstration uses conductive face balances.
    run_simulation also includes convection in these volumes and balances the
    outlet advective flux through the same finite-volume kernel as the interior.
    """
    validate_positive(k=k, dx=dx, dy=dy)
    if T.shape != (ny + 1, nx + 1):
        raise ValueError("Temperature dimensions must match nx and ny.")
    T[:, 0], T[:, -1] = Ti, T[:, -2]
    for j in range(1, nx):
        west = east = k * dy / (2 * dx)
        inward = k * dx / dy
        source = hf * dx
        update_temperature(west, east, inward, 0, source, T, 0, j, nx, ny)
        update_temperature(west, east, 0, inward, source, T, ny, j, nx, ny)


def run_simulation(
    T: np.ndarray,
    nx: int,
    ny: int,
    Rtol: float,
    dx: float,
    dy: float,
    y: np.ndarray,
    Ti: float,
    hf: float,
    rho: float,
    c: float,
    k: float,
    h: float,
    max_iterations: int = 50000,
) -> np.ndarray:
    """Solve steady energy transport in prescribed laminar channel flow.

    The channel is y in [0,h], with unit mean speed. hf is imposed wall heat
    flux in W/m², positive into the fluid at both walls. Ti fixes the inlet.
    The outlet carries advective energy with zero conductive boundary flux.
    The kernel accounts for half wall/outlet control volumes and checks their
    equations as part of convergence. Rtol is a temperature residual (K).
    It returns a new T[y,x]; nonconvergence raises ConvergenceError.
    """
    validate_count("nx", nx, minimum=2)
    validate_count("ny", ny, minimum=2)
    if np.shape(T) != (ny + 1, nx + 1):
        raise ValueError("Temperature dimensions must match nx and ny.")
    field, _, _ = solve_channel(
        T, dx, dy, y, Ti, hf, rho, c, k, h, Rtol, max_iterations
    )
    return field


def calculate_bulk_temperature(
    T: np.ndarray,
    y: np.ndarray,
    dy: float,
    h: float,
    nx: int,
) -> np.ndarray:
    """Velocity-weighted bulk temperature using the solver's face flux weights."""
    validate_positive(dy=dy, h=h)
    validate_count("nx", nx, minimum=2)
    _, y, T = validate_field(np.arange(nx + 1), y, T)
    if not np.isclose(dy * (len(y) - 1), h) or not np.allclose(
        y, np.arange(len(y)) * dy
    ):
        raise ValueError("y, dy, and h must describe the same channel grid.")
    weights, _ = channel_weights(y, h)
    return np.sum(T * weights[:, None], axis=0) / np.sum(weights)


def calculate_heat_flux(T: np.ndarray, k: float, dy: float) -> np.ndarray:
    """Estimate inward bottom-wall heat flux with a second-order derivative.

    This derivative diagnostic approaches the prescribed boundary-face flux
    under refinement; the solver imposes that flux by a control-volume balance.
    The inlet corner has a prescribed temperature rather than a wall flux.
    """
    validate_positive(k=k, dy=dy)
    if T.ndim != 2 or T.shape[0] < 3 or not np.all(np.isfinite(T)):
        raise ValueError("A finite field with at least three y nodes is required.")
    return -k * (-3 * T[0, :] + 4 * T[1, :] - T[2, :]) / (2 * dy)


def calculate_ht_and_nusselt(
    Tb: np.ndarray,
    q_w: np.ndarray,
    k: float,
    h: float,
    T: np.ndarray,
) -> tuple:
    """Compute bottom-wall h_t and Nu with hydraulic diameter 2h.

    A coefficient is undefined where wall and bulk temperatures coincide (for
    example at the inlet). Report NaN there instead of fabricating zero or inf.
    """
    validate_positive(k=k, h=h)
    T, Tb, q_w = (np.asarray(value, dtype=float) for value in (T, Tb, q_w))
    if T.ndim != 2 or Tb.shape != (T.shape[1],):
        raise ValueError("Bulk temperature must match the x grid.")
    if not all(np.all(np.isfinite(value)) for value in (T, Tb, q_w)):
        raise ValueError("Temperatures and heat flux must be finite.")
    difference = np.asarray(T[0, :] - Tb, dtype=float)
    if np.shape(q_w) != difference.shape:
        raise ValueError("Heat flux and bulk temperature must match the x grid.")
    ht = np.full_like(difference, np.nan)
    np.divide(q_w, difference, out=ht, where=np.abs(difference) > 1e-12)
    return ht, ht * (2 * h) / k


def plot_results(
    x: np.ndarray,
    y: np.ndarray,
    T: np.ndarray,
    Ti: float,
    Tb: np.ndarray,
    ht: np.ndarray,
    Nu: np.ndarray,
):
    """
    Plot simulation results.

    Args:
        x (np.ndarray): X-coordinates.
        y (np.ndarray): Y-coordinates.
        T (np.ndarray): Temperature array.
        Ti (float): Inlet temperature.
        Tb (np.ndarray): Bulk temperature values.
        ht (np.ndarray): Heat transfer coefficients.
        Nu (np.ndarray): Nusselt numbers.
    """
    plt.figure(figsize=(12, 8))
    plot_titles = [
        "Temperature Distribution",
        "Bulk Temperature",
        "Heat Transfer Coefficient",
        "Nusselt Number",
    ]

    for i in range(4):
        plt.subplot(2, 2, i + 1)
        if i == 0:
            plt.contourf(x, y, T, levels=np.linspace(Ti, np.max(T), 50))
            plt.colorbar()
        elif i == 1:
            plt.plot(x, Tb)
            plt.ylabel("Bulk Temperature (°C)")
        elif i == 2:
            plt.plot(x, ht)
            plt.ylabel("Heat Transfer Coefficient")
        elif i == 3:
            plt.plot(x, Nu)
            plt.ylabel("Nusselt Number")

        plt.xlabel("x (m)")
        plt.title(plot_titles[i])

    plt.tight_layout()
    plt.show()


def create_points(x: np.ndarray, y: np.ndarray, T: np.ndarray) -> tuple:
    """Attach T[y,x] to its physical coordinates without displacement."""
    return field_points(x, y, T)


def create_poly_data(nx: int, ny: int, points, temperature_values):
    """Connect the scalar grid with consistently oriented quads."""
    return quad_field(nx, ny, points, temperature_values)


def create_color_map(
    T: np.ndarray, color_range=None, hue_range=(0.667, 0)
) -> vtk.vtkLookupTable:
    """
    Create a color map for VTK visualization.

    Args:
        T (np.ndarray): Temperature array.
        color_range (tuple): Range for color mapping (min, max).
        hue_range (tuple): Range for hue variation.

    Returns:
        vtk.vtkLookupTable: VTK lookup table for color mapping.
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(T), np.max(T)
    color_map.SetRange(color_range if color_range else (min_val, max_val))
    color_map.SetHueRange(*hue_range)
    color_map.Build()
    return color_map


def vtk_surface_plot(x: np.ndarray, y: np.ndarray, T: np.ndarray):
    """
    Visualize the temperature distribution using VTK.

    Args:
        x (np.ndarray): X-coordinates.
        y (np.ndarray): Y-coordinates.
        T (np.ndarray): Temperature array.
        color_range (tuple): Range for color mapping (min, max).
        hue_range (tuple): Range for hue variation.
    """
    points, temperature_values = create_points(x, y, T)
    polyData = create_poly_data(len(x), len(y), points, temperature_values)
    color_map = create_color_map(T)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.InterpolateScalarsBeforeMappingOn()
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetScalarRange(np.min(T), np.max(T))

    # Surface actor setup
    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)
    surface_actor.GetProperty().LightingOff()
    surface_actor.GetProperty().SetLineWidth(5)

    # Scalar bar (legend)
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.UnconstrainedFontSizeOn()
    scalar_bar.GetTitleTextProperty().SetFontSize(16)
    scalar_bar.GetLabelTextProperty().SetFontSize(14)
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.09)
    scalar_bar.SetHeight(0.45)
    scalar_bar.SetPosition(0.88, 0.28)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(scalar_bar)
    annotation = vtk.vtkTextActor()
    annotation.SetInput(
        "Prescribed laminar flow: mean u = 1 m/s\n"
        "Physical x/y coordinates (m); temperature shown by color"
    )
    annotation.SetPosition(12, 12)
    annotation.GetTextProperty().SetFontSize(14)
    renderer.AddActor2D(annotation)

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1000, 600)
    renderWindow.SetWindowName("Channel thermal development")

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
    # Simulation parameters
    length = 20.0
    h = 1.0
    dx = 0.1
    dy = 0.1
    rho = 1.0
    c = 100.0
    k = 1.0
    hf = 1000.0
    Ti = 50.0
    Rtol = 1e-3

    # Initialize grid and temperature arrays
    nx, ny, x, y = initialize_grid(length, h, dx, dy)
    T = initialize_temperature_array(ny, nx, Ti)

    # Run the simulation
    T = run_simulation(
        T=T,
        nx=nx,
        ny=ny,
        Rtol=Rtol,
        dx=dx,
        dy=dy,
        y=y,
        Ti=Ti,
        hf=hf,
        rho=rho,
        c=c,
        k=k,
        h=h,
    )

    # Use the prescribed boundary-face heat flux for h_t and Nu. A numerical
    # wall derivative (calculate_heat_flux) is a separate refinement diagnostic.
    Tb = calculate_bulk_temperature(T, y, dy, h, nx)
    q_w = np.full(nx + 1, hf)

    # Calculate heat transfer and Nusselt number
    ht, Nu = calculate_ht_and_nusselt(Tb, q_w, k, h, T)

    # Plot the results
    plot_results(x, y, T, Ti, Tb, ht, Nu)
    vtk_surface_plot(x, y, T)


if __name__ == "__main__":
    main()
