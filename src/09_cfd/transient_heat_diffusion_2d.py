"""
Transient 2D Heat Diffusion: a beginner's first time-marching solver

A thin steel plate starts with two hot spots while its edges are held at
room temperature. We march the heat equation forward in time with the explicit
FTCS scheme (Forward in Time, Centred in Space) and watch the spots spread out,
merge, and drain through the cold edges.

Physics:
    dT/dt = alpha * (d²T/dx² + d²T/dy²)

where alpha = k / (rho * c_p) is the thermal diffusivity (m²/s).

Discretization (FTCS), for every interior node (i, j):
    T_new = T + alpha*dt * ((T_E - 2T + T_W)/dx² + (T_N - 2T + T_S)/dy²)

The scheme is explicit: each new value only needs old neighbours, so one time
step is a single NumPy expression. The price is a stability limit. With the
Fourier numbers Fo_x = alpha*dt/dx² and Fo_y = alpha*dt/dy², FTCS is stable
only when

    Fo_x + Fo_y <= 1/2        (dx = dy:  alpha*dt/dx² <= 1/4)

Take a larger step and errors grow without bound. Try it: set ``safety`` above
1 in :func:`stable_time_step` and watch a checkerboard explode.

Verification:
A single sine mode sin(pi x/L) sin(pi y/L) with T = 0 on the edges is an exact
eigenfunction of the discrete operator. Each step multiplies it by
1 - 4*Fo*sin²(pi*dx/(2L)) (per direction), and the continuous solution decays
like exp(-2*pi²*alpha*t/L²). The tests check both.

Run:
    python -m src transient_heat_diffusion_2d
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.field_visualization import field_points, quad_field

# Carbon steel: k = 45 W/(m·K), rho = 7850 kg/m³, c_p = 490 J/(kg·K)
STEEL_DIFFUSIVITY = 45.0 / (7850.0 * 490.0)  # ≈ 1.17e-5 m²/s
DEFAULT_SPOTS = (
    # (x, y, radius, temperature rise) in metres and °C
    (0.032, 0.062, 0.012, 180.0),
    (0.068, 0.036, 0.010, 140.0),
)


def stable_time_step(alpha, dx, dy, safety=0.9):
    """Largest FTCS step with Fo_x + Fo_y <= 1/2, scaled by ``safety``."""
    if alpha <= 0 or dx <= 0 or dy <= 0:
        raise ValueError("alpha, dx, and dy must be positive")
    return safety * 0.5 / (alpha * (1 / dx**2 + 1 / dy**2))


def hot_spots(x, y, background, spots=DEFAULT_SPOTS):
    """Initial field: Gaussian hot spots on a uniform background (°C)."""
    X, Y = np.meshgrid(x, y)
    T = np.full(X.shape, float(background))
    for cx, cy, radius, rise in spots:
        T += rise * np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / radius**2)
    return T


def ftcs_step(T, alpha, dt, dx, dy):
    """Advance one explicit step; boundary rows/columns keep their values."""
    new = T.copy()
    new[1:-1, 1:-1] = T[1:-1, 1:-1] + alpha * dt * (
        (T[1:-1, 2:] - 2 * T[1:-1, 1:-1] + T[1:-1, :-2]) / dx**2
        + (T[2:, 1:-1] - 2 * T[1:-1, 1:-1] + T[:-2, 1:-1]) / dy**2
    )
    return new


def simulate(
    n=80,
    length=0.1,
    alpha=STEEL_DIFFUSIVITY,
    edge_temperature=20.0,
    duration=60.0,
    snapshots=121,
    safety=0.9,
    spots=DEFAULT_SPOTS,
    progress_callback=None,
):
    """March from the hot-spot field and record ``snapshots`` frames.

    Snapshot times are spaced more densely early on, where the field changes
    fastest. Returns ``(x, y, times, steps, frames, dt)`` with frames shaped
    ``(snapshots, n+1, n+1)`` and indexed ``[time, y, x]``.
    """
    x = np.linspace(0, length, n + 1)
    y = np.linspace(0, length, n + 1)
    dx = dy = length / n
    dt = stable_time_step(alpha, dx, dy, safety)
    total = int(np.ceil(duration / dt))
    targets = np.unique(np.rint(total * np.linspace(0, 1, snapshots) ** 1.6).astype(int))
    T = hot_spots(x, y, edge_temperature, spots)
    T[0, :] = T[-1, :] = T[:, 0] = T[:, -1] = edge_temperature
    frames = []
    for step in range(total + 1):
        if step == targets[len(frames)]:
            frames.append(T.copy())
            if progress_callback is not None:
                progress_callback(step, total)
            if len(frames) == len(targets):
                break
        T = ftcs_step(T, alpha, dt, dx, dy)
    return x, y, targets * dt, targets, np.array(frames), dt


def plot_snapshots(x, y, times, frames):
    """Four moments of the cooling plate on a shared color scale."""
    picks = np.linspace(0, len(frames) - 1, 4).astype(int)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4), constrained_layout=True)
    low, high = frames.min(), frames.max()
    for axis, index in zip(axes, picks):
        image = axis.imshow(
            frames[index],
            origin="lower",
            extent=(x[0] * 100, x[-1] * 100, y[0] * 100, y[-1] * 100),
            cmap="inferno",
            vmin=low,
            vmax=high,
        )
        axis.set_title(f"t = {times[index]:.1f} s")
        axis.set_xlabel("x (cm)")
    axes[0].set_ylabel("y (cm)")
    fig.colorbar(image, ax=axes, label="Temperature (°C)", shrink=0.9)
    fig.suptitle("Explicit (FTCS) heat diffusion in a steel plate")
    plt.show()


def animate_vtk(x, y, times, frames):
    """Play the recorded frames as a warped, colored surface."""
    low, high = float(frames.min()), float(frames.max())
    span = float(x[-1] - x[0])
    height = 0.35 * span
    points, scalars = field_points(x, y, frames[0], height=0.0)
    mesh = quad_field(len(x), len(y), points, scalars)

    table = vtk.vtkLookupTable()
    table.SetRange(low, high)
    table.SetHueRange(0.66, 0.0)
    table.Build()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetLookupTable(table)
    mapper.UseLookupTableScalarRangeOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    bar = vtk.vtkScalarBarActor()
    bar.SetLookupTable(table)
    bar.SetTitle("T (°C)")
    clock = vtk.vtkTextActor()
    clock.SetPosition(20, 20)
    clock.GetTextProperty().SetFontSize(22)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor2D(bar)
    renderer.AddActor2D(clock)
    renderer.SetBackground(0.06, 0.08, 0.12)
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(1200, 850)
    window.SetWindowName("Transient 2D heat diffusion")
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

    state = {"frame": 0}

    def show(index):
        field = frames[index]
        values = mesh.GetPointData().GetScalars()
        for point, value in enumerate(field.ravel()):
            px, py, _ = mesh.GetPoint(point)
            mesh.GetPoints().SetPoint(point, px, py, (value - low) / (high - low) * height)
            values.SetValue(point, value)
        mesh.GetPoints().Modified()
        values.Modified()
        clock.SetInput(f"t = {times[index]:5.1f} s   (height and color: temperature)")

    def tick(*_):
        state["frame"] = (state["frame"] + 1) % len(frames)
        show(state["frame"])
        window.Render()

    show(0)
    camera = renderer.GetActiveCamera()
    camera.SetPosition(x[-1] * 1.6, -y[-1] * 1.4, span * 1.5)
    camera.SetFocalPoint(span / 2, span / 2, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    interactor.Initialize()
    interactor.AddObserver("TimerEvent", tick)
    interactor.CreateRepeatingTimer(60)
    window.Render()
    interactor.Start()


def main():
    x, y, times, steps, frames, dt = simulate()
    fo = STEEL_DIFFUSIVITY * dt / (x[1] - x[0]) ** 2
    print("Transient 2D heat diffusion (explicit FTCS)")
    print(f"  plate 10 cm × 10 cm, {len(x)} × {len(y)} nodes, alpha = {STEEL_DIFFUSIVITY:.3e} m²/s")
    print(f"  dt = {dt * 1000:.2f} ms, Fo = {fo:.3f} per direction (stable when <= 0.25)")
    print(f"  {steps[-1]:,} steps to t = {times[-1]:.0f} s")
    print(f"  peak temperature {frames[0].max():.1f} °C -> {frames[-1].max():.1f} °C")
    plot_snapshots(x, y, times, frames)
    animate_vtk(x, y, times, frames)


if __name__ == "__main__":
    main()
