"""
Rotating Blob Advection: see numerical diffusion with your own eyes

A blob of dye sits in a fluid that spins like a turntable (solid-body
rotation). Physically the blob should ride around the centre unchanged: pure
advection has no diffusion term. We solve it with the simplest stable scheme,
first-order upwind, and watch the blob smear and its peak drop anyway. That
smearing is numerical diffusion, a property of the scheme rather than of the
physics.

Physics (linear advection, no diffusion):
    dphi/dt + u * dphi/dx + v * dphi/dy = 0
    u = -omega * (y - yc),   v = omega * (x - xc)

The velocity is divergence-free, so the conservative and advective forms agree
and the total amount of dye is conserved while the blob stays away from the
edges.

Discretization (explicit, first-order upwind):
Each face takes phi from its upwind side:
    F_east = u * (phi_P if u > 0 else phi_E)
    phi_new = phi - dt/dx * (F_east - F_west) - dt/dy * (G_north - G_south)

Stability (CFL condition):
    CFL = dt * (|u|/dx + |v|/dy) <= 1

Under this limit, upwind is monotone: phi never rises above its initial maximum
or below its minimum. Beginners can compare the numerical peak with the exact
peak (always 1) and see the error shrink with grid refinement.

Run:
    python -m src rotating_blob_advection_2d
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.field_visualization import field_points, quad_field


def rotation_velocity(x, y, omega=2 * np.pi, center=(0.5, 0.5)):
    """Solid-body rotation, counterclockwise; returns (u, v) on the grid."""
    X, Y = np.meshgrid(x, y)
    return -omega * (Y - center[1]), omega * (X - center[0])


def gaussian_blob(x, y, center=(0.5, 0.7), radius=0.07):
    X, Y = np.meshgrid(x, y)
    return np.exp(-((X - center[0]) ** 2 + (Y - center[1]) ** 2) / radius**2)


def exact_blob(x, y, time, omega=2 * np.pi, center=(0.5, 0.5), start=(0.5, 0.7), radius=0.07):
    """The initial blob rotated by omega*t about the centre."""
    angle = omega * time
    dx, dy = start[0] - center[0], start[1] - center[1]
    moved = (
        center[0] + dx * np.cos(angle) - dy * np.sin(angle),
        center[1] + dx * np.sin(angle) + dy * np.cos(angle),
    )
    return gaussian_blob(x, y, moved, radius)


def cfl_time_step(u, v, dx, dy, cfl=0.8):
    """Explicit step for a target CFL number (stable for cfl <= 1)."""
    rate = np.max(np.abs(u) / dx + np.abs(v) / dy)
    if rate <= 0:
        raise ValueError("velocity field must be nonzero")
    return cfl / rate


def upwind_step(phi, u, v, dt, dx, dy):
    """One conservative first-order upwind step with zero inflow at the edges."""
    padded = np.pad(phi, 1)  # phi = 0 outside the domain
    # Face velocities: u varies only with y and v only with x, so node values
    # are exact face values.
    east = np.where(u > 0, padded[1:-1, 1:-1], padded[1:-1, 2:]) * u
    west = np.where(u > 0, padded[1:-1, :-2], padded[1:-1, 1:-1]) * u
    north = np.where(v > 0, padded[1:-1, 1:-1], padded[2:, 1:-1]) * v
    south = np.where(v > 0, padded[:-2, 1:-1], padded[1:-1, 1:-1]) * v
    return phi - dt / dx * (east - west) - dt / dy * (north - south)


def simulate(n=100, revolutions=1.0, cfl=0.8, snapshots=121, omega=2 * np.pi, progress_callback=None):
    """Rotate the blob and record evenly spaced frames.

    Returns ``(x, y, times, steps, frames, dt, (u, v))``.
    """
    x = np.linspace(0, 1, n + 1)
    y = np.linspace(0, 1, n + 1)
    dx = dy = 1 / n
    u, v = rotation_velocity(x, y, omega)
    period = 2 * np.pi / omega
    duration = revolutions * period
    dt = cfl_time_step(u, v, dx, dy, cfl)
    total = int(np.ceil(duration / dt))
    dt = duration / total  # land exactly on the final time
    targets = np.unique(np.rint(np.linspace(0, total, snapshots)).astype(int))
    phi = gaussian_blob(x, y)
    frames = []
    for step in range(total + 1):
        if step == targets[len(frames)]:
            frames.append(phi.copy())
            if progress_callback is not None:
                progress_callback(step, total)
            if len(frames) == len(targets):
                break
        phi = upwind_step(phi, u, v, dt, dx, dy)
    return x, y, targets * dt, targets, np.array(frames), dt, (u, v)


def plot_comparison(x, y, times, frames):
    """Initial blob, the upwind result after one turn, and peak decay."""
    exact = exact_blob(x, y, times[-1])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for axis, field, title in (
        (axes[0], frames[0], "Start"),
        (axes[1], frames[-1], "Upwind after one turn"),
    ):
        image = axis.imshow(field, origin="lower", extent=(0, 1, 0, 1), cmap="magma", vmin=0, vmax=1)
        axis.contour(x, y, exact, levels=[0.5], colors="cyan", linewidths=1)
        axis.set_title(title + "  (cyan: exact 0.5 contour)")
    fig.colorbar(image, ax=axes[:2], label="phi")
    axes[2].plot(times, frames.max(axis=(1, 2)), label="Upwind peak")
    axes[2].axhline(1.0, color="k", linestyle="--", label="Exact peak")
    axes[2].set_xlabel("t (s)")
    axes[2].set_ylabel("peak phi")
    axes[2].legend()
    axes[2].set_title("Numerical diffusion lowers the peak")
    plt.show()


def animate_vtk(x, y, times, frames):
    """Play the blob's journey as a warped, colored surface."""
    points, scalars = field_points(x, y, frames[0], height=0.0)
    mesh = quad_field(len(x), len(y), points, scalars)
    table = vtk.vtkLookupTable()
    table.SetRange(0, 1)
    table.SetHueRange(0.7, 0.0)
    table.Build()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetLookupTable(table)
    mapper.UseLookupTableScalarRangeOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    clock = vtk.vtkTextActor()
    clock.SetPosition(20, 20)
    clock.GetTextProperty().SetFontSize(22)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor2D(clock)
    renderer.SetBackground(0.06, 0.08, 0.12)
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(1200, 850)
    window.SetWindowName("Rotating blob advection (first-order upwind)")
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    state = {"frame": 0}

    def show(index):
        field = frames[index]
        values = mesh.GetPointData().GetScalars()
        for point, value in enumerate(field.ravel()):
            px, py, _ = mesh.GetPoint(point)
            mesh.GetPoints().SetPoint(point, px, py, 0.3 * value)
            values.SetValue(point, value)
        mesh.GetPoints().Modified()
        values.Modified()
        clock.SetInput(f"t = {times[index]:.2f} s   peak = {field.max():.3f} (exact: 1)")

    def tick(*_):
        state["frame"] = (state["frame"] + 1) % len(frames)
        show(state["frame"])
        window.Render()

    show(0)
    camera = renderer.GetActiveCamera()
    camera.SetPosition(1.4, -1.0, 1.4)
    camera.SetFocalPoint(0.5, 0.5, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    interactor.Initialize()
    interactor.AddObserver("TimerEvent", tick)
    interactor.CreateRepeatingTimer(60)
    window.Render()
    interactor.Start()


def main():
    x, y, times, steps, frames, dt, _ = simulate()
    print("Rotating blob advection (first-order upwind)")
    print(f"  grid {len(x)} × {len(y)}, dt = {dt * 1000:.3f} ms, {steps[-1]:,} steps for one turn")
    print(f"  peak: exact 1.000, upwind {frames[-1].max():.3f}")
    print(f"  total dye: start {frames[0].sum():.4f}, end {frames[-1].sum():.4f}")
    plot_comparison(x, y, times, frames)
    animate_vtk(x, y, times, frames)


if __name__ == "__main__":
    main()
