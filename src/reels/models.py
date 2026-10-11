"""Headless numerical adapters: reel fields come from the actual example solvers."""

from contextlib import redirect_stdout
from dataclasses import dataclass, field
import hashlib
import importlib
import io
import json
from pathlib import Path
import tempfile

import numpy as np

from src.common.heat_transfer import (
    solve_channel,
    solve_quarter_conduction,
    solve_transport,
)
from src.common.numerics import tdma_solve


@dataclass
class Model:
    key: str
    kind: str
    x: np.ndarray
    y: np.ndarray
    fields: np.ndarray  # (snapshots, y, x); rods have one y row
    unit: str
    field_name: str
    params: dict
    note: str
    metrics: list[tuple[str, str]]
    vectors: np.ndarray | None = None  # (snapshots, 2, y, x)
    solid: np.ndarray | None = None
    steps: np.ndarray | None = None
    chart: dict = field(default_factory=dict)
    cache_hit: bool = False


def example_module(key):
    return importlib.import_module(f"src.09_cfd.{key}")


def _solve(key, quality, progress):
    fine = quality == "studio"
    if key == "obstacle_flow_simulation":
        module = example_module(key)
        # Re is based on actual obstacle height. Coarse grids need more viscosity.
        params = dict(
            nx=240 if fine else 180,
            ny=100 if fine else 76,
            max_iter=2400 if fine else 1600,
            reynolds=50.0,
            inlet_speed=0.05,
            obstacle_length=30 if fine else 24,
            obstacle_height=18 if fine else 14,
            snapshot_interval=20,
        )
        snapshots, solid, _, _ = module.simulate_flow(
            **params, progress_callback=lambda current, total: progress(current / total)
        )
        vectors = np.asarray(snapshots)
        fields = np.hypot(vectors[:, 0], vectors[:, 1])
        fields[:, solid] = np.nan
        steps = np.arange(
            params["snapshot_interval"],
            params["max_iter"] + 1,
            params["snapshot_interval"],
        )
        return Model(
            key,
            "flow",
            np.arange(params["nx"], dtype=float),
            np.arange(params["ny"], dtype=float),
            fields,
            "lattice units",
            "SPEED",
            params,
            "Simulated snapshots · interpolated playback · periodic y",
            [("Re", "50"), ("Recorded steps", str(params["max_iter"]))],
            vectors=vectors,
            solid=solid,
            steps=steps,
        )
    if key == "fluid_flow_simulator":
        module = example_module(key)
        params = dict(
            nx=41 if fine else 31,
            ny=41 if fine else 31,
            reynolds=100.0,
            u_lid=1.0,
            max_iter=100000,
            print_interval=100000,
            tolerance=1e-6,
        )
        with redirect_stdout(io.StringIO()):
            x, y, u, v, p = module.run_simulation(**params)
        progress(1)
        div = module.compute_divergence(u, v, x[1] - x[0], y[1] - y[0])
        return Model(
            key,
            "flow",
            x,
            y,
            np.hypot(u, v)[None],
            "m/s",
            "SPEED",
            params,
            "Steady velocity field · animated tracers are illustrative",
            [("Re", "100"), ("Max divergence", f"{div:.1e}")],
            vectors=np.stack((u, v))[None],
            chart={
                "x": y,
                "series": [u[:, len(x) // 2]],
                "labels": ["Centerline u"],
                "xlabel": "y (m)",
                "ylabel": "u (m/s)",
            },
        )
    if key == "transient_heat_diffusion_2d":
        module = example_module(key)
        # Steel sheds most of its heat within ~25 s; later frames would be uniformly cold.
        params = dict(n=96 if fine else 64, duration=24.0, snapshots=181 if fine else 121)
        x, y, times, steps, frames, dt = module.simulate(
            **params, progress_callback=lambda current, total: progress(current / total)
        )
        fo = module.STEEL_DIFFUSIVITY * dt / (x[1] - x[0]) ** 2
        return Model(
            key,
            "surface",
            x * 100,
            y * 100,
            frames,
            "°C",
            "TEMPERATURE",
            {**params, "dt": float(dt), "fourier": float(fo), "length_unit": "cm"},
            "Explicit time steps · x/y in cm · height: temperature",
            [("Fourier no.", f"{fo:.3f}"), ("Time steps", f"{int(steps[-1]):,}")],
            steps=steps,
            chart={
                "x": times,
                "series": [frames.max(axis=(1, 2)), frames.mean(axis=(1, 2))],
                "labels": ["Peak", "Plate mean"],
                "xlabel": "t (s)",
                "ylabel": "T (°C)",
            },
        )
    if key == "rotating_blob_advection_2d":
        module = example_module(key)
        params = dict(n=160 if fine else 110, revolutions=1.0, snapshots=181 if fine else 121)
        x, y, times, steps, frames, dt, (u, v) = module.simulate(
            **params, progress_callback=lambda current, total: progress(current / total)
        )
        exact = np.array([module.exact_blob(x, y, t).max() for t in times])
        return Model(
            key,
            "surface",
            x,
            y,
            frames,
            "scalar",
            "DYE CONCENTRATION",
            {**params, "dt": float(dt), "cfl": 0.8},
            "Upwind solution · height: concentration · dots: tracers carried by the flow",
            [("Final peak", f"{frames[-1].max():.2f} (exact 1)"), ("CFL", "0.8")],
            vectors=np.stack((u, v))[None],
            steps=steps,
            chart={
                "x": times,
                "series": [frames.max(axis=(1, 2)), exact],
                "labels": ["Upwind peak", "Exact peak"],
                "xlabel": "t (s)",
                "ylabel": "peak",
            },
        )
    if key == "enhanced_heat_transfer_solver_2d":
        params = dict(
            n=64 if fine else 40, L=1.0, k=50.0, S=50000.0, boundary_temperature=100.0
        )
        x, y, T, iterations = solve_quarter_conduction(**params, tolerance=1e-7)
        qy, qx = np.gradient(T, y, x)
        return Model(
            key,
            "surface",
            x,
            y,
            T[None],
            "°C",
            "TEMPERATURE",
            {**params, "iterations": iterations},
            "Steady solution · normalized height · arrows: heat flux",
            [("Peak temperature", f"{T.max():.1f} °C"), ("Generation", "50 kW/m³")],
            vectors=np.stack((-50 * qx, -50 * qy))[None],
            chart={
                "x": x,
                "series": [T[0], T[len(y) // 2]],
                "labels": ["Insulated edge", "Midsection"],
                "xlabel": "x (m)",
                "ylabel": "T (°C)",
            },
        )
    if key == "heat_convection_solver_2d":
        params = dict(
            nx=72 if fine else 48,
            ny=72 if fine else 48,
            Lx=1.0,
            Ly=1.0,
            Gamma=0.1,
            rho=1.0,
            u=2.0,
            v=1.0,
            phi_left=0.0,
            phi_right=1.0,
            phi_bottom=0.0,
            phi_top=1.0,
        )
        x, y, phi, iterations, px, py = solve_transport(**params, tolerance=1e-8)
        velocity = np.stack((np.full_like(phi, 2.0), np.full_like(phi, 1.0)))[None]
        return Model(
            key,
            "surface",
            x,
            y,
            phi[None],
            "scalar",
            "SCALAR VALUE",
            {**params, "iterations": iterations},
            "Steady scalar · normalized height · prescribed velocity arrows",
            [("Pe x / y", f"{px:.0f} / {py:.0f}"), ("Flow u / v", "2 / 1 m/s")],
            vectors=velocity,
            chart={
                "x": x,
                "series": [phi[len(y) // 2]],
                "labels": ["Midsection"],
                "xlabel": "x (m)",
                "ylabel": "Scalar",
            },
        )
    if key == "heat_conduction_solver_2d":
        nx, ny = (160, 32) if fine else (100, 20)
        x, y = np.linspace(0, 30, nx + 1), np.linspace(0, 1, ny + 1)
        params = dict(
            length=30.0,
            height=1.0,
            nx=nx,
            ny=ny,
            Ti=20.0,
            hf=1.0,
            rho=1.0,
            c=10.0,
            k=1.0,
        )
        T, iterations, residual = solve_channel(
            np.full((ny + 1, nx + 1), 20.0),
            30 / nx,
            1 / ny,
            y,
            20,
            1,
            1,
            10,
            1,
            1,
            1e-8,
        )
        module = example_module(key)
        bulk = module.calculate_bulk_temperature(T, y, 1 / ny, 1, nx)
        _, nu = module.calculate_ht_and_nusselt(bulk, np.ones(nx + 1), 1, 1, T)
        u = module.velocity(y, 1)
        velocity = np.stack((np.broadcast_to(u[:, None], T.shape), np.zeros_like(T)))[
            None
        ]
        return Model(
            key,
            "channel",
            x,
            y,
            T[None],
            "°C",
            "TEMPERATURE",
            {**params, "iterations": iterations, "residual": residual},
            "Steady thermal field · prescribed flow · y expanded for display",
            [("Developed Nu", f"{nu[-10]:.3f}"), ("Mean flow", "1 m/s")],
            vectors=velocity,
            chart={
                "x": x,
                "series": [T[0], bulk],
                "labels": ["Wall", "Bulk"],
                "xlabel": "x (m)",
                "ylabel": "T (°C)",
            },
        )
    module = example_module(key)
    n, L = (160 if fine else 100), 1.0
    x = np.linspace(0, L, n + 1)
    params = dict(n=n, L=L)
    if key == "fixed_end_heat_transfer_1d":
        params.update(k=50.0, S=100000.0, T_left=100.0, T_right=200.0)
        matrix, b = module.setup_problem(**params)
        T = tdma_solve(*matrix, b)
        reference = module.analytical_solution(x, L, 50, 100000, 100, 200)
        qleft, qright, total_flux = module.calculate_energy_balance(T, L / n, 50)
        metrics = [
            ("Peak temperature", f"{T.max():.1f} °C"),
            ("End flux sum", f"{total_flux:.0f} W/m²"),
        ]
        chart = {
            "x": x,
            "series": [T, reference],
            "labels": ["Numerical", "Analytical"],
            "xlabel": "x (m)",
            "ylabel": "T (°C)",
        }
    elif key == "convective_end_heat_transfer_1d":
        params.update(k=50.0, h_l=100.0, h_r=50.0, T_l=100.0, T_r=300.0)
        matrix = module.initialize_matrices(n, L / n, 50, 100, 50, 100, 300)
        T = tdma_solve(*matrix)
        reference = module.analytical_solution(x, L, 50, 100, 50, 100, 300)
        metrics = [
            ("Left surface", f"{T[0]:.0f} °C"),
            ("Right surface", f"{T[-1]:.0f} °C"),
        ]
        chart = {
            "x": x,
            "series": [T, reference],
            "labels": ["Numerical", "Analytical"],
            "xlabel": "x (m)",
            "ylabel": "T (°C)",
        }
    elif key == "enhanced_heat_transfer_solver_1d":
        params.update(k_c=50.0, beta=0.001, h_l=100.0, h_r=50.0, T_l=100.0, T_r=300.0)
        x, T, iterations, conductivity = module.solve_heat_transfer_equation(
            **params, Rtol=1e-11
        )
        params["iterations"] = iterations
        metrics = [
            ("k at cold end", f"{conductivity[0]:.1f} W/mK"),
            ("k at hot end", f"{conductivity[-1]:.1f} W/mK"),
        ]
        chart = {
            "x": x,
            "series": [conductivity],
            "labels": ["Conductivity"],
            "xlabel": "x (m)",
            "ylabel": "k (W/mK)",
        }
    elif key == "heat_convection_solver_1d":
        params.update(rho=1.0, u=2.5, Gamma=0.1, phi_0=0.0, phi_L=1.0, scheme="UDS")
        a, lower, upper = module.setup_matrix(n, L / n, 0.1, 2.5, True)
        a[[0, -1]] = 1
        lower[[0, -1]] = 0
        upper[[0, -1]] = 0
        b = np.zeros(n + 1)
        b[-1] = 1
        T = tdma_solve(a, lower, upper, b)
        reference = module.analytical_solution(x, L, 25, 0, 1)
        metrics = [("Pe", "25"), ("Max error", f"{np.max(abs(T - reference)):.3f}")]
        chart = {
            "x": x,
            "series": [T, reference],
            "labels": ["Upwind", "Analytical"],
            "xlabel": "x (m)",
            "ylabel": "Scalar",
        }
    else:
        raise ValueError(f"No numerical reel adapter for {key}")
    unit = "scalar" if key == "heat_convection_solver_1d" else "°C"
    return Model(
        key,
        "rod",
        x,
        np.array([0.0]),
        T[None, None, :],
        unit,
        "SCALAR VALUE" if unit == "scalar" else "TEMPERATURE",
        params,
        "Steady solution · rod radius is illustrative",
        metrics,
        chart=chart,
    )


def source_fingerprint(key):
    """Invalidate caches whenever the solver or its shared kernels change."""
    root = Path(__file__).resolve().parents[1]
    paths = [
        root / "09_cfd" / f"{key}.py",
        Path(__file__),
        root / "common" / "numerics.py",
        root / "common" / "heat_transfer.py",
    ]
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.read_bytes())
    return digest.hexdigest()


def prepare_model(key, quality="fast", cache_dir=Path(".cache/reels"), progress=None):
    """Cache numerical arrays, never a fabricated animation of a steady field."""
    progress = progress or (lambda fraction: None)
    fingerprint = source_fingerprint(key)
    cache = (
        None
        if cache_dir is None
        else Path(cache_dir) / f"{key}-{quality}-{fingerprint[:16]}.npz"
    )
    if cache is not None and cache.exists():
        try:
            with np.load(cache, allow_pickle=False) as saved:
                meta = json.loads(str(saved["metadata"]))
                arrays = {name: saved[name] for name in ("x", "y", "fields")}
                arrays.update(
                    {
                        name: saved[name] if saved[name].size else None
                        for name in ("vectors", "solid", "steps")
                    }
                )
                chart = meta.pop("chart")
                if chart:
                    chart["x"] = saved["chart_x"]
                    chart["series"] = list(saved["chart_series"])
                return Model(**meta, **arrays, chart=chart, cache_hit=True)
        except (ValueError, OSError, KeyError):
            pass  # A partial/stale cache must never replace computing a field.
    model = _solve(key, quality, progress)
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        meta = {
            name: getattr(model, name)
            for name in (
                "key",
                "kind",
                "unit",
                "field_name",
                "params",
                "note",
                "metrics",
            )
        }
        meta["chart"] = {
            k: v for k, v in model.chart.items() if k not in ("x", "series")
        }
        payload = {
            name: getattr(model, name)
            if getattr(model, name) is not None
            else np.array([])
            for name in ("x", "y", "fields", "vectors", "solid", "steps")
        }
        if model.chart:
            payload.update(
                chart_x=model.chart["x"], chart_series=np.array(model.chart["series"])
            )
        with tempfile.NamedTemporaryFile(
            dir=cache.parent, suffix=".npz", delete=False
        ) as handle:
            temporary = Path(handle.name)
        try:
            np.savez_compressed(temporary, metadata=json.dumps(meta), **payload)
            temporary.replace(cache)
        finally:
            temporary.unlink(missing_ok=True)
    return model
