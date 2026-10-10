"""NumPy-only finite-volume kernels used by the 2D teaching examples.

Fields use ``field[y_index, x_index]``. Diffusion and convection coefficients
include the area of the corresponding face (per unit depth). Boundary nodes
with a specified flux use half control volumes rather than copied temperatures.
"""

import numpy as np

from .numerics import (
    ConvergenceError,
    validate_count,
    validate_iteration_controls,
    validate_positive,
)


def validate_field(x, y, field):
    """Validate coordinates and a scalar field before computation or display."""
    x, y, field = (np.asarray(value, dtype=float) for value in (x, y, field))
    if x.ndim != 1 or y.ndim != 1 or len(x) < 2 or len(y) < 2:
        raise ValueError("x and y must be coordinate vectors with at least two points.")
    if field.shape != (len(y), len(x)):
        raise ValueError("Field shape must be (len(y), len(x)).")
    if not all(np.all(np.isfinite(value)) for value in (x, y, field)):
        raise ValueError("Coordinates and field values must be finite.")
    if np.any(np.diff(x) <= 0) or np.any(np.diff(y) <= 0):
        raise ValueError("Coordinates must be strictly increasing.")
    return x, y, field


def _relax(field, region, west, east, south, north, source, tolerance, max_iter):
    """Red-black Gauss-Seidel with a checked, diagonal-scaled equation residual.

    All coefficients have the shape of the unknown region. Fixed boundary
    values surround that region. Zero coefficients remove flux boundary faces.
    A residual in field units is used, so tolerance does not change with the
    arbitrary scale of the conductances. Failure to converge raises explicitly.
    """
    validate_iteration_controls(tolerance, max_iter)
    diagonal = west + east + south + north
    if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
        raise ValueError("Every control volume must have a finite positive diagonal.")
    row, column = region
    y0, y1, x0, x1 = row.start, row.stop, column.start, column.stop
    padded = np.pad(field, 1)
    unknown = padded[y0 + 1 : y1 + 1, x0 + 1 : x1 + 1]
    neighbors = (
        padded[y0 + 1 : y1 + 1, x0:x1],
        padded[y0 + 1 : y1 + 1, x0 + 2 : x1 + 2],
        padded[y0:y1, x0 + 1 : x1 + 1],
        padded[y0 + 2 : y1 + 2, x0 + 1 : x1 + 1],
    )
    checkerboard = np.indices(unknown.shape).sum(axis=0) % 2

    def right_hand_side():
        return source + sum(
            coefficient * neighbor
            for coefficient, neighbor in zip((west, east, south, north), neighbors)
        )

    for iteration in range(max_iter + 1):
        residual = float(np.max(np.abs(right_hand_side() / diagonal - unknown)))
        if not np.isfinite(residual):
            raise ConvergenceError(
                "Nonfinite heat-transfer residual", iteration, residual
            )
        if residual <= tolerance:
            return padded[1:-1, 1:-1].copy(), iteration, residual
        if iteration == max_iter:
            break
        for color in (0, 1):
            values = right_hand_side() / diagonal
            mask = checkerboard == color
            unknown[mask] = values[mask]
    raise ConvergenceError("Heat-transfer iteration limit reached", max_iter, residual)


def _boundary_values(value, coordinates, name):
    """Accept a constant, samples, or a callable for a Dirichlet boundary."""
    if callable(value):
        value = value(coordinates)
    try:
        values = np.broadcast_to(np.asarray(value, dtype=float), coordinates.shape)
    except ValueError as error:
        raise ValueError(
            f"{name} boundary must match its coordinate vector."
        ) from error
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} boundary values must be finite.")
    return values


def solve_transport(
    nx,
    ny,
    Lx,
    Ly,
    Gamma,
    rho,
    u,
    v,
    phi_left,
    phi_right,
    phi_bottom,
    phi_top,
    tolerance=1e-6,
    max_iter=50000,
):
    """Solve rho*(u*phi_x+v*phi_y)=Gamma*Laplacian(phi) with fixed edges.

    Constant velocity, positive diffusion, and first-order upwind convection
    are used. Top/bottom values take precedence at conflicting corner nodes;
    those corners do not enter the interior five-point stencil.
    """
    validate_count("nx", nx, minimum=2)
    validate_count("ny", ny, minimum=2)
    validate_positive(Lx=Lx, Ly=Ly, Gamma=Gamma, rho=rho)
    if not np.all(np.isfinite([u, v])):
        raise ValueError("Velocities must be finite.")
    x, y = np.linspace(0, Lx, nx + 1), np.linspace(0, Ly, ny + 1)
    edges = (
        _boundary_values(phi_left, y, "Left"),
        _boundary_values(phi_right, y, "Right"),
        _boundary_values(phi_bottom, x, "Bottom"),
        _boundary_values(phi_top, x, "Top"),
    )
    field = np.full((ny + 1, nx + 1), np.mean(np.concatenate(edges)))
    field[:, 0], field[:, -1], field[0, :], field[-1, :] = edges
    dx, dy = Lx / nx, Ly / ny
    # Face-integrated fluxes: vertical faces have length dy; horizontal dx.
    Fx, Fy = rho * u * dy, rho * v * dx
    Dx, Dy = Gamma * dy / dx, Gamma * dx / dy
    shape = (ny - 1, nx - 1)
    coefficients = [
        np.full(shape, Dx + max(Fx, 0)),
        np.full(shape, Dx + max(-Fx, 0)),
        np.full(shape, Dy + max(Fy, 0)),
        np.full(shape, Dy + max(-Fy, 0)),
    ]
    field, iterations, _ = _relax(
        field,
        (slice(1, ny), slice(1, nx)),
        *coefficients,
        np.zeros(shape),
        tolerance,
        max_iter,
    )
    return x, y, field, iterations, abs(rho * u) * Lx / Gamma, abs(rho * v) * Ly / Gamma


def solve_quarter_conduction(
    n, L, k, S, boundary_temperature, tolerance, max_iter=50000
):
    """Solve k*Laplacian(T)+S=0; insulated left/bottom, fixed top/right.

    Half control volumes at each insulated edge and a quarter volume at their
    shared corner preserve flux balances and give second-order symmetry stencils.
    """
    validate_count("n", n, minimum=2)
    validate_positive(L=L, k=k)
    if not np.all(np.isfinite([S, boundary_temperature])):
        raise ValueError("Source and boundary temperature must be finite.")
    x = np.linspace(0, L, n + 1)
    dx = L / n
    widths = np.full(n, dx)
    widths[0] *= 0.5
    west = np.broadcast_to(k * widths[:, None] / dx, (n, n)).copy()
    east = west.copy()
    south = np.broadcast_to(k * widths[None, :] / dx, (n, n)).copy()
    north = south.copy()
    west[:, 0], south[0, :] = 0, 0
    source = S * widths[:, None] * widths[None, :]
    field = np.full((n + 1, n + 1), float(boundary_temperature))
    field, iterations, _ = _relax(
        field,
        (slice(0, n), slice(0, n)),
        west,
        east,
        south,
        north,
        source,
        tolerance,
        max_iter,
    )
    return x, x.copy(), field, iterations


def channel_velocity(y, height):
    """Fully developed parallel-plate velocity, mean speed 1 m/s, y in [0,h]."""
    validate_positive(height=height)
    y = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y)) or np.any(y < 0) or np.any(y > height):
        raise ValueError("Channel coordinates must lie between 0 and height.")
    eta = y / height
    return 6 * eta * (1 - eta)


def channel_weights(y, height):
    """Exact integrals of velocity and height over each nodal control volume."""
    validate_positive(height=height)
    y = np.asarray(y, dtype=float)
    if (
        y.ndim != 1
        or len(y) < 2
        or not np.all(np.isfinite(y))
        or np.any(np.diff(y) <= 0)
        or y[0] != 0
        or y[-1] != height
    ):
        raise ValueError("Channel coordinates must increase from 0 to height.")
    faces = np.r_[0, (y[:-1] + y[1:]) / 2, height]
    eta = faces / height
    integral = height * (3 * eta**2 - 2 * eta**3)
    return np.diff(integral), np.diff(faces)


def solve_channel(
    initial,
    dx,
    dy,
    y,
    Ti,
    heat_flux,
    rho,
    capacity,
    k,
    height,
    tolerance,
    max_iter=50000,
):
    """Steady thermal development in prescribed laminar parallel-plate flow.

    The velocity is prescribed, not a solved Navier-Stokes field. Both walls
    receive heat_flux W/m² (positive into the fluid); the inlet is Ti. At the
    outlet, conduction through the boundary face is zero while advective heat
    exits. Wall and outlet half control volumes enforce those flux balances.
    """
    validate_positive(dx=dx, dy=dy, rho=rho, capacity=capacity, k=k, height=height)
    if not np.all(np.isfinite([Ti, heat_flux])):
        raise ValueError("Inlet temperature and heat flux must be finite.")
    field = np.asarray(initial, dtype=float).copy()
    if field.ndim != 2 or min(field.shape) < 3 or not np.all(np.isfinite(field)):
        raise ValueError("Initial field must be finite and have at least 3x3 nodes.")
    ny, nx = np.array(field.shape) - 1
    y = np.asarray(y, dtype=float)
    if (
        y.shape != (ny + 1,)
        or not np.allclose(y, np.linspace(0, height, ny + 1), rtol=1e-12, atol=1e-12)
        or not np.isclose(ny * dy, height, rtol=1e-12, atol=1e-12)
    ):
        raise ValueError(
            "y, dy, height, and temperature rows must describe the same grid."
        )
    velocity_integrals, heights = channel_weights(y, height)
    widths = np.full(nx, dx)
    widths[-1] *= 0.5
    shape = (ny + 1, nx)
    diffusion_x = np.broadcast_to(k * heights[:, None] / dx, shape).copy()
    flux = rho * capacity * velocity_integrals[:, None]
    west = diffusion_x + flux
    east = diffusion_x.copy()
    east[:, -1] = 0
    south = np.broadcast_to(k * widths[None, :] / dy, shape).copy()
    north = south.copy()
    south[0, :], north[-1, :] = 0, 0
    source = np.zeros(shape)
    source[0, :] = source[-1, :] = heat_flux * widths
    field[:, 0] = Ti
    return _relax(
        field,
        (slice(0, ny + 1), slice(1, nx + 1)),
        west,
        east,
        south,
        north,
        source,
        tolerance,
        max_iter,
    )
