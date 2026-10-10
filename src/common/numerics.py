"""Shared numerical helpers for the heat-transfer and flow examples."""

from numbers import Integral

import numpy as np


class ConvergenceError(RuntimeError):
    """Raised when an iterative numerical method reaches its iteration limit."""

    def __init__(self, message: str, iterations: int, residual: float):
        self.iterations = int(iterations)
        self.residual = float(residual)
        super().__init__(
            f"{message} after {self.iterations} iterations "
            f"(residual={self.residual:.6g})"
        )


def _finite_vector(name, values):
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a one-dimensional array of finite numbers") from error
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if array.size == 0:
        raise ValueError(f"{name} must not be empty")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array.copy()


def tdma_solve(main, lower, upper, rhs) -> np.ndarray:
    """Solve a tridiagonal system without modifying the supplied arrays.

    ``lower[i]`` multiplies ``x[i - 1]`` and ``upper[i]`` multiplies
    ``x[i + 1]``. The first lower and last upper entries are unused.
    """
    main = _finite_vector("main", main)
    lower = _finite_vector("lower", lower)
    upper = _finite_vector("upper", upper)
    rhs = _finite_vector("rhs", rhs)

    size = main.size
    if any(vector.size != size for vector in (lower, upper, rhs)):
        raise ValueError("main, lower, upper, and rhs must have the same length")

    scale = max(
        float(np.max(np.abs(main))),
        float(np.max(np.abs(lower))),
        float(np.max(np.abs(upper))),
    )
    pivot_tolerance = np.finfo(float).eps * max(1, size) * scale

    def check_pivot(pivot, index):
        if not np.isfinite(pivot) or abs(pivot) <= pivot_tolerance:
            raise np.linalg.LinAlgError(
                f"Unsafe or zero TDMA pivot at index {index}: {pivot!r}"
            )

    for index in range(1, size):
        check_pivot(main[index - 1], index - 1)
        factor = lower[index] / main[index - 1]
        main[index] -= factor * upper[index - 1]
        rhs[index] -= factor * rhs[index - 1]
        if not np.isfinite(main[index]) or not np.isfinite(rhs[index]):
            raise np.linalg.LinAlgError(
                f"Non-finite value during TDMA elimination at index {index}"
            )

    check_pivot(main[-1], size - 1)
    solution = np.empty(size, dtype=float)
    solution[-1] = rhs[-1] / main[-1]
    for index in range(size - 2, -1, -1):
        check_pivot(main[index], index)
        solution[index] = (rhs[index] - upper[index] * solution[index + 1]) / main[index]
    if not np.all(np.isfinite(solution)):
        raise np.linalg.LinAlgError("TDMA solution contains non-finite values")
    return solution


def validate_positive(**parameters):
    """Require each named scalar parameter to be finite and greater than zero."""
    for name, value in parameters.items():
        if isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} must be finite and greater than zero")
        try:
            numeric_value = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{name} must be finite and greater than zero") from error
        if not np.isfinite(numeric_value) or numeric_value <= 0:
            raise ValueError(f"{name} must be finite and greater than zero")


def validate_finite(**parameters):
    """Require each named scalar parameter to be a finite number."""
    for name, value in parameters.items():
        if isinstance(value, (bool, np.bool_)):
            raise ValueError(f"{name} must be finite")
        try:
            numeric_value = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"{name} must be finite") from error
        if not np.isfinite(numeric_value):
            raise ValueError(f"{name} must be finite")


def validate_count(name, value, minimum=1):
    """Require an integer count (excluding booleans) at least ``minimum``."""
    if isinstance(minimum, bool) or not isinstance(minimum, Integral) or minimum < 0:
        raise ValueError("minimum must be a non-negative integer")
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer greater than or equal to {minimum}")
    return int(value)


def validate_iteration_controls(tolerance, max_iter):
    """Validate and normalize tolerance and iteration-limit arguments."""
    validate_positive(tolerance=tolerance)
    max_iter = validate_count("max_iter", max_iter)
    return float(tolerance), max_iter
