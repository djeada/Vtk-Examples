"""VTK meshes with the same y-row/x-column convention as the heat solvers."""

import numpy as np
import vtk

from .heat_transfer import validate_field
from .numerics import validate_count


def field_points(x, y, field, name="Temperature", height=0.0):
    """Map each field[y,x] to its physical (x,y), with optional display height.

    Height is a normalized scalar displacement, not a physical z coordinate.
    The attached point scalars always retain the original field values.
    """
    x, y, field = validate_field(x, y, field)
    if not np.isfinite(height) or height < 0:
        raise ValueError("Display height must be finite and nonnegative.")
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    scalars = vtk.vtkDoubleArray()
    scalars.SetName(name)
    span = float(np.ptp(field))
    displacement = (field - field.min()) * height / span if span else field * 0
    for j, py in enumerate(y):
        for i, px in enumerate(x):
            points.InsertNextPoint(px, py, displacement[j, i])
            scalars.InsertNextValue(field[j, i])
    return points, scalars


def quad_field(nx, ny, points, scalars):
    """Connect a row-major grid with counterclockwise, positive-z quads."""
    validate_count("nx", nx, minimum=2)
    validate_count("ny", ny, minimum=2)
    if points.GetNumberOfPoints() != nx * ny or scalars.GetNumberOfTuples() != nx * ny:
        raise ValueError("Point and scalar counts must match the grid dimensions.")
    cells = vtk.vtkCellArray()
    for j in range(ny - 1):
        for i in range(nx - 1):
            ids = (j * nx + i, j * nx + i + 1, (j + 1) * nx + i + 1, (j + 1) * nx + i)
            cells.InsertNextCell(4, ids)
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(cells)
    mesh.GetPointData().SetScalars(scalars)
    return mesh
