"""
VTK Cell Types Demo: Interactive Visualization of All VTK Cell Types

This module provides an interactive demonstration of all VTK cell types using
a PyQt6-based interface with a combo box for selecting different cell types.

Each cell type is displayed clearly in the center of the view with:
- Visible vertices (as spheres)
- Visible edges for better understanding of topology
- Cell information displayed as text overlay

Cell Types Covered:
-------------------
Linear Cells:
- Vertex (VTK_VERTEX) - 0D: Single point
- PolyVertex (VTK_POLY_VERTEX) - 0D: Multiple points
- Line (VTK_LINE) - 1D: Two-point line segment
- PolyLine (VTK_POLY_LINE) - 1D: Connected line segments
- Triangle (VTK_TRIANGLE) - 2D: Three-point triangle
- TriangleStrip (VTK_TRIANGLE_STRIP) - 2D: Connected triangles
- Polygon (VTK_POLYGON) - 2D: N-sided polygon
- Quad (VTK_QUAD) - 2D: Four-point quadrilateral
- Pixel (VTK_PIXEL) - 2D: Axis-aligned quad (image element)
- Tetra (VTK_TETRA) - 3D: Four-point tetrahedron
- Voxel (VTK_VOXEL) - 3D: Axis-aligned hexahedron
- Hexahedron (VTK_HEXAHEDRON) - 3D: Eight-point hexahedron
- Wedge (VTK_WEDGE) - 3D: Six-point wedge (triangular prism)
- Pyramid (VTK_PYRAMID) - 3D: Five-point pyramid
- PentagonalPrism (VTK_PENTAGONAL_PRISM) - 3D: Ten-point pentagonal prism
- HexagonalPrism (VTK_HEXAGONAL_PRISM) - 3D: Twelve-point hexagonal prism

Quadratic Cells:
- QuadraticEdge (VTK_QUADRATIC_EDGE) - 1D: Three-point quadratic line
- QuadraticTriangle (VTK_QUADRATIC_TRIANGLE) - 2D: Six-point quadratic triangle
- QuadraticQuad (VTK_QUADRATIC_QUAD) - 2D: Eight-point quadratic quad
- QuadraticTetra (VTK_QUADRATIC_TETRA) - 3D: Ten-point quadratic tetrahedron
- QuadraticHexahedron (VTK_QUADRATIC_HEXAHEDRON) - 3D: Twenty-point quadratic hex
- QuadraticWedge (VTK_QUADRATIC_WEDGE) - 3D: Fifteen-point quadratic wedge
- QuadraticPyramid (VTK_QUADRATIC_PYRAMID) - 3D: Thirteen-point quadratic pyramid

Usage:
------
Run this script directly to launch the interactive viewer:
    python cell_types_demo.py

Select different cell types from the combo box to visualize each one.
"""

import math
import sys

import vtk
from PyQt6.QtCore import Qt

# Visual constants for cell rendering
VERTEX_SPHERE_RADIUS = 0.08  # Radius of spheres showing cell vertices
CAMERA_AZIMUTH = 30  # Initial camera azimuth angle in degrees
CAMERA_ELEVATION = 20  # Initial camera elevation angle in degrees
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QVBoxLayout,
    QWidget,
)
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


def create_vertex():
    """Create a single vertex (0D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)

    vertex = vtk.vtkVertex()
    vertex.GetPointIds().SetId(0, 0)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())

    return ugrid, "Vertex (VTK_VERTEX)", "0D cell: A single point in space"


def create_poly_vertex():
    """Create a poly vertex (multiple vertices, 0D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-0.5, -0.5, 0)
    points.InsertNextPoint(0.5, -0.5, 0)
    points.InsertNextPoint(0, 0.5, 0)
    points.InsertNextPoint(0, 0, 0.5)

    poly_vertex = vtk.vtkPolyVertex()
    poly_vertex.GetPointIds().SetNumberOfIds(4)
    for i in range(4):
        poly_vertex.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(poly_vertex.GetCellType(), poly_vertex.GetPointIds())

    return (
        ugrid,
        "PolyVertex (VTK_POLY_VERTEX)",
        "0D cell: Multiple disconnected points",
    )


def create_line():
    """Create a line (1D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1, 0, 0)
    points.InsertNextPoint(1, 0, 0)

    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, 0)
    line.GetPointIds().SetId(1, 1)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(line.GetCellType(), line.GetPointIds())

    return ugrid, "Line (VTK_LINE)", "1D cell: Two-point line segment"


def create_poly_line():
    """Create a poly line (connected line segments, 1D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1, 0, 0)
    points.InsertNextPoint(-0.3, 0.5, 0)
    points.InsertNextPoint(0.3, -0.5, 0)
    points.InsertNextPoint(1, 0, 0)

    poly_line = vtk.vtkPolyLine()
    poly_line.GetPointIds().SetNumberOfIds(4)
    for i in range(4):
        poly_line.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(poly_line.GetCellType(), poly_line.GetPointIds())

    return ugrid, "PolyLine (VTK_POLY_LINE)", "1D cell: Connected line segments"


def create_triangle():
    """Create a triangle (2D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1, -0.577, 0)
    points.InsertNextPoint(1, -0.577, 0)
    points.InsertNextPoint(0, 0.577, 0)

    triangle = vtk.vtkTriangle()
    for i in range(3):
        triangle.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(triangle.GetCellType(), triangle.GetPointIds())

    return ugrid, "Triangle (VTK_TRIANGLE)", "2D cell: Three-point planar triangle"


def create_triangle_strip():
    """Create a triangle strip (connected triangles, 2D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1.5, -0.5, 0)
    points.InsertNextPoint(-1.5, 0.5, 0)
    points.InsertNextPoint(-0.5, -0.5, 0)
    points.InsertNextPoint(-0.5, 0.5, 0)
    points.InsertNextPoint(0.5, -0.5, 0)
    points.InsertNextPoint(0.5, 0.5, 0)
    points.InsertNextPoint(1.5, -0.5, 0)
    points.InsertNextPoint(1.5, 0.5, 0)

    strip = vtk.vtkTriangleStrip()
    strip.GetPointIds().SetNumberOfIds(8)
    for i in range(8):
        strip.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(strip.GetCellType(), strip.GetPointIds())

    return (
        ugrid,
        "TriangleStrip (VTK_TRIANGLE_STRIP)",
        "2D cell: Strip of connected triangles",
    )


def create_polygon():
    """Create a polygon (N-sided, 2D cell) - pentagon."""
    points = vtk.vtkPoints()
    n_sides = 5
    radius = 1.0

    for i in range(n_sides):
        angle = 2 * math.pi * i / n_sides - math.pi / 2
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.InsertNextPoint(x, y, 0)

    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(n_sides)
    for i in range(n_sides):
        polygon.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(polygon.GetCellType(), polygon.GetPointIds())

    return (
        ugrid,
        "Polygon (VTK_POLYGON)",
        "2D cell: N-sided planar polygon (pentagon shown)",
    )


def create_quad():
    """Create a quadrilateral (2D cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(-1, 1, 0)

    quad = vtk.vtkQuad()
    for i in range(4):
        quad.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(quad.GetCellType(), quad.GetPointIds())

    return ugrid, "Quad (VTK_QUAD)", "2D cell: Four-point quadrilateral"


def create_pixel():
    """Create a pixel (axis-aligned quad, 2D cell)."""
    points = vtk.vtkPoints()
    # Pixel uses different ordering: (i,j), (i+1,j), (i,j+1), (i+1,j+1)
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(-1, 1, 0)
    points.InsertNextPoint(1, 1, 0)

    pixel = vtk.vtkPixel()
    for i in range(4):
        pixel.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(pixel.GetCellType(), pixel.GetPointIds())

    return (
        ugrid,
        "Pixel (VTK_PIXEL)",
        "2D cell: Axis-aligned quadrilateral (image element)",
    )


def create_tetra():
    """Create a tetrahedron (3D cell)."""
    points = vtk.vtkPoints()
    # Regular tetrahedron
    points.InsertNextPoint(1, 0, -0.707)
    points.InsertNextPoint(-1, 0, -0.707)
    points.InsertNextPoint(0, 1, 0.707)
    points.InsertNextPoint(0, -1, 0.707)

    tetra = vtk.vtkTetra()
    for i in range(4):
        tetra.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

    return ugrid, "Tetra (VTK_TETRA)", "3D cell: Four-point tetrahedron"


def create_voxel():
    """Create a voxel (axis-aligned hexahedron, 3D cell)."""
    points = vtk.vtkPoints()
    # Voxel ordering: (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k),
    #                 (i,j,k+1), (i+1,j,k+1), (i,j+1,k+1), (i+1,j+1,k+1)
    points.InsertNextPoint(-1, -1, -1)
    points.InsertNextPoint(1, -1, -1)
    points.InsertNextPoint(-1, 1, -1)
    points.InsertNextPoint(1, 1, -1)
    points.InsertNextPoint(-1, -1, 1)
    points.InsertNextPoint(1, -1, 1)
    points.InsertNextPoint(-1, 1, 1)
    points.InsertNextPoint(1, 1, 1)

    voxel = vtk.vtkVoxel()
    for i in range(8):
        voxel.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(voxel.GetCellType(), voxel.GetPointIds())

    return ugrid, "Voxel (VTK_VOXEL)", "3D cell: Axis-aligned hexahedron"


def create_hexahedron():
    """Create a hexahedron (3D cell)."""
    points = vtk.vtkPoints()
    # Hexahedron ordering: bottom face CCW, then top face CCW
    points.InsertNextPoint(-1, -1, -1)
    points.InsertNextPoint(1, -1, -1)
    points.InsertNextPoint(1, 1, -1)
    points.InsertNextPoint(-1, 1, -1)
    points.InsertNextPoint(-1, -1, 1)
    points.InsertNextPoint(1, -1, 1)
    points.InsertNextPoint(1, 1, 1)
    points.InsertNextPoint(-1, 1, 1)

    hexa = vtk.vtkHexahedron()
    for i in range(8):
        hexa.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(hexa.GetCellType(), hexa.GetPointIds())

    return (
        ugrid,
        "Hexahedron (VTK_HEXAHEDRON)",
        "3D cell: Eight-point hexahedron (cube)",
    )


def create_wedge():
    """Create a wedge/triangular prism (3D cell)."""
    points = vtk.vtkPoints()
    # Bottom triangle
    points.InsertNextPoint(-1, 0, -1)
    points.InsertNextPoint(1, 0, -1)
    points.InsertNextPoint(0, 1.5, -1)
    # Top triangle
    points.InsertNextPoint(-1, 0, 1)
    points.InsertNextPoint(1, 0, 1)
    points.InsertNextPoint(0, 1.5, 1)

    wedge = vtk.vtkWedge()
    for i in range(6):
        wedge.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(wedge.GetCellType(), wedge.GetPointIds())

    return ugrid, "Wedge (VTK_WEDGE)", "3D cell: Six-point triangular prism"


def create_pyramid():
    """Create a pyramid (3D cell)."""
    points = vtk.vtkPoints()
    # Base (quad)
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(-1, 1, 0)
    # Apex
    points.InsertNextPoint(0, 0, 1.5)

    pyramid = vtk.vtkPyramid()
    for i in range(5):
        pyramid.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(pyramid.GetCellType(), pyramid.GetPointIds())

    return ugrid, "Pyramid (VTK_PYRAMID)", "3D cell: Five-point pyramid with quad base"


def create_pentagonal_prism():
    """Create a pentagonal prism (3D cell)."""
    points = vtk.vtkPoints()
    n_sides = 5
    radius = 1.0

    # Bottom pentagon
    for i in range(n_sides):
        angle = 2 * math.pi * i / n_sides - math.pi / 2
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.InsertNextPoint(x, y, -1)

    # Top pentagon
    for i in range(n_sides):
        angle = 2 * math.pi * i / n_sides - math.pi / 2
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.InsertNextPoint(x, y, 1)

    prism = vtk.vtkPentagonalPrism()
    for i in range(10):
        prism.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(prism.GetCellType(), prism.GetPointIds())

    return (
        ugrid,
        "PentagonalPrism (VTK_PENTAGONAL_PRISM)",
        "3D cell: Ten-point pentagonal prism",
    )


def create_hexagonal_prism():
    """Create a hexagonal prism (3D cell)."""
    points = vtk.vtkPoints()
    n_sides = 6
    radius = 1.0

    # Bottom hexagon
    for i in range(n_sides):
        angle = 2 * math.pi * i / n_sides
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.InsertNextPoint(x, y, -1)

    # Top hexagon
    for i in range(n_sides):
        angle = 2 * math.pi * i / n_sides
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.InsertNextPoint(x, y, 1)

    prism = vtk.vtkHexagonalPrism()
    for i in range(12):
        prism.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(prism.GetCellType(), prism.GetPointIds())

    return (
        ugrid,
        "HexagonalPrism (VTK_HEXAGONAL_PRISM)",
        "3D cell: Twelve-point hexagonal prism",
    )


def create_quadratic_edge():
    """Create a quadratic edge (1D quadratic cell)."""
    points = vtk.vtkPoints()
    points.InsertNextPoint(-1, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(0, 0.3, 0)  # Mid-edge point (offset to show curvature)

    edge = vtk.vtkQuadraticEdge()
    for i in range(3):
        edge.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(edge.GetCellType(), edge.GetPointIds())

    return (
        ugrid,
        "QuadraticEdge (VTK_QUADRATIC_EDGE)",
        "1D cell: Three-point quadratic line",
    )


def create_quadratic_triangle():
    """Create a quadratic triangle (2D quadratic cell)."""
    points = vtk.vtkPoints()
    # Corner vertices
    points.InsertNextPoint(-1, -0.577, 0)
    points.InsertNextPoint(1, -0.577, 0)
    points.InsertNextPoint(0, 1.155, 0)
    # Mid-edge vertices
    points.InsertNextPoint(0, -0.577, 0)
    points.InsertNextPoint(0.5, 0.289, 0)
    points.InsertNextPoint(-0.5, 0.289, 0)

    tri = vtk.vtkQuadraticTriangle()
    for i in range(6):
        tri.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(tri.GetCellType(), tri.GetPointIds())

    return (
        ugrid,
        "QuadraticTriangle (VTK_QUADRATIC_TRIANGLE)",
        "2D cell: Six-point quadratic triangle",
    )


def create_quadratic_quad():
    """Create a quadratic quadrilateral (2D quadratic cell)."""
    points = vtk.vtkPoints()
    # Corner vertices
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(-1, 1, 0)
    # Mid-edge vertices
    points.InsertNextPoint(0, -1, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(0, 1, 0)
    points.InsertNextPoint(-1, 0, 0)

    quad = vtk.vtkQuadraticQuad()
    for i in range(8):
        quad.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(quad.GetCellType(), quad.GetPointIds())

    return (
        ugrid,
        "QuadraticQuad (VTK_QUADRATIC_QUAD)",
        "2D cell: Eight-point quadratic quadrilateral",
    )


def create_quadratic_tetra():
    """Create a quadratic tetrahedron (3D quadratic cell)."""
    points = vtk.vtkPoints()
    # Corner vertices
    points.InsertNextPoint(1, 0, -0.707)
    points.InsertNextPoint(-1, 0, -0.707)
    points.InsertNextPoint(0, 1, 0.707)
    points.InsertNextPoint(0, -1, 0.707)
    # Mid-edge vertices
    points.InsertNextPoint(0, 0, -0.707)  # Edge 0-1
    points.InsertNextPoint(-0.5, 0.5, 0)  # Edge 1-2
    points.InsertNextPoint(0.5, 0.5, 0)  # Edge 0-2
    points.InsertNextPoint(0.5, -0.5, 0)  # Edge 0-3
    points.InsertNextPoint(-0.5, -0.5, 0)  # Edge 1-3
    points.InsertNextPoint(0, 0, 0.707)  # Edge 2-3

    tetra = vtk.vtkQuadraticTetra()
    for i in range(10):
        tetra.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

    return (
        ugrid,
        "QuadraticTetra (VTK_QUADRATIC_TETRA)",
        "3D cell: Ten-point quadratic tetrahedron",
    )


def create_quadratic_hexahedron():
    """Create a quadratic hexahedron (3D quadratic cell)."""
    points = vtk.vtkPoints()
    # 8 corner vertices
    points.InsertNextPoint(-1, -1, -1)
    points.InsertNextPoint(1, -1, -1)
    points.InsertNextPoint(1, 1, -1)
    points.InsertNextPoint(-1, 1, -1)
    points.InsertNextPoint(-1, -1, 1)
    points.InsertNextPoint(1, -1, 1)
    points.InsertNextPoint(1, 1, 1)
    points.InsertNextPoint(-1, 1, 1)
    # 12 mid-edge vertices
    points.InsertNextPoint(0, -1, -1)
    points.InsertNextPoint(1, 0, -1)
    points.InsertNextPoint(0, 1, -1)
    points.InsertNextPoint(-1, 0, -1)
    points.InsertNextPoint(0, -1, 1)
    points.InsertNextPoint(1, 0, 1)
    points.InsertNextPoint(0, 1, 1)
    points.InsertNextPoint(-1, 0, 1)
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(-1, 1, 0)

    hexa = vtk.vtkQuadraticHexahedron()
    for i in range(20):
        hexa.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(hexa.GetCellType(), hexa.GetPointIds())

    return (
        ugrid,
        "QuadraticHexahedron (VTK_QUADRATIC_HEXAHEDRON)",
        "3D cell: Twenty-point quadratic hexahedron",
    )


def create_quadratic_wedge():
    """Create a quadratic wedge (3D quadratic cell)."""
    points = vtk.vtkPoints()
    # 6 corner vertices
    points.InsertNextPoint(-1, 0, -1)
    points.InsertNextPoint(1, 0, -1)
    points.InsertNextPoint(0, 1.5, -1)
    points.InsertNextPoint(-1, 0, 1)
    points.InsertNextPoint(1, 0, 1)
    points.InsertNextPoint(0, 1.5, 1)
    # 9 mid-edge vertices
    points.InsertNextPoint(0, 0, -1)  # Edge 0-1
    points.InsertNextPoint(0.5, 0.75, -1)  # Edge 1-2
    points.InsertNextPoint(-0.5, 0.75, -1)  # Edge 2-0
    points.InsertNextPoint(0, 0, 1)  # Edge 3-4
    points.InsertNextPoint(0.5, 0.75, 1)  # Edge 4-5
    points.InsertNextPoint(-0.5, 0.75, 1)  # Edge 5-3
    points.InsertNextPoint(-1, 0, 0)  # Edge 0-3
    points.InsertNextPoint(1, 0, 0)  # Edge 1-4
    points.InsertNextPoint(0, 1.5, 0)  # Edge 2-5

    wedge = vtk.vtkQuadraticWedge()
    for i in range(15):
        wedge.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(wedge.GetCellType(), wedge.GetPointIds())

    return (
        ugrid,
        "QuadraticWedge (VTK_QUADRATIC_WEDGE)",
        "3D cell: Fifteen-point quadratic wedge",
    )


def create_quadratic_pyramid():
    """Create a quadratic pyramid (3D quadratic cell)."""
    points = vtk.vtkPoints()
    # 5 corner vertices
    points.InsertNextPoint(-1, -1, 0)
    points.InsertNextPoint(1, -1, 0)
    points.InsertNextPoint(1, 1, 0)
    points.InsertNextPoint(-1, 1, 0)
    points.InsertNextPoint(0, 0, 1.5)
    # 8 mid-edge vertices (base edges then apex edges)
    points.InsertNextPoint(0, -1, 0)  # Edge 0-1
    points.InsertNextPoint(1, 0, 0)  # Edge 1-2
    points.InsertNextPoint(0, 1, 0)  # Edge 2-3
    points.InsertNextPoint(-1, 0, 0)  # Edge 3-0
    points.InsertNextPoint(-0.5, -0.5, 0.75)  # Edge 0-4
    points.InsertNextPoint(0.5, -0.5, 0.75)  # Edge 1-4
    points.InsertNextPoint(0.5, 0.5, 0.75)  # Edge 2-4
    points.InsertNextPoint(-0.5, 0.5, 0.75)  # Edge 3-4

    pyramid = vtk.vtkQuadraticPyramid()
    for i in range(13):
        pyramid.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(pyramid.GetCellType(), pyramid.GetPointIds())

    return (
        ugrid,
        "QuadraticPyramid (VTK_QUADRATIC_PYRAMID)",
        "3D cell: Thirteen-point quadratic pyramid",
    )


# Dictionary mapping cell type names to their creation functions
CELL_CREATORS = {
    "Vertex": create_vertex,
    "PolyVertex": create_poly_vertex,
    "Line": create_line,
    "PolyLine": create_poly_line,
    "Triangle": create_triangle,
    "TriangleStrip": create_triangle_strip,
    "Polygon (Pentagon)": create_polygon,
    "Quad": create_quad,
    "Pixel": create_pixel,
    "Tetra": create_tetra,
    "Voxel": create_voxel,
    "Hexahedron": create_hexahedron,
    "Wedge": create_wedge,
    "Pyramid": create_pyramid,
    "PentagonalPrism": create_pentagonal_prism,
    "HexagonalPrism": create_hexagonal_prism,
    "QuadraticEdge": create_quadratic_edge,
    "QuadraticTriangle": create_quadratic_triangle,
    "QuadraticQuad": create_quadratic_quad,
    "QuadraticTetra": create_quadratic_tetra,
    "QuadraticHexahedron": create_quadratic_hexahedron,
    "QuadraticWedge": create_quadratic_wedge,
    "QuadraticPyramid": create_quadratic_pyramid,
}


class CellTypesDemo(QMainWindow):
    """
    Interactive VTK Cell Types Demo using PyQt6.

    This widget displays a combo box for selecting different VTK cell types
    and renders the selected cell in a VTK view with visible vertices and edges.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("VTK Cell Types Demo")
        self.resize(900, 700)

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create header with combo box
        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Select Cell Type:"))

        self.combo_box = QComboBox()
        self.combo_box.addItems(CELL_CREATORS.keys())
        self.combo_box.setMinimumWidth(200)
        self.combo_box.currentTextChanged.connect(self.on_cell_type_changed)
        header_layout.addWidget(self.combo_box)

        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Create info label
        self.info_label = QLabel()
        self.info_label.setStyleSheet(
            "QLabel { background-color: #2b2b2b; color: #ffffff; padding: 10px; "
            "border-radius: 5px; font-size: 12px; }"
        )
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.info_label)

        # Create VTK widget
        self.vtk_widget = QVTKRenderWindowInteractor(central_widget)
        layout.addWidget(self.vtk_widget, stretch=1)

        # Set up renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.1, 0.1, 0.15)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)

        # Initialize actors
        self.surface_actor = vtk.vtkActor()
        self.point_actor = vtk.vtkActor()
        self.renderer.AddActor(self.surface_actor)
        self.renderer.AddActor(self.point_actor)

        # Add axes widget
        self.axes = vtk.vtkAxesActor()
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(self.axes)
        self.axes_widget.SetInteractor(self.vtk_widget)
        self.axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
        self.axes_widget.SetEnabled(1)
        self.axes_widget.InteractiveOff()

        # Initialize interactor
        self.vtk_widget.Initialize()

        # Display first cell type
        self.on_cell_type_changed(self.combo_box.currentText())

    def on_cell_type_changed(self, cell_type_name):
        """Handle combo box selection change."""
        if cell_type_name not in CELL_CREATORS:
            return

        # Create the selected cell type
        ugrid, display_name, description = CELL_CREATORS[cell_type_name]()

        # Update info label
        num_points = ugrid.GetNumberOfPoints()
        cell = ugrid.GetCell(0)
        num_edges = cell.GetNumberOfEdges()
        num_faces = cell.GetNumberOfFaces()

        info_text = (
            f"<b>{display_name}</b><br/>"
            f"{description}<br/>"
            f"Points: {num_points} | Edges: {num_edges} | Faces: {num_faces}"
        )
        self.info_label.setText(info_text)

        # Create surface mapper
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(ugrid)

        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)
        self.surface_actor.GetProperty().SetOpacity(0.8)
        self.surface_actor.GetProperty().SetEdgeVisibility(1)
        self.surface_actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
        self.surface_actor.GetProperty().SetLineWidth(2)

        # Create point glyph visualization
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(VERTEX_SPHERE_RADIUS)
        sphere_source.SetThetaResolution(16)
        sphere_source.SetPhiResolution(16)

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(ugrid)
        glyph.SetSourceConnection(sphere_source.GetOutputPort())
        glyph.SetScaleModeToDataScalingOff()

        point_mapper = vtk.vtkPolyDataMapper()
        point_mapper.SetInputConnection(glyph.GetOutputPort())

        self.point_actor.SetMapper(point_mapper)
        self.point_actor.GetProperty().SetColor(0.9, 0.2, 0.2)

        # Reset camera to show the cell clearly
        self.renderer.ResetCamera()
        camera = self.renderer.GetActiveCamera()
        camera.Azimuth(CAMERA_AZIMUTH)
        camera.Elevation(CAMERA_ELEVATION)
        self.renderer.ResetCameraClippingRange()

        # Render the scene
        self.vtk_widget.GetRenderWindow().Render()

    def closeEvent(self, event):
        """Clean up VTK resources on close."""
        self.vtk_widget.Finalize()
        super().closeEvent(event)


def main():
    """Launch the VTK Cell Types Demo application."""
    app = QApplication(sys.argv)
    window = CellTypesDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
