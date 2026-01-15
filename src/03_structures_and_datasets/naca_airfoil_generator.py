"""
NACA Airfoil Generator - Comprehensive Module for Aerodynamic Profile Generation

This module provides a comprehensive toolkit for generating, visualizing, and exporting
NACA airfoil profiles for computational fluid dynamics (CFD) and aerodynamic analysis.

Supported Airfoil Types:
------------------------
1. NACA 4-digit (e.g., NACA 0012, NACA 2412, NACA 4415)
   - First digit (M): Maximum camber as percentage of chord
   - Second digit (P): Position of maximum camber (tenths of chord)
   - Last two digits (T): Maximum thickness as percentage of chord

2. NACA 5-digit (e.g., NACA 23012, NACA 24112)
   - First digit: Design lift coefficient * 3/20
   - Second digit: Location of maximum camber (2*second digit = position in tenths)
   - Third digit: 0 = standard, 1 = reflex camber line
   - Last two digits: Maximum thickness as percentage of chord

Features:
---------
- 2D profile generation with parametric control
- 3D wing extrusion (constant, tapered, swept wings)
- Interactive VTK visualization (2D and 3D)
- Matplotlib 2D plotting with annotations
- Multiple export formats (VTP, STL, OBJ, CSV, DAT)
- Pressure coefficient (Cp) visualization
- Surface area and property calculations
- Command-line interface for batch processing

Usage:
------
    # Create a NACA 2412 airfoil
    generator = NACAGenerator("2412")

    # Get 2D profile points
    x, y_upper, y_lower = generator.get_profile_points(num_points=100)

    # Create VTK PolyData
    polydata = generator.create_polydata_2d()

    # Create 3D wing
    wing = generator.create_wing_3d(span=5.0, taper_ratio=0.5, sweep_angle=15.0)

    # Export to various formats
    generator.export("airfoil.vtp", format="vtp")
    generator.export("airfoil.stl", format="stl")

    # Visualize
    generator.visualize_2d_matplotlib()
    generator.visualize_3d_vtk()

References:
-----------
- NACA Report No. 460 (1933) - 4-digit airfoils
- NACA Report No. 537 (1935) - 5-digit airfoils
- https://en.wikipedia.org/wiki/NACA_airfoil
"""

import argparse
import math
import os
import sys
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple

import numpy as np

try:
    import vtk

    HAS_VTK = True
except ImportError:
    HAS_VTK = False

try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


class AirfoilType(Enum):
    """Enumeration of supported NACA airfoil types."""

    NACA_4_DIGIT = "4-digit"
    NACA_5_DIGIT = "5-digit"


class ExportFormat(Enum):
    """Supported export file formats."""

    VTP = "vtp"  # VTK PolyData XML format
    STL = "stl"  # Stereolithography format
    OBJ = "obj"  # Wavefront OBJ format
    CSV = "csv"  # Comma-separated values
    DAT = "dat"  # XFOIL/XFLR5 format
    PLY = "ply"  # Polygon file format


# NACA 4-digit thickness distribution coefficients
# Reference: NACA Report No. 460
NACA_COEFF = {
    "a0": 0.2969,  # sqrt(x) coefficient
    "a1": 0.1260,  # x coefficient
    "a2": 0.3516,  # x^2 coefficient
    "a3": 0.2843,  # x^3 coefficient
    "a4": 0.1015,  # x^4 coefficient (use 0.1036 for closed trailing edge)
}

# NACA 5-digit camber line coefficients
# Reference: NACA Report No. 537
NACA_5_DIGIT_COEFFS = {
    210: {"r": 0.0580, "k1": 361.4},
    220: {"r": 0.1260, "k1": 51.64},
    230: {"r": 0.2025, "k1": 15.957},
    240: {"r": 0.2900, "k1": 6.643},
    250: {"r": 0.3910, "k1": 3.230},
}


@dataclass
class AirfoilParameters:
    """
    Data class for storing NACA airfoil parameters.

    Attributes:
        naca_code: The NACA designation string (e.g., "2412", "23012")
        airfoil_type: Type of NACA airfoil (4-digit or 5-digit)
        max_camber: Maximum camber as fraction of chord (M)
        camber_position: Position of maximum camber as fraction of chord (P)
        thickness: Maximum thickness as fraction of chord (T)
        chord: Chord length (default 1.0)
        is_symmetric: Whether the airfoil is symmetric (M=0)
        is_reflex: For 5-digit airfoils, whether it has reflex camber
        design_cl: Design lift coefficient (5-digit only)
    """

    naca_code: str
    airfoil_type: AirfoilType
    max_camber: float
    camber_position: float
    thickness: float
    chord: float = 1.0
    is_symmetric: bool = False
    is_reflex: bool = False
    design_cl: float = 0.0

    def __post_init__(self):
        self.is_symmetric = self.max_camber == 0


@dataclass
class WingParameters:
    """
    Data class for 3D wing geometry parameters.

    Attributes:
        span: Wing span (tip to tip)
        taper_ratio: Ratio of tip chord to root chord (0 to 1)
        sweep_angle: Leading edge sweep angle in degrees
        dihedral_angle: Dihedral angle in degrees
        twist_angle: Tip twist angle in degrees (washout)
        num_sections: Number of spanwise sections for mesh
    """

    span: float = 10.0
    taper_ratio: float = 1.0
    sweep_angle: float = 0.0
    dihedral_angle: float = 0.0
    twist_angle: float = 0.0
    num_sections: int = 20


class NACAGenerator:
    """
    Comprehensive NACA Airfoil Generator.

    This class provides methods for generating NACA 4-digit and 5-digit airfoil
    profiles, creating 2D and 3D geometries, and exporting to various formats.

    Example:
        >>> generator = NACAGenerator("2412", chord=1.0, num_points=100)
        >>> x, y_upper, y_lower = generator.get_profile_points()
        >>> generator.visualize_2d_matplotlib()
        >>> generator.export("airfoil.vtp")
    """

    def __init__(
        self,
        naca_code: str,
        chord: float = 1.0,
        num_points: int = 100,
        cosine_spacing: bool = True,
        closed_te: bool = False,
    ):
        """
        Initialize the NACA airfoil generator.

        Args:
            naca_code: NACA designation (e.g., "0012", "2412", "23012")
            chord: Chord length (default 1.0)
            num_points: Number of points along each surface
            cosine_spacing: Use cosine spacing for better LE/TE resolution
            closed_te: Close trailing edge to zero thickness
        """
        self.naca_code = naca_code.strip()

        if chord <= 0:
            raise ValueError("Chord length must be positive")
        self.chord = chord

        if num_points < 2:
            raise ValueError("Number of points must be at least 2")
        self.num_points = num_points

        self.cosine_spacing = cosine_spacing
        self.closed_te = closed_te

        # Parse NACA code and extract parameters
        self.params = self._parse_naca_code()

        # Pre-computed profile data
        self._x_coords: Optional[np.ndarray] = None
        self._y_upper: Optional[np.ndarray] = None
        self._y_lower: Optional[np.ndarray] = None
        self._y_camber: Optional[np.ndarray] = None
        self._thickness: Optional[np.ndarray] = None
        self._x_upper: Optional[np.ndarray] = None
        self._x_lower: Optional[np.ndarray] = None

        # Generate profile on initialization
        self._generate_profile()

    def _parse_naca_code(self) -> AirfoilParameters:
        """
        Parse NACA code and extract airfoil parameters.

        Returns:
            AirfoilParameters: Parsed airfoil parameters

        Raises:
            ValueError: If NACA code is invalid
        """
        code = self.naca_code.upper()

        # Remove common prefixes (all uppercase for comparison)
        prefixes = ["NACA-", "NACA ", "NACA"]
        for prefix in prefixes:
            if code.startswith(prefix):
                code = code[len(prefix) :]
                break

        code = code.strip()

        # Update stored naca_code to the clean version
        self.naca_code = code

        if len(code) == 4 and code.isdigit():
            return self._parse_4_digit(code)
        elif len(code) == 5 and code.isdigit():
            return self._parse_5_digit(code)
        else:
            raise ValueError(
                f"Invalid NACA code: '{self.naca_code}'. "
                "Expected 4-digit (e.g., '0012', '2412') or "
                "5-digit (e.g., '23012') format."
            )

    def _parse_4_digit(self, code: str) -> AirfoilParameters:
        """
        Parse NACA 4-digit airfoil code.

        NACA MPTT format:
        - M: Maximum camber (percentage of chord)
        - P: Position of maximum camber (tenths of chord)
        - TT: Maximum thickness (percentage of chord)
        """
        m = int(code[0]) / 100.0  # Max camber
        p = int(code[1]) / 10.0  # Camber position
        t = int(code[2:4]) / 100.0  # Thickness

        return AirfoilParameters(
            naca_code=code,
            airfoil_type=AirfoilType.NACA_4_DIGIT,
            max_camber=m,
            camber_position=p,
            thickness=t,
            chord=self.chord,
        )

    def _parse_5_digit(self, code: str) -> AirfoilParameters:
        """
        Parse NACA 5-digit airfoil code.

        NACA LPSTT format:
        - L: Design lift coefficient = L * (3/20)
        - P: Position of maximum camber = P/2 (tenths of chord)
        - S: 0 = standard, 1 = reflex
        - TT: Maximum thickness (percentage of chord)
        """
        l_digit = int(code[0])
        p_digit = int(code[1])
        s_digit = int(code[2])
        t = int(code[3:5]) / 100.0

        cl = l_digit * 3.0 / 20.0
        p = p_digit / 20.0
        is_reflex = s_digit == 1

        return AirfoilParameters(
            naca_code=code,
            airfoil_type=AirfoilType.NACA_5_DIGIT,
            max_camber=cl,  # Store design Cl
            camber_position=p,
            thickness=t,
            chord=self.chord,
            is_reflex=is_reflex,
            design_cl=cl,
        )

    def _generate_x_spacing(self) -> np.ndarray:
        """
        Generate x-coordinates with optional cosine spacing.

        Cosine spacing provides better resolution at leading and trailing edges
        where curvature is highest.

        Returns:
            np.ndarray: Array of x-coordinates from 0 to chord
        """
        if self.cosine_spacing:
            # Cosine spacing: more points at LE and TE
            beta = np.linspace(0, np.pi, self.num_points)
            x = self.chord * (1 - np.cos(beta)) / 2
        else:
            # Uniform spacing
            x = np.linspace(0, self.chord, self.num_points)
        return x

    def _compute_thickness(self, x: np.ndarray) -> np.ndarray:
        """
        Compute NACA thickness distribution.

        Uses the standard NACA thickness equation:
        yt = 5*t*c * (a0*sqrt(x/c) - a1*(x/c) - a2*(x/c)^2 + a3*(x/c)^3 - a4*(x/c)^4)

        Args:
            x: Array of x-coordinates

        Returns:
            np.ndarray: Half-thickness values at each x
        """
        t = self.params.thickness
        c = self.chord
        x_norm = np.clip(x / c, 0, 1)

        # Trailing edge coefficient (0.1015 for open TE, 0.1036 for closed)
        a4 = 0.1036 if self.closed_te else NACA_COEFF["a4"]

        yt = (
            5
            * t
            * c
            * (
                NACA_COEFF["a0"] * np.sqrt(x_norm)
                - NACA_COEFF["a1"] * x_norm
                - NACA_COEFF["a2"] * x_norm**2
                + NACA_COEFF["a3"] * x_norm**3
                - a4 * x_norm**4
            )
        )

        return yt

    def _compute_camber_4_digit(
        self, x: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute camber line and gradient for NACA 4-digit airfoils.

        Args:
            x: Array of x-coordinates

        Returns:
            Tuple of (camber values, camber line gradient)
        """
        m = self.params.max_camber
        p = self.params.camber_position
        c = self.chord
        x_norm = x / c

        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)

        if m == 0 or p == 0:
            # Symmetric airfoil
            return yc, dyc_dx

        # Front section (x < p*c)
        front = x_norm <= p
        if np.any(front):
            yc[front] = (
                (m / p**2) * (2 * p * x_norm[front] - x_norm[front] ** 2) * c
            )
            dyc_dx[front] = (2 * m / p**2) * (p - x_norm[front])

        # Rear section (x > p*c)
        rear = x_norm > p
        if np.any(rear):
            yc[rear] = (
                (m / (1 - p) ** 2)
                * (1 - 2 * p + 2 * p * x_norm[rear] - x_norm[rear] ** 2)
                * c
            )
            dyc_dx[rear] = (2 * m / (1 - p) ** 2) * (p - x_norm[rear])

        return yc, dyc_dx

    def _compute_camber_5_digit(
        self, x: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute camber line and gradient for NACA 5-digit airfoils.

        Args:
            x: Array of x-coordinates

        Returns:
            Tuple of (camber values, camber line gradient)
        """
        c = self.chord
        x_norm = x / c

        # Get coefficients for this camber line
        key = int(self.naca_code[0:2] + "0")
        if key not in NACA_5_DIGIT_COEFFS:
            key = 230  # Default to most common

        r = NACA_5_DIGIT_COEFFS[key]["r"]
        k1 = NACA_5_DIGIT_COEFFS[key]["k1"]

        yc = np.zeros_like(x)
        dyc_dx = np.zeros_like(x)

        # Front section
        front = x_norm <= r
        if np.any(front):
            xf = x_norm[front]
            yc[front] = (
                (k1 / 6) * (xf**3 - 3 * r * xf**2 + r**2 * (3 - r) * xf) * c
            )
            dyc_dx[front] = (k1 / 6) * (
                3 * xf**2 - 6 * r * xf + r**2 * (3 - r)
            )

        # Rear section
        rear = x_norm > r
        if np.any(rear):
            yc[rear] = (k1 * r**3 / 6) * (1 - x_norm[rear]) * c
            dyc_dx[rear] = -(k1 * r**3 / 6)

        return yc, dyc_dx

    def _generate_profile(self):
        """Generate the complete airfoil profile."""
        x = self._generate_x_spacing()
        yt = self._compute_thickness(x)

        if self.params.airfoil_type == AirfoilType.NACA_4_DIGIT:
            yc, dyc_dx = self._compute_camber_4_digit(x)
        else:
            yc, dyc_dx = self._compute_camber_5_digit(x)

        # Compute surface coordinates
        theta = np.arctan(dyc_dx)

        self._x_coords = x
        self._y_camber = yc
        self._thickness = yt

        # Upper surface
        self._y_upper = yc + yt * np.cos(theta)
        self._x_upper = x - yt * np.sin(theta)

        # Lower surface
        self._y_lower = yc - yt * np.cos(theta)
        self._x_lower = x + yt * np.sin(theta)

    def get_profile_points(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get the airfoil profile coordinates.

        Returns:
            Tuple of (x_coords, y_upper, y_lower)
        """
        return self._x_coords.copy(), self._y_upper.copy(), self._y_lower.copy()

    def get_camber_line(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the camber line coordinates.

        Returns:
            Tuple of (x_coords, y_camber)
        """
        return self._x_coords.copy(), self._y_camber.copy()

    def get_closed_profile(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get a closed profile suitable for polygon creation.

        Returns ordered points going around the airfoil clockwise:
        TE (upper) -> LE -> TE (lower)

        The upper surface is reversed (TE to LE), then the lower surface
        continues from LE to TE. The first point of the lower surface is
        skipped because it's the same as the last point of the reversed
        upper surface (both are at the leading edge).

        Returns:
            Tuple of (x_coords, y_coords) for closed profile
        """
        # Upper surface from TE to LE (reversed)
        x_upper = self._x_upper[::-1]
        y_upper = self._y_upper[::-1]

        # Lower surface from LE to TE
        # Skip first point (index 0) which is at LE - same as last point of upper
        x_lower = self._x_lower[1:]
        y_lower = self._y_lower[1:]

        x = np.concatenate([x_upper, x_lower])
        y = np.concatenate([y_upper, y_lower])

        return x, y

    def compute_area(self) -> float:
        """
        Compute the cross-sectional area of the airfoil.

        Returns:
            float: Cross-sectional area
        """
        x, y = self.get_closed_profile()

        # Shoelace formula for polygon area
        n = len(x)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += x[i] * y[j] - x[j] * y[i]

        return abs(area) / 2.0

    def compute_perimeter(self) -> float:
        """
        Compute the perimeter (arc length) of the airfoil.

        Returns:
            float: Perimeter length
        """
        x, y = self.get_closed_profile()

        # Sum of segment lengths
        dx = np.diff(x)
        dy = np.diff(y)
        segments = np.sqrt(dx**2 + dy**2)

        # Add final segment back to start
        final_segment = np.sqrt((x[-1] - x[0]) ** 2 + (y[-1] - y[0]) ** 2)

        return np.sum(segments) + final_segment

    def create_polydata_2d(self) -> "vtk.vtkPolyData":
        """
        Create VTK PolyData for the 2D airfoil profile.

        Returns:
            vtkPolyData: 2D airfoil geometry
        """
        if not HAS_VTK:
            raise ImportError("VTK is required for creating PolyData")

        x, y = self.get_closed_profile()
        n_points = len(x)

        # Create VTK points
        points = vtk.vtkPoints()
        for xi, yi in zip(x, y):
            points.InsertNextPoint(xi, yi, 0.0)

        # Create polygon
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(n_points)
        for i in range(n_points):
            polygon.GetPointIds().SetId(i, i)

        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polygon)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetPolys(cells)

        return polydata

    def create_line_polydata(self) -> "vtk.vtkPolyData":
        """
        Create VTK PolyData as lines (outline) for the airfoil.

        Returns:
            vtkPolyData: Airfoil outline as lines
        """
        if not HAS_VTK:
            raise ImportError("VTK is required for creating PolyData")

        x, y = self.get_closed_profile()
        n_points = len(x)

        points = vtk.vtkPoints()
        for xi, yi in zip(x, y):
            points.InsertNextPoint(xi, yi, 0.0)

        # Create polyline
        lines = vtk.vtkCellArray()
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(n_points + 1)  # +1 to close
        for i in range(n_points):
            polyline.GetPointIds().SetId(i, i)
        polyline.GetPointIds().SetId(n_points, 0)  # Close the loop
        lines.InsertNextCell(polyline)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(lines)

        return polydata

    def create_wing_3d(
        self,
        wing_params: Optional[WingParameters] = None,
        span: float = 10.0,
        taper_ratio: float = 1.0,
        sweep_angle: float = 0.0,
        dihedral_angle: float = 0.0,
        twist_angle: float = 0.0,
        num_sections: int = 20,
    ) -> "vtk.vtkPolyData":
        """
        Create a 3D wing by extruding the airfoil profile.

        Args:
            wing_params: WingParameters object (overrides other args if provided)
            span: Wing semi-span (one side)
            taper_ratio: Tip chord / root chord ratio
            sweep_angle: Leading edge sweep in degrees
            dihedral_angle: Dihedral angle in degrees
            twist_angle: Twist at tip in degrees
            num_sections: Number of spanwise sections

        Returns:
            vtkPolyData: 3D wing surface mesh
        """
        if not HAS_VTK:
            raise ImportError("VTK is required for creating 3D wing")

        if wing_params is not None:
            span = wing_params.span
            taper_ratio = wing_params.taper_ratio
            sweep_angle = wing_params.sweep_angle
            dihedral_angle = wing_params.dihedral_angle
            twist_angle = wing_params.twist_angle
            num_sections = wing_params.num_sections

        # Get base profile
        x_profile, y_profile = self.get_closed_profile()
        n_profile = len(x_profile)

        # Create points for all sections
        all_points = vtk.vtkPoints()

        sweep_rad = math.radians(sweep_angle)
        dihedral_rad = math.radians(dihedral_angle)

        for i in range(num_sections):
            # Spanwise position
            t = i / (num_sections - 1) if num_sections > 1 else 0
            z = t * span

            # Local chord (taper)
            local_chord = self.chord * (1 - t * (1 - taper_ratio))

            # Leading edge position (sweep and dihedral)
            le_x = t * span * math.tan(sweep_rad)
            le_y = t * span * math.tan(dihedral_rad)

            # Local twist
            local_twist = t * twist_angle
            twist_rad = math.radians(local_twist)
            cos_twist = math.cos(twist_rad)
            sin_twist = math.sin(twist_rad)

            # Scale and position profile
            scale = local_chord / self.chord
            for j in range(n_profile):
                # Scale and rotate (twist)
                x_local = x_profile[j] * scale
                y_local = y_profile[j] * scale

                # Apply twist around leading edge
                x_twisted = x_local * cos_twist + y_local * sin_twist
                y_twisted = -x_local * sin_twist + y_local * cos_twist

                # Position in space
                x_final = x_twisted + le_x
                y_final = y_twisted + le_y
                z_final = z

                all_points.InsertNextPoint(x_final, y_final, z_final)

        # Create surface mesh (quads connecting sections)
        polys = vtk.vtkCellArray()

        for i in range(num_sections - 1):
            for j in range(n_profile):
                # Current section indices
                curr_start = i * n_profile
                next_start = (i + 1) * n_profile

                # Quad vertices (wrap around at end)
                j_next = (j + 1) % n_profile

                p0 = curr_start + j
                p1 = curr_start + j_next
                p2 = next_start + j_next
                p3 = next_start + j

                quad = vtk.vtkQuad()
                quad.GetPointIds().SetId(0, p0)
                quad.GetPointIds().SetId(1, p1)
                quad.GetPointIds().SetId(2, p2)
                quad.GetPointIds().SetId(3, p3)
                polys.InsertNextCell(quad)

        # Create polydata
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(all_points)
        polydata.SetPolys(polys)

        # Compute normals for proper rendering
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(polydata)
        normals.ComputePointNormalsOn()
        normals.ComputeCellNormalsOn()
        normals.SplittingOff()
        normals.ConsistencyOn()
        normals.Update()

        return normals.GetOutput()

    def add_pressure_field(
        self,
        polydata: "vtk.vtkPolyData",
        field_name: str = "Pressure",
        mode: str = "x_position",
    ) -> "vtk.vtkPolyData":
        """
        Add a scalar pressure-like field to the polydata.

        Args:
            polydata: VTK PolyData to modify
            field_name: Name for the scalar field
            mode: Computation mode ("x_position", "cp_estimate", "distance")

        Returns:
            Modified polydata with scalar field
        """
        if not HAS_VTK:
            raise ImportError("VTK is required")

        if self.chord <= 0:
            raise ValueError("Chord length must be positive for pressure field computation")

        scalars = vtk.vtkFloatArray()
        scalars.SetName(field_name)
        scalars.SetNumberOfComponents(1)

        n_points = polydata.GetNumberOfPoints()

        for i in range(n_points):
            x, y, z = polydata.GetPoint(i)

            if mode == "x_position":
                # Chordwise position (simulates stagnation at LE)
                value = x / self.chord
            elif mode == "cp_estimate":
                # Simplified pressure coefficient estimate
                # Cp ≈ 1 - (V/V_inf)^2, approximated by position
                x_norm = x / self.chord
                value = 1.0 - 4 * x_norm * (1 - x_norm)
            elif mode == "distance":
                # Distance from chord line
                value = abs(y) / self.chord
            else:
                value = x / self.chord

            scalars.InsertNextValue(value)

        polydata.GetPointData().AddArray(scalars)
        polydata.GetPointData().SetActiveScalars(field_name)

        return polydata

    def export(
        self,
        filename: str,
        format: Optional[str] = None,
        polydata: Optional["vtk.vtkPolyData"] = None,
        binary: bool = True,
    ) -> bool:
        """
        Export the airfoil to a file.

        Args:
            filename: Output filename
            format: Export format (auto-detected from extension if None)
            polydata: PolyData to export (creates 2D if None)
            binary: Use binary format where applicable

        Returns:
            bool: True if export succeeded
        """
        if format is None:
            _, ext = os.path.splitext(filename)
            format = ext.lower().lstrip(".")

        format = format.lower()

        # Create default 2D polydata if not provided
        if polydata is None and format in ["vtp", "stl", "obj", "ply"]:
            if HAS_VTK:
                polydata = self.create_polydata_2d()
            else:
                raise ImportError("VTK is required for mesh export formats")

        if format == "vtp":
            return self._export_vtp(filename, polydata, binary)
        elif format == "stl":
            return self._export_stl(filename, polydata, binary)
        elif format == "obj":
            return self._export_obj(filename, polydata)
        elif format == "ply":
            return self._export_ply(filename, polydata, binary)
        elif format == "csv":
            return self._export_csv(filename)
        elif format == "dat":
            return self._export_dat(filename)
        else:
            raise ValueError(f"Unsupported export format: {format}")

    def _export_vtp(
        self, filename: str, polydata: "vtk.vtkPolyData", binary: bool
    ) -> bool:
        """Export to VTK XML PolyData format."""
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        if binary:
            writer.SetDataModeToBinary()
        else:
            writer.SetDataModeToAscii()
        return writer.Write() == 1

    def _export_stl(
        self, filename: str, polydata: "vtk.vtkPolyData", binary: bool
    ) -> bool:
        """Export to STL format (requires triangulation)."""
        # Triangulate first
        triangulate = vtk.vtkTriangleFilter()
        triangulate.SetInputData(polydata)
        triangulate.Update()

        writer = vtk.vtkSTLWriter()
        writer.SetFileName(filename)
        writer.SetInputData(triangulate.GetOutput())
        if binary:
            writer.SetFileTypeToBinary()
        else:
            writer.SetFileTypeToASCII()
        return writer.Write() == 1

    def _export_obj(self, filename: str, polydata: "vtk.vtkPolyData") -> bool:
        """Export to Wavefront OBJ format."""
        writer = vtk.vtkOBJWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        return writer.Write() == 1

    def _export_ply(
        self, filename: str, polydata: "vtk.vtkPolyData", binary: bool
    ) -> bool:
        """Export to PLY format."""
        writer = vtk.vtkPLYWriter()
        writer.SetFileName(filename)
        writer.SetInputData(polydata)
        if binary:
            writer.SetFileTypeToBinary()
        else:
            writer.SetFileTypeToASCII()
        return writer.Write() == 1

    def _export_csv(self, filename: str) -> bool:
        """Export coordinates to CSV format."""
        x, y_upper, y_lower = self.get_profile_points()

        with open(filename, "w") as f:
            f.write("x,y_upper,y_lower,y_camber,thickness\n")
            for i in range(len(x)):
                f.write(
                    f"{x[i]:.8f},{y_upper[i]:.8f},{y_lower[i]:.8f},"
                    f"{self._y_camber[i]:.8f},{self._thickness[i]:.8f}\n"
                )
        return True

    def _export_dat(self, filename: str) -> bool:
        """Export to XFOIL/XFLR5 DAT format."""
        if self.chord <= 0:
            raise ValueError("Chord length must be positive for DAT export")

        x, y = self.get_closed_profile()

        with open(filename, "w") as f:
            f.write(f"NACA {self.naca_code}\n")
            for xi, yi in zip(x, y):
                f.write(f"  {xi/self.chord:.7f}  {yi/self.chord:.7f}\n")
        return True

    def visualize_2d_matplotlib(
        self,
        show_camber: bool = True,
        show_thickness: bool = False,
        show_points: bool = False,
        show_grid: bool = True,
        title: Optional[str] = None,
        figsize: Tuple[int, int] = (12, 6),
        save_path: Optional[str] = None,
        dpi: int = 150,
    ) -> Optional["plt.Figure"]:
        """
        Visualize the airfoil using matplotlib.

        Args:
            show_camber: Show camber line
            show_thickness: Show thickness distribution
            show_points: Show discrete points
            show_grid: Show grid
            title: Plot title
            figsize: Figure size in inches
            save_path: Path to save figure (if provided)
            dpi: Resolution for saved figure

        Returns:
            matplotlib Figure object or None if saving
        """
        if not HAS_MATPLOTLIB:
            raise ImportError("Matplotlib is required for 2D visualization")

        fig, ax = plt.subplots(figsize=figsize)

        # Get profile data
        x = self._x_coords
        y_upper = self._y_upper
        y_lower = self._y_lower

        # Plot airfoil surface
        ax.fill_between(
            x, y_lower, y_upper, alpha=0.3, color="blue", label="Airfoil"
        )
        ax.plot(x, y_upper, "b-", linewidth=2, label="Upper surface")
        ax.plot(x, y_lower, "b-", linewidth=2, label="Lower surface")

        # Camber line
        if show_camber and not self.params.is_symmetric:
            ax.plot(
                x, self._y_camber, "r--", linewidth=1.5, label="Camber line"
            )

        # Chord line
        ax.axhline(y=0, color="k", linestyle="-", linewidth=0.5, alpha=0.5)

        # Thickness distribution
        if show_thickness:
            ax_twin = ax.twinx()
            ax_twin.plot(
                x, self._thickness, "g:", linewidth=1.5, label="Thickness"
            )
            ax_twin.set_ylabel("Thickness", color="green")
            ax_twin.tick_params(axis="y", labelcolor="green")

        # Show points
        if show_points:
            ax.plot(x, y_upper, "bo", markersize=3)
            ax.plot(x, y_lower, "bo", markersize=3)

        # Formatting
        ax.set_xlabel("x / chord")
        ax.set_ylabel("y / chord")
        ax.set_aspect("equal")
        ax.grid(show_grid, alpha=0.3)
        ax.legend(loc="upper right")

        if title is None:
            title = f"NACA {self.naca_code} Airfoil"
        ax.set_title(title)

        # Add annotations
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
        info_text = (
            f"Chord: {self.chord:.2f}\n"
            f"Thickness: {self.params.thickness*100:.1f}%\n"
            f"Camber: {self.params.max_camber*100:.1f}%\n"
            f"Area: {self.compute_area():.4f}"
        )
        ax.text(
            0.02,
            0.98,
            info_text,
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="top",
            bbox=props,
        )

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
            plt.close()
            return None

        return fig

    def visualize_3d_vtk(
        self,
        polydata: Optional["vtk.vtkPolyData"] = None,
        wing_params: Optional[WingParameters] = None,
        show_edges: bool = True,
        show_axes: bool = True,
        background_color: Tuple[float, float, float] = (0.1, 0.1, 0.15),
        window_size: Tuple[int, int] = (1200, 800),
        window_title: Optional[str] = None,
        save_screenshot: Optional[str] = None,
    ):
        """
        Visualize the airfoil or wing using VTK.

        Args:
            polydata: PolyData to visualize (creates 3D wing if None)
            wing_params: Wing parameters for 3D generation
            show_edges: Show mesh edges
            show_axes: Show orientation axes
            background_color: Background color (RGB)
            window_size: Window size in pixels
            window_title: Window title
            save_screenshot: Path to save screenshot
        """
        if not HAS_VTK:
            raise ImportError("VTK is required for 3D visualization")

        # Create geometry if not provided
        if polydata is None:
            if wing_params is not None:
                polydata = self.create_wing_3d(wing_params)
            else:
                polydata = self.create_wing_3d(span=5.0)

        # Add pressure field
        polydata = self.add_pressure_field(polydata, "Pressure", "cp_estimate")

        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.SetScalarModeToUsePointData()
        scalar_range = polydata.GetPointData().GetScalars().GetRange()
        mapper.SetScalarRange(scalar_range)

        # Color lookup table (blue-white-red for pressure)
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(256)
        lut.SetHueRange(0.667, 0.0)  # Blue to red
        lut.Build()
        mapper.SetLookupTable(lut)

        # Create actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        if show_edges:
            actor.GetProperty().SetEdgeVisibility(1)
            actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
            actor.GetProperty().SetLineWidth(0.5)

        actor.GetProperty().SetInterpolationToPhong()

        # Create renderer
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(*background_color)
        renderer.AddActor(actor)

        # Add scalar bar (color legend)
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(lut)
        scalar_bar.SetTitle("Cp")
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetWidth(0.08)
        scalar_bar.SetHeight(0.4)
        scalar_bar.GetPositionCoordinate().SetValue(0.9, 0.3)
        renderer.AddActor2D(scalar_bar)

        # Add title text
        if window_title is None:
            window_title = f"NACA {self.naca_code} Wing"

        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(window_title)
        text_actor.GetTextProperty().SetFontSize(24)
        text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)
        text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        text_actor.GetPositionCoordinate().SetValue(0.02, 0.95)
        renderer.AddActor2D(text_actor)

        # Create render window
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(*window_size)
        render_window.SetWindowName(window_title)

        # Create interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)

        # Add axes widget
        if show_axes:
            axes = vtk.vtkAxesActor()
            axes_widget = vtk.vtkOrientationMarkerWidget()
            axes_widget.SetOrientationMarker(axes)
            axes_widget.SetInteractor(interactor)
            axes_widget.SetViewport(0.0, 0.0, 0.15, 0.2)
            axes_widget.SetEnabled(1)
            axes_widget.InteractiveOff()

        # Set camera for good initial view
        camera = renderer.GetActiveCamera()
        camera.SetPosition(1.5, 0.5, 1.0)
        camera.SetFocalPoint(0.5, 0, 0.5)
        camera.SetViewUp(0, 1, 0)
        renderer.ResetCamera()

        # Initialize and render
        interactor.Initialize()
        render_window.Render()

        # Save screenshot if requested
        if save_screenshot:
            window_to_image = vtk.vtkWindowToImageFilter()
            window_to_image.SetInput(render_window)
            window_to_image.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetFileName(save_screenshot)
            writer.SetInputConnection(window_to_image.GetOutputPort())
            writer.Write()

        # Start interaction
        interactor.Start()

    def __repr__(self) -> str:
        return (
            f"NACAGenerator(naca_code='{self.naca_code}', "
            f"chord={self.chord}, num_points={self.num_points})"
        )

    def __str__(self) -> str:
        return (
            f"NACA {self.naca_code} Airfoil\n"
            f"  Type: {self.params.airfoil_type.value}\n"
            f"  Chord: {self.chord}\n"
            f"  Thickness: {self.params.thickness*100:.1f}%\n"
            f"  Max Camber: {self.params.max_camber*100:.1f}%\n"
            f"  Camber Pos: {self.params.camber_position*100:.1f}% chord\n"
            f"  Symmetric: {self.params.is_symmetric}\n"
            f"  Area: {self.compute_area():.6f}\n"
            f"  Perimeter: {self.compute_perimeter():.6f}"
        )


def create_comparison_visualization(
    airfoil_codes: List[str], chord: float = 1.0, num_points: int = 100
) -> Optional["plt.Figure"]:
    """
    Create a comparison plot of multiple NACA airfoils.

    Args:
        airfoil_codes: List of NACA codes to compare
        chord: Common chord length
        num_points: Number of points per airfoil

    Returns:
        matplotlib Figure object
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("Matplotlib is required for comparison visualization")

    if not airfoil_codes:
        raise ValueError("At least one airfoil code is required for comparison")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    n_airfoils = len(airfoil_codes)
    colors = plt.cm.tab10(np.linspace(0, 1, max(n_airfoils, 1)))

    # Profile comparison (top-left)
    ax1 = axes[0, 0]
    for i, code in enumerate(airfoil_codes):
        gen = NACAGenerator(code, chord=chord, num_points=num_points)
        x, y = gen.get_closed_profile()
        ax1.plot(x, y, color=colors[i], linewidth=2, label=f"NACA {code}")

    ax1.set_xlabel("x / chord")
    ax1.set_ylabel("y / chord")
    ax1.set_title("Airfoil Profile Comparison")
    ax1.set_aspect("equal")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Thickness comparison (top-right)
    ax2 = axes[0, 1]
    for i, code in enumerate(airfoil_codes):
        gen = NACAGenerator(code, chord=chord, num_points=num_points)
        ax2.plot(
            gen._x_coords,
            gen._thickness,
            color=colors[i],
            linewidth=2,
            label=f"NACA {code}",
        )

    ax2.set_xlabel("x / chord")
    ax2.set_ylabel("Thickness")
    ax2.set_title("Thickness Distribution")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Camber comparison (bottom-left)
    ax3 = axes[1, 0]
    for i, code in enumerate(airfoil_codes):
        gen = NACAGenerator(code, chord=chord, num_points=num_points)
        ax3.plot(
            gen._x_coords,
            gen._y_camber,
            color=colors[i],
            linewidth=2,
            label=f"NACA {code}",
        )

    ax3.set_xlabel("x / chord")
    ax3.set_ylabel("Camber")
    ax3.set_title("Camber Line Comparison")
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Properties table (bottom-right)
    ax4 = axes[1, 1]
    ax4.axis("off")

    table_data = []
    for code in airfoil_codes:
        gen = NACAGenerator(code, chord=chord, num_points=num_points)
        table_data.append(
            [
                f"NACA {code}",
                f"{gen.params.thickness*100:.1f}%",
                f"{gen.params.max_camber*100:.1f}%",
                f"{gen.compute_area():.4f}",
                f"{gen.compute_perimeter():.4f}",
            ]
        )

    table = ax4.table(
        cellText=table_data,
        colLabels=["Airfoil", "Thickness", "Camber", "Area", "Perimeter"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    ax4.set_title("Airfoil Properties", y=0.8)

    plt.tight_layout()
    return fig


def main():
    """
    Main function demonstrating the NACA airfoil generator capabilities.
    """
    parser = argparse.ArgumentParser(
        description="NACA Airfoil Generator - Generate and visualize NACA airfoils",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 0012                     # Generate NACA 0012 (symmetric)
  %(prog)s 2412                     # Generate NACA 2412 (cambered)
  %(prog)s 23012                    # Generate NACA 23012 (5-digit)
  %(prog)s 2412 --export wing.stl   # Export to STL
  %(prog)s 0015 --wing --span 10    # Create 3D wing
  %(prog)s 2412 --compare 0012 4415 # Compare multiple airfoils
        """,
    )

    parser.add_argument(
        "naca_code",
        nargs="?",
        default="2412",
        help="NACA airfoil designation (e.g., 0012, 2412, 23012)",
    )
    parser.add_argument(
        "--chord",
        "-c",
        type=float,
        default=1.0,
        help="Chord length (default: 1.0)",
    )
    parser.add_argument(
        "--points",
        "-n",
        type=int,
        default=100,
        help="Number of points per surface (default: 100)",
    )
    parser.add_argument(
        "--export",
        "-e",
        type=str,
        help="Export to file (vtp, stl, obj, csv, dat)",
    )
    parser.add_argument(
        "--visualize",
        "-v",
        action="store_true",
        help="Show 2D matplotlib visualization",
    )
    parser.add_argument(
        "--vtk", action="store_true", help="Show 3D VTK visualization"
    )
    parser.add_argument(
        "--wing",
        "-w",
        action="store_true",
        help="Generate 3D wing instead of 2D profile",
    )
    parser.add_argument(
        "--span",
        type=float,
        default=5.0,
        help="Wing span for 3D generation (default: 5.0)",
    )
    parser.add_argument(
        "--taper",
        type=float,
        default=1.0,
        help="Wing taper ratio (default: 1.0)",
    )
    parser.add_argument(
        "--sweep",
        type=float,
        default=0.0,
        help="Wing sweep angle in degrees (default: 0.0)",
    )
    parser.add_argument(
        "--compare", nargs="+", help="Compare with other NACA codes"
    )
    parser.add_argument(
        "--info", "-i", action="store_true", help="Print airfoil information"
    )
    parser.add_argument(
        "--screenshot", type=str, help="Save visualization to image file"
    )

    args = parser.parse_args()

    # Create generator
    try:
        generator = NACAGenerator(
            args.naca_code, chord=args.chord, num_points=args.points
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Print info
    if args.info or not any(
        [args.export, args.visualize, args.vtk, args.compare]
    ):
        print(generator)
        print()

    # Export
    if args.export:
        success = generator.export(args.export)
        if success:
            print(f"Exported to: {args.export}")
        else:
            print(f"Failed to export to: {args.export}", file=sys.stderr)
            return 1

    # Comparison visualization
    if args.compare:
        all_codes = [args.naca_code] + args.compare
        fig = create_comparison_visualization(all_codes, args.chord, args.points)
        if args.screenshot:
            plt.savefig(args.screenshot, dpi=150, bbox_inches="tight")
            print(f"Saved comparison to: {args.screenshot}")
        else:
            plt.show()

    # 2D visualization
    elif args.visualize:
        fig = generator.visualize_2d_matplotlib()
        if args.screenshot:
            plt.savefig(args.screenshot, dpi=150, bbox_inches="tight")
            print(f"Saved 2D plot to: {args.screenshot}")
        else:
            plt.show()

    # 3D VTK visualization
    elif args.vtk or args.wing:
        wing_params = WingParameters(
            span=args.span, taper_ratio=args.taper, sweep_angle=args.sweep
        )
        generator.visualize_3d_vtk(
            wing_params=wing_params, save_screenshot=args.screenshot
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
