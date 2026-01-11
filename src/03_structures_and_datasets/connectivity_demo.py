"""
VTK Connectivity Demo: Understanding How Points Connect to Form Cells

This module provides an interactive demonstration of the fundamental concept
of CONNECTIVITY in VTK - how points are connected to form different structures.

What is Connectivity?
---------------------
Connectivity is the relationship between points that defines the mesh structure.
The SAME SET OF POINTS can be interpreted completely differently based on how
they are connected:

- No connectivity: Just a point cloud (individual vertices)
- Linear connectivity: Points form lines/edges
- Surface connectivity: Points form triangles/polygons
- Volume connectivity: Points form 3D cells (tetrahedra, hexahedra)

Key Concept:
-----------
Points alone have NO structure. Connectivity DEFINES the structure.

    *    *    *    *         <- Points without connectivity (point cloud)

    *----*----*----*         <- Points with linear connectivity (polyline)

    *----*                   <- Points with surface connectivity (polygon)
    |\\  |
    | \\ |
    |  \\|
    *----*

This is why connectivity is fundamental to VTK:
- It determines how your data is interpreted
- It affects rendering (points vs lines vs surfaces)
- It impacts computations (normals, interpolation, etc.)

Real-World Relevance:
--------------------
- CFD: Mesh cells are defined by point connectivity
- FEA: Element types are determined by connectivity patterns
- Medical Imaging: Surface extraction creates new connectivity
- Point Cloud Processing: Creating meshes requires defining connectivity

Usage:
------
Run this script directly to launch the interactive viewer:
    python connectivity_demo.py

Select different connectivity modes from the combo box to see how
the SAME POINTS can form completely different structures!
"""

import sys

import vtk
from PyQt6.QtCore import Qt
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

# Visual constants
CAMERA_AZIMUTH = 30
CAMERA_ELEVATION = 20
VERTEX_SPHERE_RADIUS = 0.1


def create_grid_points():
    """
    Create a 3x3 grid of points that will be connected in different ways.

    The points are arranged in a regular grid pattern:

        6----7----8
        |    |    |
        3----4----5
        |    |    |
        0----1----2

    These 9 points can be connected in various ways to demonstrate
    how connectivity defines structure.

    Returns:
        vtkPoints object containing 9 points in a 3x3 grid
    """
    points = vtk.vtkPoints()

    # Create a 3x3 grid of points
    for j in range(3):
        for i in range(3):
            x = i - 1  # Center the grid: x from -1 to 1
            y = j - 1  # Center the grid: y from -1 to 1
            z = 0
            points.InsertNextPoint(x, y, z)

    return points


def create_point_cloud(points):
    """
    Create a dataset with NO connectivity - just individual vertices.

    This represents a point cloud where each point exists independently.
    There is no relationship between points - they are just locations in space.

    In VTK, this is represented by vtkVertex cells, where each cell
    contains exactly one point.

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Add each point as a separate vertex cell
    for i in range(points.GetNumberOfPoints()):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())

    description = (
        "NO CONNECTIVITY: Each point is an isolated vertex.\n"
        "Points exist independently with no relationships between them.\n"
        "This represents a point cloud - just positions in space."
    )

    return ugrid, description


def create_lines(points):
    """
    Create a dataset with LINEAR connectivity - points form lines.

    Points are connected sequentially to form a continuous path.
    This demonstrates 1D connectivity where each interior point
    connects to exactly 2 neighbors (before and after).

    The path traces: 0 → 1 → 2 → 5 → 8 → 7 → 6 → 3 → 4
    (a snake-like pattern through the grid)

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Create a polyline that snakes through the grid
    # Pattern: bottom row left-to-right, up, top row right-to-left,
    # down to center
    point_order = [0, 1, 2, 5, 8, 7, 6, 3, 4]

    poly_line = vtk.vtkPolyLine()
    poly_line.GetPointIds().SetNumberOfIds(len(point_order))
    for i, pt_id in enumerate(point_order):
        poly_line.GetPointIds().SetId(i, pt_id)

    ugrid.InsertNextCell(poly_line.GetCellType(), poly_line.GetPointIds())

    description = (
        "LINEAR CONNECTIVITY: Points form a continuous path (polyline).\n"
        "Each point connects to its neighbors in the path.\n"
        "The same points now represent a 1D structure (curve/path)."
    )

    return ugrid, description


def create_triangles(points):
    """
    Create a dataset with SURFACE connectivity - points form triangles.

    The 9 points are connected to form 8 triangles that tile the grid.
    This demonstrates 2D surface connectivity where points are shared
    between adjacent triangles.

        6----7----8
        |\\ 5|\\7 |
        |4\\ | 6\\|
        3----4----5
        |\\1 |\\3 |
        |0\\ | 2\\|
        0----1----2

    Each triangle shares edges and vertices with its neighbors.

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Define 8 triangles that cover the 3x3 grid
    # Each row of 2 grid cells has 4 triangles
    triangles = [
        # Bottom row
        (0, 1, 3),  # Triangle 0
        (1, 4, 3),  # Triangle 1
        (1, 2, 4),  # Triangle 2
        (2, 5, 4),  # Triangle 3
        # Top row
        (3, 4, 6),  # Triangle 4
        (4, 7, 6),  # Triangle 5
        (4, 5, 7),  # Triangle 6
        (5, 8, 7),  # Triangle 7
    ]

    for tri_points in triangles:
        triangle = vtk.vtkTriangle()
        for i, pt_id in enumerate(tri_points):
            triangle.GetPointIds().SetId(i, pt_id)
        ugrid.InsertNextCell(triangle.GetCellType(), triangle.GetPointIds())

    description = (
        "SURFACE CONNECTIVITY: Points form triangles (2D cells).\n"
        "Points are shared between adjacent triangles.\n"
        "The same points now form a continuous surface mesh."
    )

    return ugrid, description


def create_quads(points):
    """
    Create a dataset with QUAD connectivity - points form quadrilaterals.

    The 9 points form 4 quadrilateral cells that tile the grid.

        6----7----8
        | Q2 | Q3 |
        3----4----5
        | Q0 | Q1 |
        0----1----2

    Quads share edges and vertices with their neighbors, similar to
    triangles but with 4 points per cell instead of 3.

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Define 4 quads that cover the 3x3 grid
    quads = [
        (0, 1, 4, 3),  # Quad 0: bottom-left
        (1, 2, 5, 4),  # Quad 1: bottom-right
        (3, 4, 7, 6),  # Quad 2: top-left
        (4, 5, 8, 7),  # Quad 3: top-right
    ]

    for quad_points in quads:
        quad = vtk.vtkQuad()
        for i, pt_id in enumerate(quad_points):
            quad.GetPointIds().SetId(i, pt_id)
        ugrid.InsertNextCell(quad.GetCellType(), quad.GetPointIds())

    description = (
        "QUAD CONNECTIVITY: Points form quadrilaterals (2D cells).\n"
        "4 quads share the center point (point 4).\n"
        "Same points, different cell topology than triangles."
    )

    return ugrid, description


def create_polygon(points):
    """
    Create a dataset where all boundary points form a single polygon.

    The outer 8 points form a single closed polygon, while the center
    point is included as a separate vertex.

        6----7----8
        |         |
        3    4    5
        |         |
        0----1----2

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Outer boundary points in counter-clockwise order
    boundary_order = [0, 1, 2, 5, 8, 7, 6, 3]

    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(len(boundary_order))
    for i, pt_id in enumerate(boundary_order):
        polygon.GetPointIds().SetId(i, pt_id)

    ugrid.InsertNextCell(polygon.GetCellType(), polygon.GetPointIds())

    # Add center point as isolated vertex
    vertex = vtk.vtkVertex()
    vertex.GetPointIds().SetId(0, 4)
    ugrid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())

    description = (
        "POLYGON CONNECTIVITY: Boundary points form one closed polygon.\n"
        "8 points form the perimeter, center point is isolated.\n"
        "Polygons can have any number of sides."
    )

    return ugrid, description


def create_3d_cells(points):
    """
    Create a dataset with 3D connectivity by extruding points into 3D.

    This modifies the z-coordinates to create a 3D structure with
    two layers of points, then connects them as hexahedral cells.

    Args:
        points: vtkPoints object (will be modified to 3D)

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    # Create new points with 3D structure (two layers)
    points_3d = vtk.vtkPoints()

    # Bottom layer (z = -0.5)
    for j in range(3):
        for i in range(3):
            x = i - 1
            y = j - 1
            z = -0.5
            points_3d.InsertNextPoint(x, y, z)

    # Top layer (z = 0.5)
    for j in range(3):
        for i in range(3):
            x = i - 1
            y = j - 1
            z = 0.5
            points_3d.InsertNextPoint(x, y, z)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points_3d)

    # Create 4 hexahedral cells
    # Bottom layer indices: 0-8, Top layer indices: 9-17
    hexahedra = [
        # (bottom-front-left, bottom-front-right, bottom-back-right,
        #  bottom-back-left, top-front-left, top-front-right, top-back-right,
        #  top-back-left)
        (0, 1, 4, 3, 9, 10, 13, 12),   # Hex 0
        (1, 2, 5, 4, 10, 11, 14, 13),  # Hex 1
        (3, 4, 7, 6, 12, 13, 16, 15),  # Hex 2
        (4, 5, 8, 7, 13, 14, 17, 16),  # Hex 3
    ]

    for hex_points in hexahedra:
        hexa = vtk.vtkHexahedron()
        for i, pt_id in enumerate(hex_points):
            hexa.GetPointIds().SetId(i, pt_id)
        ugrid.InsertNextCell(hexa.GetCellType(), hexa.GetPointIds())

    description = (
        "3D CONNECTIVITY: Points form hexahedral (3D) cells.\n"
        "18 points in two layers, connected as 4 hexahedra.\n"
        "Volumetric cells for CFD/FEA simulations."
    )

    return ugrid, description


# Connectivity modes and their creation functions
CONNECTIVITY_MODES = {
    "No Connectivity (Point Cloud)": create_point_cloud,
    "Linear (Polyline)": create_lines,
    "Surface (Triangles)": create_triangles,
    "Surface (Quads)": create_quads,
    "Polygon (Boundary)": create_polygon,
    "Volume (Hexahedra)": create_3d_cells,
}


class ConnectivityDemo(QMainWindow):
    """
    Interactive VTK Connectivity Demo using PyQt6.

    This widget demonstrates how the same set of points can form
    completely different structures based on their connectivity.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("VTK Connectivity Demo: How Points Form Structures")
        self.resize(1000, 750)

        # Create the base points
        self.base_points = create_grid_points()

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create header with combo box
        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Connectivity Mode:"))

        self.combo_box = QComboBox()
        self.combo_box.addItems(CONNECTIVITY_MODES.keys())
        self.combo_box.setMinimumWidth(250)
        self.combo_box.currentTextChanged.connect(self.on_mode_changed)
        header_layout.addWidget(self.combo_box)

        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Create info label
        self.info_label = QLabel()
        self.info_label.setStyleSheet(
            "QLabel { background-color: #2b2b2b; color: #ffffff; padding: 10px; "
            "border-radius: 5px; font-size: 12px; }"
        )
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.info_label.setWordWrap(True)
        self.info_label.setMinimumHeight(100)
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
        self.edge_actor = vtk.vtkActor()
        self.point_actor = vtk.vtkActor()
        self.renderer.AddActor(self.surface_actor)
        self.renderer.AddActor(self.edge_actor)
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

        # Display first mode
        self.on_mode_changed(self.combo_box.currentText())

    def on_mode_changed(self, mode):
        """Handle combo box selection change."""
        self.update_visualization(mode)

    def update_visualization(self, mode):
        """Update the visualization based on the selected mode."""
        if mode not in CONNECTIVITY_MODES:
            return

        # Create new points for each mode (some modes modify points)
        points = create_grid_points()

        # Create the dataset with the selected connectivity
        create_func = CONNECTIVITY_MODES[mode]
        ugrid, description = create_func(points)

        # Update info label
        num_points = ugrid.GetNumberOfPoints()
        num_cells = ugrid.GetNumberOfCells()

        # Get cell types present
        cell_types = set()
        for i in range(num_cells):
            cell_types.add(ugrid.GetCellTypeName(ugrid.GetCellType(i)))
        cell_type_str = ", ".join(sorted(cell_types))

        info_text = (
            f"<b>Mode: {mode}</b><br/>"
            f"{description}<br/><br/>"
            f"<b>Points:</b> {num_points} | <b>Cells:</b> {num_cells} | "
            f"<b>Cell Types:</b> {cell_type_str}"
        )
        self.info_label.setText(info_text)

        # Create surface mapper with edge visibility
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(ugrid)

        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)
        self.surface_actor.GetProperty().SetOpacity(0.7)
        self.surface_actor.GetProperty().SetEdgeVisibility(1)
        self.surface_actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.4)
        self.surface_actor.GetProperty().SetLineWidth(3)

        # Create explicit edge extraction for better visibility
        edges = vtk.vtkExtractEdges()
        edges.SetInputData(ugrid)

        edge_mapper = vtk.vtkPolyDataMapper()
        edge_mapper.SetInputConnection(edges.GetOutputPort())

        self.edge_actor.SetMapper(edge_mapper)
        self.edge_actor.GetProperty().SetColor(0.1, 0.1, 0.4)
        self.edge_actor.GetProperty().SetLineWidth(4)

        # Create point visualization using spheres
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

        # Reset camera
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


def print_educational_summary():
    """Print educational information about connectivity in VTK."""
    print("\n" + "=" * 70)
    print("VTK CONNECTIVITY: Understanding How Points Form Structures")
    print("=" * 70)
    print("""
┌─────────────────────────────────────────────────────────────────────┐
│                    CONNECTIVITY IN VTK                              │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  What is Connectivity?                                              │
│  ─────────────────────                                              │
│  Connectivity defines HOW POINTS ARE CONNECTED to form cells.       │
│  The SAME POINTS can form completely different structures           │
│  depending on their connectivity:                                   │
│                                                                     │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  CONNECTIVITY TYPES:                                                │
│                                                                     │
│  • No Connectivity (0D): Points exist independently (point cloud)   │
│                                                                     │
│    *    *    *                                                      │
│    *    *    *                                                      │
│    *    *    *                                                      │
│                                                                     │
│  • Linear Connectivity (1D): Points form lines/paths                │
│                                                                     │
│    *────*────*                                                      │
│              │                                                      │
│    *    *    *                                                      │
│    │                                                                │
│    *────*────*                                                      │
│                                                                     │
│  • Surface Connectivity (2D): Points form triangles/quads           │
│                                                                     │
│    *────*────*                                                      │
│    │╲   │╲   │                                                      │
│    │ ╲  │ ╲  │                                                      │
│    *────*────*                                                      │
│    │╲   │╲   │                                                      │
│    │ ╲  │ ╲  │                                                      │
│    *────*────*                                                      │
│                                                                     │
│  • Volume Connectivity (3D): Points form hexahedra/tetrahedra       │
│                                                                     │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  WHY THIS MATTERS:                                                  │
│                                                                     │
│  • CFD/FEA: Connectivity defines the mesh structure                 │
│  • Rendering: Lines vs surfaces vs volumes                          │
│  • Interpolation: How values are interpolated between points        │
│  • Computation: Normals, gradients, and other derived quantities    │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

In this demonstration:
- 9 points arranged in a 3x3 grid
- Select different connectivity modes to see how they form structures
- Notice: SAME POINTS, DIFFERENT STRUCTURES based on connectivity!
""")


def main():
    """Launch the VTK Connectivity Demo application."""
    print_educational_summary()

    print("\nStarting interactive demonstration...")
    print("Use the combo box to select different connectivity modes.")
    print("Watch how the SAME POINTS form different structures!")
    print("-" * 70)

    app = QApplication(sys.argv)
    window = ConnectivityDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
