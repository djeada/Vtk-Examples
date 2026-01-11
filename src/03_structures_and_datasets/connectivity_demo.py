"""
VTK Connectivity Demo: Understanding the Multiple Levels of Connectivity

This module provides an interactive demonstration of CONNECTIVITY in VTK,
covering multiple levels of connectivity that define mesh structures.

Levels of Connectivity:
-----------------------

1️⃣ POINT-CELL CONNECTIVITY (Basic Level)
   "Which points belong to this cell?"
   - A triangle uses 3 points
   - A quad uses 4 points
   - A hexahedron uses 8 points
   This is the foundation - defining cells from points.

2️⃣ CELL-CELL CONNECTIVITY (Adjacency)
   "Which cells touch this cell?"
   Cells may share:
   - A face (strong adjacency)
   - An edge (medium adjacency)
   - A vertex (weak adjacency)
   Important for: neighborhood queries, gradient calculations, PDE solvers

3️⃣ LOGICAL CONNECTIVITY (Index Structure)
   "Is there a predictable pattern to neighbors?"
   - Structured grids: neighbors follow (i±1, j±1, k±1) pattern
   - Unstructured grids: arbitrary connectivity, requires lookup
   This is why structured grids are faster - connectivity is implicit!

4️⃣ TOPOLOGICAL CONSTRAINTS
   "What rules govern the mesh structure?"
   - StructuredGrid: fixed topology, cannot add/remove cells
   - UnstructuredGrid: arbitrary topology, flexible but slower
   - Mixed cell types: allowed in unstructured, not in structured

Key Insight:
-----------
Points alone have NO structure. Connectivity at ALL LEVELS defines:
- How data is interpreted (points vs surfaces vs volumes)
- How cells relate to each other (neighbors, boundaries)
- How efficiently queries can be performed (implicit vs explicit)

Usage:
------
Run this script directly to launch the interactive viewer:
    python connectivity_demo.py

The demo shows:
- Point-Cell connectivity: How points form different cell types
- Cell-Cell adjacency: How cells share faces/edges/vertices
- The difference between structured and unstructured connectivity
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
    # VTK hexahedron node ordering (8 nodes):
    #   - Nodes 0-3: Bottom face, counter-clockwise when viewed from above
    #   - Nodes 4-7: Top face, counter-clockwise when viewed from above
    #   - Node 4 is directly above node 0, 5 above 1, etc.
    #
    # Grid layout:
    #   Bottom layer (z=-0.5), indices 0-8:     Top layer (z=0.5), indices 9-17:
    #     6----7----8                             15---16---17
    #     |    |    |                             |    |    |
    #     3----4----5                             12---13---14
    #     |    |    |                             |    |    |
    #     0----1----2                             9----10---11
    #
    hexahedra = [
        # (front-left-bot, front-right-bot, back-right-bot, back-left-bot,
        #  front-left-top, front-right-top, back-right-top, back-left-top)
        (0, 1, 4, 3, 9, 10, 13, 12),   # Hex 0: bottom-left in XY
        (1, 2, 5, 4, 10, 11, 14, 13),  # Hex 1: bottom-right in XY
        (3, 4, 7, 6, 12, 13, 16, 15),  # Hex 2: top-left in XY
        (4, 5, 8, 7, 13, 14, 17, 16),  # Hex 3: top-right in XY
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


def create_cell_adjacency_demo(points):
    """
    Create a dataset demonstrating CELL-CELL CONNECTIVITY (adjacency).

    This shows how cells relate to each other through shared geometry:
    - Face adjacency: Cells sharing a complete face
    - Edge adjacency: Cells sharing an edge
    - Vertex adjacency: Cells sharing only a vertex

    The 4 quads share the center point, demonstrating vertex adjacency.
    Adjacent quads share edges, demonstrating edge adjacency.

        6----7----8
        | Q2 | Q3 |     Q0-Q1: share edge (1,4)
        3----4----5     Q0-Q2: share edge (3,4)
        | Q0 | Q1 |     All 4: share vertex (4)
        0----1----2

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

    # Add cell data to color cells and show adjacency info
    cell_colors = vtk.vtkIntArray()
    cell_colors.SetName("CellID")
    cell_colors.SetNumberOfComponents(1)
    for i in range(4):
        cell_colors.InsertNextValue(i)
    ugrid.GetCellData().SetScalars(cell_colors)

    description = (
        "CELL-CELL CONNECTIVITY (Adjacency): How cells relate to neighbors.\n"
        "• Q0↔Q1: Share edge (points 1-4) - EDGE ADJACENCY\n"
        "• Q0↔Q2: Share edge (points 3-4) - EDGE ADJACENCY\n"
        "• All 4 quads share vertex (point 4) - VERTEX ADJACENCY\n"
        "Cell neighbors matter for gradient calculations & PDE solvers."
    )

    return ugrid, description


def create_structured_connectivity_demo(points):
    """
    Create a STRUCTURED GRID demonstrating implicit/logical connectivity.

    In structured grids, connectivity is IMPLICIT based on indices:
    - Cell (i,j) has neighbors at (i±1, j) and (i, j±1)
    - No need to store explicit neighbor lists
    - This makes structured grids faster for neighbor lookups

    The grid has logical structure where neighbors follow index patterns:
        Cell[i,j] neighbors: Cell[i-1,j], Cell[i+1,j], Cell[i,j-1], Cell[i,j+1]

    Args:
        points: vtkPoints object (not used - creates new structured grid)

    Returns:
        tuple: (vtkStructuredGrid, description string)
    """
    # Create a 3x3 structured grid
    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(3, 3, 1)

    # Create points for the structured grid
    pts = vtk.vtkPoints()
    for j in range(3):
        for i in range(3):
            x = i - 1
            y = j - 1
            z = 0
            pts.InsertNextPoint(x, y, z)

    sgrid.SetPoints(pts)

    # Add cell data showing logical indices
    cell_indices = vtk.vtkIntArray()
    cell_indices.SetName("LogicalIndex")
    cell_indices.SetNumberOfComponents(1)
    # 2x2 cells in a 3x3 point grid
    for j in range(2):
        for i in range(2):
            cell_indices.InsertNextValue(j * 2 + i)
    sgrid.GetCellData().SetScalars(cell_indices)

    description = (
        "STRUCTURED CONNECTIVITY: Implicit neighbor pattern.\n"
        "• Cell(i,j) neighbors are at (i±1,j) and (i,j±1)\n"
        "• NO explicit neighbor storage needed - derived from indices\n"
        "• FAST lookups: O(1) neighbor access vs O(n) for unstructured\n"
        "This is why structured grids are preferred when geometry allows."
    )

    return sgrid, description


def create_unstructured_connectivity_demo(points):
    """
    Create an UNSTRUCTURED GRID demonstrating explicit connectivity.

    In unstructured grids, connectivity must be EXPLICITLY stored:
    - No predictable pattern for neighbors
    - Arbitrary cell arrangements and types
    - More flexible but slower neighbor lookups

    This example uses mixed cell types to show unstructured flexibility:
    - Triangles and quads mixed together
    - No regular index pattern for neighbors

    Args:
        points: vtkPoints object with the point coordinates

    Returns:
        tuple: (vtkUnstructuredGrid, description string)
    """
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    # Mix of triangles and quads - this is only possible in unstructured grids!
    # Use triangles on left, quad on right to show mixed topology

    # Left side: 2 triangles
    tri1 = vtk.vtkTriangle()
    tri1.GetPointIds().SetId(0, 0)
    tri1.GetPointIds().SetId(1, 1)
    tri1.GetPointIds().SetId(2, 3)
    ugrid.InsertNextCell(tri1.GetCellType(), tri1.GetPointIds())

    tri2 = vtk.vtkTriangle()
    tri2.GetPointIds().SetId(0, 1)
    tri2.GetPointIds().SetId(1, 4)
    tri2.GetPointIds().SetId(2, 3)
    ugrid.InsertNextCell(tri2.GetCellType(), tri2.GetPointIds())

    tri3 = vtk.vtkTriangle()
    tri3.GetPointIds().SetId(0, 3)
    tri3.GetPointIds().SetId(1, 4)
    tri3.GetPointIds().SetId(2, 6)
    ugrid.InsertNextCell(tri3.GetCellType(), tri3.GetPointIds())

    tri4 = vtk.vtkTriangle()
    tri4.GetPointIds().SetId(0, 4)
    tri4.GetPointIds().SetId(1, 7)
    tri4.GetPointIds().SetId(2, 6)
    ugrid.InsertNextCell(tri4.GetCellType(), tri4.GetPointIds())

    # Right side: 2 quads
    quad1 = vtk.vtkQuad()
    quad1.GetPointIds().SetId(0, 1)
    quad1.GetPointIds().SetId(1, 2)
    quad1.GetPointIds().SetId(2, 5)
    quad1.GetPointIds().SetId(3, 4)
    ugrid.InsertNextCell(quad1.GetCellType(), quad1.GetPointIds())

    quad2 = vtk.vtkQuad()
    quad2.GetPointIds().SetId(0, 4)
    quad2.GetPointIds().SetId(1, 5)
    quad2.GetPointIds().SetId(2, 8)
    quad2.GetPointIds().SetId(3, 7)
    ugrid.InsertNextCell(quad2.GetCellType(), quad2.GetPointIds())

    # Add cell type indicator
    cell_types = vtk.vtkIntArray()
    cell_types.SetName("CellType")
    cell_types.SetNumberOfComponents(1)
    cell_types.InsertNextValue(0)  # Triangle
    cell_types.InsertNextValue(0)  # Triangle
    cell_types.InsertNextValue(0)  # Triangle
    cell_types.InsertNextValue(0)  # Triangle
    cell_types.InsertNextValue(1)  # Quad
    cell_types.InsertNextValue(1)  # Quad
    ugrid.GetCellData().SetScalars(cell_types)

    description = (
        "UNSTRUCTURED CONNECTIVITY: Explicit, arbitrary topology.\n"
        "• MIXED CELL TYPES: 4 triangles + 2 quads in same mesh\n"
        "• No predictable neighbor pattern - must store explicitly\n"
        "• FLEXIBLE: Can represent any geometry, any topology\n"
        "Trade-off: Flexibility vs. speed of structured grids."
    )

    return ugrid, description


# Connectivity modes and their creation functions
CONNECTIVITY_MODES = {
    "1. Point-Cell: No Connectivity (Point Cloud)": create_point_cloud,
    "2. Point-Cell: Linear (Polyline)": create_lines,
    "3. Point-Cell: Surface (Triangles)": create_triangles,
    "4. Point-Cell: Surface (Quads)": create_quads,
    "5. Point-Cell: Volume (Hexahedra)": create_3d_cells,
    "6. Cell-Cell: Adjacency Demo": create_cell_adjacency_demo,
    "7. Logical: Structured Grid (Implicit)": create_structured_connectivity_demo,
    "8. Logical: Unstructured Grid (Mixed Types)": create_unstructured_connectivity_demo,
}


class ConnectivityDemo(QMainWindow):
    """
    Interactive VTK Connectivity Demo using PyQt6.

    This widget demonstrates multiple levels of connectivity in VTK:
    - Point-Cell: How points form cells
    - Cell-Cell: How cells relate to neighbors (adjacency)
    - Logical: Structured (implicit) vs Unstructured (explicit) connectivity
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("VTK Connectivity Demo: Multiple Levels of Connectivity")
        self.resize(1100, 800)

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
        self.combo_box.setMinimumWidth(350)
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
        dataset, description = create_func(points)

        # Update info label
        num_points = dataset.GetNumberOfPoints()
        num_cells = dataset.GetNumberOfCells()

        # Get cell types present
        cell_types = set()
        for i in range(num_cells):
            cell_types.add(dataset.GetCellTypeName(dataset.GetCellType(i)))
        cell_type_str = ", ".join(sorted(cell_types))

        # Add dataset type info
        dataset_type = dataset.GetClassName()

        info_text = (
            f"<b>Mode: {mode}</b><br/>"
            f"{description}<br/><br/>"
            f"<b>Dataset:</b> {dataset_type} | "
            f"<b>Points:</b> {num_points} | <b>Cells:</b> {num_cells} | "
            f"<b>Cell Types:</b> {cell_type_str}"
        )
        self.info_label.setText(info_text)

        # Create surface mapper with edge visibility
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(dataset)

        # For adjacency and mixed type demos, use cell colors
        if dataset.GetCellData().GetScalars():
            surface_mapper.SetScalarModeToUseCellData()
            surface_mapper.SetColorModeToMapScalars()
            # Create a distinct color lookup table
            lut = vtk.vtkLookupTable()
            lut.SetNumberOfColors(8)
            lut.SetHueRange(0.0, 0.7)
            lut.Build()
            surface_mapper.SetLookupTable(lut)

        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)
        self.surface_actor.GetProperty().SetOpacity(0.7)
        self.surface_actor.GetProperty().SetEdgeVisibility(1)
        self.surface_actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.4)
        self.surface_actor.GetProperty().SetLineWidth(3)

        # Create explicit edge extraction for better visibility
        edges = vtk.vtkExtractEdges()
        edges.SetInputData(dataset)

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
        glyph.SetInputData(dataset)
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
    print("\n" + "=" * 75)
    print("VTK CONNECTIVITY: Understanding Multiple Levels of Connectivity")
    print("=" * 75)
    print("""
┌─────────────────────────────────────────────────────────────────────────┐
│                    LEVELS OF CONNECTIVITY IN VTK                        │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  1️⃣  POINT-CELL CONNECTIVITY (Basic Level)                              │
│  ─────────────────────────────────────────                              │
│  "Which points belong to this cell?"                                    │
│  • Triangle: 3 points   • Quad: 4 points   • Hex: 8 points              │
│  This defines the fundamental cell structure.                           │
│                                                                         │
│  2️⃣  CELL-CELL CONNECTIVITY (Adjacency)                                  │
│  ─────────────────────────────────────                                  │
│  "Which cells touch this cell?"                                         │
│  Cells may share: faces, edges, or vertices                             │
│  Important for: gradients, PDE solvers, neighborhood queries            │
│                                                                         │
│  3️⃣  LOGICAL CONNECTIVITY (Index Structure)                              │
│  ─────────────────────────────────────────                              │
│  "Is there a predictable pattern to neighbors?"                         │
│  • Structured grids: neighbors at (i±1, j±1, k±1) - IMPLICIT            │
│  • Unstructured grids: arbitrary - must store EXPLICITLY                │
│  This is why structured grids are FASTER!                               │
│                                                                         │
│  4️⃣  TOPOLOGICAL CONSTRAINTS                                             │
│  ────────────────────────────                                           │
│  "What rules govern the mesh?"                                          │
│  • StructuredGrid: fixed topology, regular pattern                      │
│  • UnstructuredGrid: arbitrary topology, mixed cell types               │
│                                                                         │
│  ─────────────────────────────────────────────────────────────────────  │
│                                                                         │
│  WHY THIS MATTERS:                                                      │
│                                                                         │
│  • CFD/FEA: All connectivity levels affect solver behavior              │
│  • Performance: Structured = fast, Unstructured = flexible              │
│  • Algorithms: Gradient, interpolation depend on cell neighbors         │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘

This demonstration covers:
- Point-Cell: How points form cells (vertices, lines, triangles, quads, hexes)
- Cell-Cell: How cells share faces/edges/vertices (adjacency)
- Logical: Structured (implicit) vs Unstructured (explicit) connectivity
""")


def main():
    """Launch the VTK Connectivity Demo application."""
    print_educational_summary()

    print("\nStarting interactive demonstration...")
    print("Explore all levels of connectivity:")
    print("  - Point-Cell: How points form different cell types")
    print("  - Cell-Cell: How cells share faces/edges/vertices")
    print("  - Logical: Structured (implicit) vs Unstructured (explicit)")
    print("-" * 75)

    app = QApplication(sys.argv)
    window = ConnectivityDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
