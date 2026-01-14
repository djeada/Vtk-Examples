"""
VTK Data Types Demo: Interactive Visualization of VTK Dataset Types

This module provides an interactive demonstration of all VTK dataset types using
a PyQt6-based interface with a combo box for selecting different data types.

Each dataset type is displayed clearly in the center of the view with:
- Visible vertices (as spheres)
- Visible edges for better understanding of topology
- Connectivity graph between adjacent cells
- Data type and field information displayed as text overlay

Dataset Types Hierarchy (Primitive -> Powerful):
-----------------------------------------------
The data types form a hierarchy from most constrained to most flexible:

1. vtkImageData (MOST CONSTRAINED / MOST EFFICIENT):
   - Uniform spacing in all directions
   - No explicit point storage (computed from origin + spacing)
   - Fastest access, lowest memory footprint
   - Use when: Regular grids, medical imaging, volumetric data

2. vtkRectilinearGrid (+VARIABLE SPACING):
   - Adds non-uniform spacing via 1D coordinate arrays
   - Still axis-aligned cells only
   - Use when: Stretched grids, boundary layer clustering

3. vtkStructuredGrid (+CURVED GEOMETRY):
   - Adds curvilinear coordinates (points at any position)
   - Grid lines can be curved (body-fitted meshes)
   - Stores all 3D points explicitly
   - Use when: O-grids around airfoils, pipe flows, curved boundaries

4. vtkPolyData (SURFACE SPECIALIST):
   - Vertices, lines, polygons, and triangle strips
   - Optimized for surface representation
   - Explicit connectivity but surface elements only
   - Use when: STL/OBJ files, surface visualization, rendering

5. vtkUnstructuredGrid (MOST POWERFUL / HIGHEST COST):
   - Any cell type (tetra, hexa, wedge, pyramid, polyhedra)
   - Mixed cell types in same mesh
   - Arbitrary connectivity
   - Use when: Complex FEM/CFD geometries, adaptive meshes

Usage:
------
Run this script directly to launch the interactive viewer:
    python data_types_demo.py

Select different data types from the combo box to visualize each one.
"""

import math
import sys

import vtk
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QVBoxLayout,
    QWidget,
)
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


# Visual constants for rendering
VERTEX_SPHERE_RADIUS = 0.08  # Radius of spheres showing vertices
CAMERA_AZIMUTH = 30  # Initial camera azimuth angle in degrees
CAMERA_ELEVATION = 20  # Initial camera elevation angle in degrees
SCALAR_FIELD_NAME = "SyntheticField"
CELL_TYPE_FIELD_NAME = "CellType"
COLOR_MODE_POINT = "Point Scalar"
COLOR_MODE_CELL = "Cell Type"
CONNECTIVITY_LABEL = "Connectivity"
CELL_LABELS_LABEL = "Cell Labels"


def add_point_scalars(dataset, name=SCALAR_FIELD_NAME):
    """Attach a synthetic point scalar field so datasets render with meaningful coloring."""
    scalars = vtk.vtkFloatArray()
    scalars.SetName(name)

    for point_id in range(dataset.GetNumberOfPoints()):
        x, y, z = dataset.GetPoint(point_id)
        value = math.sqrt(x * x + y * y + z * z) + 0.35 * math.sin(3.0 * x) * math.cos(
            2.0 * y
        )
        scalars.InsertNextValue(value)

    dataset.GetPointData().SetScalars(scalars)


def add_cell_type_scalars(dataset, name=CELL_TYPE_FIELD_NAME):
    """Attach a cell-data array that stores the VTK cell type id for each cell."""
    cell_types = vtk.vtkIntArray()
    cell_types.SetName(name)

    for cell_id in range(dataset.GetNumberOfCells()):
        cell_types.InsertNextValue(dataset.GetCellType(cell_id))

    dataset.GetCellData().SetScalars(cell_types)


def get_cell_type_summary(dataset):
    """Return a summary of cell types present in the dataset."""
    counts = {}
    for cell_id in range(dataset.GetNumberOfCells()):
        cell_type_id = dataset.GetCellType(cell_id)
        cell_name = vtk.vtkCellTypes.GetClassNameFromTypeId(cell_type_id)
        counts[cell_name] = counts.get(cell_name, 0) + 1
    return counts


def format_cell_type_summary(counts):
    """Format a cell type summary dict for UI display."""
    if not counts:
        return "None"
    parts = [f"{name} x{count}" for name, count in sorted(counts.items())]
    return ", ".join(parts)


def build_cell_type_lookup(cell_type_ids):
    """Build a categorical lookup table for cell types."""
    if not cell_type_ids:
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(1)
        lut.Build()
        return lut, (0, 0), []

    max_id = max(cell_type_ids)
    unique_ids = sorted(set(cell_type_ids))

    lut = vtk.vtkLookupTable()
    lut.SetIndexedLookup(1)
    lut.SetNumberOfTableValues(max_id + 1)

    colors = vtk.vtkColorSeries()
    colors.SetColorScheme(vtk.vtkColorSeries.BREWER_QUALITATIVE_SET3)
    color_count = max(colors.GetNumberOfColors(), 1)

    for idx, cell_type_id in enumerate(unique_ids):
        color = colors.GetColor(idx % color_count)
        lut.SetTableValue(
            cell_type_id,
            color.GetRed() / 255.0,
            color.GetGreen() / 255.0,
            color.GetBlue() / 255.0,
            1.0,
        )
        lut.SetAnnotation(
            cell_type_id, vtk.vtkCellTypes.GetClassNameFromTypeId(cell_type_id)
        )

    lut.Build()
    return lut, (0, max_id), unique_ids


def build_cell_adjacency_polydata(dataset):
    """Build a line-based connectivity graph between adjacent cells."""
    num_cells = dataset.GetNumberOfCells()
    if num_cells < 2:
        return None, 0

    if hasattr(dataset, "BuildLinks"):
        dataset.BuildLinks()

    centers = vtk.vtkCellCenters()
    centers.SetInputData(dataset)
    centers.Update()
    center_output = centers.GetOutput()
    center_points = vtk.vtkPoints()
    center_points.DeepCopy(center_output.GetPoints())

    lines = vtk.vtkCellArray()
    cell_points = vtk.vtkIdList()
    point_cells = vtk.vtkIdList()
    edge_count = 0

    for cell_id in range(num_cells):
        dataset.GetCellPoints(cell_id, cell_points)
        neighbors = set()
        for idx in range(cell_points.GetNumberOfIds()):
            point_id = cell_points.GetId(idx)
            dataset.GetPointCells(point_id, point_cells)
            for neighbor_idx in range(point_cells.GetNumberOfIds()):
                neighbor_id = point_cells.GetId(neighbor_idx)
                if neighbor_id != cell_id:
                    neighbors.add(neighbor_id)
        for neighbor_id in neighbors:
            if neighbor_id > cell_id:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, cell_id)
                line.GetPointIds().SetId(1, neighbor_id)
                lines.InsertNextCell(line)
                edge_count += 1

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(center_points)
    polydata.SetLines(lines)
    return polydata, edge_count


def build_cell_id_label_actor(dataset):
    """Create a 2D label actor for cell IDs."""
    centers = vtk.vtkCellCenters()
    centers.SetInputData(dataset)
    centers.Update()
    center_output = centers.GetOutput()

    labels = vtk.vtkIntArray()
    labels.SetName("Label")
    labels.SetNumberOfComponents(1)
    num_cells = dataset.GetNumberOfCells()
    labels.SetNumberOfTuples(num_cells)

    for i in range(num_cells):
        labels.SetTuple1(i, i)

    center_output.GetPointData().SetScalars(labels)

    label_mapper = vtk.vtkLabeledDataMapper()
    label_mapper.SetInputData(center_output)
    label_mapper.SetLabelModeToLabelScalars()
    label_mapper.GetLabelTextProperty().SetColor(1.0, 1.0, 1.0)
    label_mapper.GetLabelTextProperty().SetFontSize(12)

    label_actor = vtk.vtkActor2D()
    label_actor.SetMapper(label_mapper)
    return label_actor


def get_structured_dimensions(dataset):
    """Return structured dimensions as a tuple when available."""
    if hasattr(dataset, "GetDimensions"):
        try:
            return dataset.GetDimensions()
        except TypeError:
            dims = [0, 0, 0]
            dataset.GetDimensions(dims)
            return tuple(dims)
    if hasattr(dataset, "GetExtent"):
        extent = dataset.GetExtent()
        return (
            extent[1] - extent[0] + 1,
            extent[3] - extent[2] + 1,
            extent[5] - extent[4] + 1,
        )
    return None


def create_image_data():
    """
    Create a vtkImageData representing a single voxel (3D cell).

    vtkImageData is the simplest and most memory-efficient grid type:
    - Uniform spacing in all directions
    - Origin and spacing define point positions
    - No explicit point storage needed

    Returns:
        tuple: (vtkImageData, display_name, description)
    """
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(5, 4, 3)
    image_data.SetSpacing(0.35, 0.35, 0.5)
    image_data.SetOrigin(-0.7, -0.6, -0.5)  # Center the grid at origin

    add_point_scalars(image_data)

    return (
        image_data,
        "vtkImageData (Most Constrained)",
        "⚡ SIMPLEST: Uniform spacing only. No explicit points needed--positions computed "
        "from origin + spacing. Fastest access, lowest memory. Ideal for regular volumetric grids.",
    )


def create_rectilinear_grid():
    """
    Create a vtkRectilinearGrid representing a single cell with varying spacing.

    vtkRectilinearGrid allows non-uniform spacing:
    - Axis-aligned cells
    - Spacing varies in each direction via 1D coordinate arrays
    - More flexible than ImageData, less than StructuredGrid

    Returns:
        tuple: (vtkRectilinearGrid, display_name, description)
    """
    # Create coordinate arrays with varying spacing
    x_coords = vtk.vtkFloatArray()
    for value in (-0.9, -0.4, 0.2, 0.8):
        x_coords.InsertNextValue(value)

    y_coords = vtk.vtkFloatArray()
    for value in (-0.7, -0.2, 0.5):
        y_coords.InsertNextValue(value)

    z_coords = vtk.vtkFloatArray()
    for value in (-0.6, -0.1, 0.4):
        z_coords.InsertNextValue(value)

    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(
        x_coords.GetNumberOfTuples(),
        y_coords.GetNumberOfTuples(),
        z_coords.GetNumberOfTuples(),
    )
    rgrid.SetXCoordinates(x_coords)
    rgrid.SetYCoordinates(y_coords)
    rgrid.SetZCoordinates(z_coords)

    add_point_scalars(rgrid)

    return (
        rgrid,
        "vtkRectilinearGrid (+Variable Spacing)",
        "📊 ADDS: Non-uniform spacing via 1D coordinate arrays. Still axis-aligned cells. "
        "Great for stretched meshes and boundary-layer clustering.",
    )


def create_structured_grid():
    """
    Create a vtkStructuredGrid representing a single curved cell.

    vtkStructuredGrid has curvilinear coordinates:
    - Regular i-j-k topology
    - Points can be at any position
    - Grid lines can be curved (non-axis-aligned)

    Returns:
        tuple: (vtkStructuredGrid, display_name, description)
    """
    grid = vtk.vtkStructuredGrid()
    nx, ny, nz = 4, 3, 3
    grid.SetDimensions(nx, ny, nz)

    points = vtk.vtkPoints()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                fx = i / (nx - 1)
                fy = j / (ny - 1)
                fz = k / (nz - 1)
                x = (fx - 0.5) * 1.6
                y = (fy - 0.5) * 1.2
                z = (fz - 0.5) * 1.0
                x += 0.15 * math.sin(math.pi * fy) * math.cos(math.pi * fz)
                y += 0.12 * math.sin(math.pi * fx) * math.cos(math.pi * fz)
                z += 0.25 * math.sin(math.pi * fx) * math.cos(math.pi * fy)
                points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)

    add_point_scalars(grid)

    return (
        grid,
        "vtkStructuredGrid (+Curved Geometry)",
        "🔄 ADDS: Curvilinear coordinates--points at ANY position with curved grid lines. "
        "Enables body-fitted meshes around airfoils, cylinders. Stores all 3D points explicitly.",
    )


def create_poly_data():
    """
    Create a vtkPolyData representing a simple polygon (quad).

    vtkPolyData is the primary format for surface meshes:
    - Contains vertices, lines, polygons, triangle strips
    - Shared points with explicit connectivity
    - Used for STL, OBJ, and surface visualization

    Returns:
        tuple: (vtkPolyData, display_name, description)
    """
    points = vtk.vtkPoints()
    points.InsertNextPoint(-0.9, -0.6, 0.0)  # 0
    points.InsertNextPoint(-0.1, -0.7, 0.2)  # 1
    points.InsertNextPoint(0.8, -0.2, 0.0)  # 2
    points.InsertNextPoint(-0.6, 0.4, 0.1)  # 3
    points.InsertNextPoint(0.2, 0.6, -0.1)  # 4
    points.InsertNextPoint(0.7, 0.8, 0.0)  # 5
    points.InsertNextPoint(-0.1, 0.9, 0.2)  # 6

    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, 0)
    quad.GetPointIds().SetId(1, 1)
    quad.GetPointIds().SetId(2, 2)
    quad.GetPointIds().SetId(3, 3)

    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 3)
    triangle.GetPointIds().SetId(1, 2)
    triangle.GetPointIds().SetId(2, 4)

    polys = vtk.vtkCellArray()
    polys.InsertNextCell(quad)
    polys.InsertNextCell(triangle)

    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(3)
    polyline.GetPointIds().SetId(0, 4)
    polyline.GetPointIds().SetId(1, 5)
    polyline.GetPointIds().SetId(2, 6)

    lines = vtk.vtkCellArray()
    lines.InsertNextCell(polyline)

    vertex = vtk.vtkVertex()
    vertex.GetPointIds().SetId(0, 6)

    verts = vtk.vtkCellArray()
    verts.InsertNextCell(vertex)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polys)
    polydata.SetLines(lines)
    polydata.SetVerts(verts)

    add_point_scalars(polydata)

    return (
        polydata,
        "vtkPolyData (Surface Specialist)",
        "🎨 SURFACE-FOCUSED: Optimized for polygonal surfaces (verts, lines, polys, strips). "
        "Standard for STL/OBJ. Explicit connectivity but limited to surface elements only.",
    )


def create_unstructured_grid():
    """
    Create a vtkUnstructuredGrid representing a single hexahedron.

    vtkUnstructuredGrid has maximum flexibility:
    - Explicit connectivity stored per cell
    - Any cell type can be used
    - Can mix different cell types

    Returns:
        tuple: (vtkUnstructuredGrid, display_name, description)
    """
    points = vtk.vtkPoints()
    # Hexahedron (0-7)
    points.InsertNextPoint(-1.0, -0.4, -0.3)
    points.InsertNextPoint(-0.2, -0.4, -0.3)
    points.InsertNextPoint(-0.2, 0.4, -0.3)
    points.InsertNextPoint(-1.0, 0.4, -0.3)
    points.InsertNextPoint(-1.0, -0.4, 0.3)
    points.InsertNextPoint(-0.2, -0.4, 0.3)
    points.InsertNextPoint(-0.2, 0.4, 0.3)
    points.InsertNextPoint(-1.0, 0.4, 0.3)
    # Pyramid apex (8)
    points.InsertNextPoint(-0.6, 0.0, 0.9)
    # Tetra apex (9)
    points.InsertNextPoint(0.4, 0.1, 0.0)
    # Wedge extrude (10-12)
    points.InsertNextPoint(-1.0, -0.4, -0.9)
    points.InsertNextPoint(-0.2, -0.4, -0.9)
    points.InsertNextPoint(-0.2, 0.4, -0.9)

    hexa = vtk.vtkHexahedron()
    for i in range(8):
        hexa.GetPointIds().SetId(i, i)

    pyramid = vtk.vtkPyramid()
    for i, point_id in enumerate((4, 5, 6, 7, 8)):
        pyramid.GetPointIds().SetId(i, point_id)

    tetra = vtk.vtkTetra()
    for i, point_id in enumerate((2, 6, 5, 9)):
        tetra.GetPointIds().SetId(i, point_id)

    wedge = vtk.vtkWedge()
    for i, point_id in enumerate((0, 1, 2, 10, 11, 12)):
        wedge.GetPointIds().SetId(i, point_id)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(hexa.GetCellType(), hexa.GetPointIds())
    ugrid.InsertNextCell(pyramid.GetCellType(), pyramid.GetPointIds())
    ugrid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())
    ugrid.InsertNextCell(wedge.GetCellType(), wedge.GetPointIds())

    add_point_scalars(ugrid)

    return (
        ugrid,
        "vtkUnstructuredGrid (Most Powerful)",
        "🚀 MAXIMUM FLEXIBILITY: Any cell type, mixed cells, arbitrary connectivity. "
        "Essential for complex FEM/CFD. Highest memory cost but handles ANY geometry.",
    )


# Dictionary mapping data type names to their creation functions
DATA_TYPE_CREATORS = {
    "vtkImageData": create_image_data,
    "vtkRectilinearGrid": create_rectilinear_grid,
    "vtkStructuredGrid": create_structured_grid,
    "vtkPolyData": create_poly_data,
    "vtkUnstructuredGrid": create_unstructured_grid,
}


class DataTypesDemo(QMainWindow):
    """
    Interactive VTK Data Types Demo using PyQt6.

    This widget displays a combo box for selecting different VTK data types
    and renders the selected data type in a VTK view with visible vertices and edges.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("VTK Data Types Demo")
        self.resize(900, 700)

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create header with combo box
        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Select Data Type:"))

        self.combo_box = QComboBox()
        self.combo_box.addItems(DATA_TYPE_CREATORS.keys())
        self.combo_box.setMinimumWidth(200)
        self.combo_box.currentTextChanged.connect(self.on_data_type_changed)
        header_layout.addWidget(self.combo_box)

        header_layout.addSpacing(12)
        header_layout.addWidget(QLabel("Color By:"))
        self.color_by_combo = QComboBox()
        self.color_by_combo.addItems([COLOR_MODE_POINT, COLOR_MODE_CELL])
        self.color_by_combo.currentTextChanged.connect(self.on_render_option_changed)
        header_layout.addWidget(self.color_by_combo)

        header_layout.addStretch()
        layout.addLayout(header_layout)

        # View options
        view_layout = QHBoxLayout()
        self.show_surface_cb = QCheckBox("Surface")
        self.show_surface_cb.setChecked(True)
        self.show_surface_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_surface_cb)

        self.show_edges_cb = QCheckBox("Edges")
        self.show_edges_cb.setChecked(True)
        self.show_edges_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_edges_cb)

        self.show_points_cb = QCheckBox("Points")
        self.show_points_cb.setChecked(True)
        self.show_points_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_points_cb)

        self.show_connectivity_cb = QCheckBox(CONNECTIVITY_LABEL)
        self.show_connectivity_cb.setChecked(True)
        self.show_connectivity_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_connectivity_cb)

        self.show_labels_cb = QCheckBox(CELL_LABELS_LABEL)
        self.show_labels_cb.setChecked(False)
        self.show_labels_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_labels_cb)

        self.show_scalar_bar_cb = QCheckBox("Scalar Bar")
        self.show_scalar_bar_cb.setChecked(True)
        self.show_scalar_bar_cb.toggled.connect(self.on_render_option_changed)
        view_layout.addWidget(self.show_scalar_bar_cb)

        view_layout.addStretch()
        layout.addLayout(view_layout)

        # Create info label
        self.info_label = QLabel()
        self.info_label.setStyleSheet(
            "QLabel { background-color: #2b2b2b; color: #ffffff; padding: 10px; "
            "border-radius: 5px; font-size: 12px; }"
        )
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.info_label.setWordWrap(True)
        layout.addWidget(self.info_label)

        # Create VTK widget
        self.vtk_widget = QVTKRenderWindowInteractor(central_widget)
        layout.addWidget(self.vtk_widget, stretch=1)

        # Set up renderer
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(0.1, 0.1, 0.15)
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)

        # Color map for scalar fields
        self.lookup_table = vtk.vtkLookupTable()
        self.lookup_table.SetNumberOfTableValues(256)
        self.lookup_table.SetHueRange(0.67, 0.0)
        self.lookup_table.Build()

        # Initialize actors
        self.surface_actor = vtk.vtkActor()
        self.point_actor = vtk.vtkActor()
        self.edge_actor = vtk.vtkActor()
        self.connectivity_actor = vtk.vtkActor()
        self.renderer.AddActor(self.surface_actor)
        self.renderer.AddActor(self.edge_actor)
        self.renderer.AddActor(self.point_actor)
        self.renderer.AddActor(self.connectivity_actor)

        # Scalar bar for field visualization
        self.scalar_bar = vtk.vtkScalarBarActor()
        self.scalar_bar.SetLookupTable(self.lookup_table)
        self.scalar_bar.SetTitle("Point Scalar")
        self.scalar_bar.SetNumberOfLabels(4)
        self.renderer.AddViewProp(self.scalar_bar)

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

        self.current_data = None
        self.current_display_name = ""
        self.current_description = ""
        self.label_actor = None

        # Display first data type
        self.on_data_type_changed(self.combo_box.currentText())

    def on_data_type_changed(self, data_type_name):
        """Handle combo box selection change."""
        if data_type_name not in DATA_TYPE_CREATORS:
            return

        # Create the selected data type
        data, display_name, description = DATA_TYPE_CREATORS[data_type_name]()
        add_cell_type_scalars(data)

        self.current_data = data
        self.current_display_name = display_name
        self.current_description = description
        self.update_scene(reset_camera=True)

    def on_render_option_changed(self):
        """Handle changes to rendering options without recreating datasets."""
        if self.current_data is None:
            return
        self.update_scene(reset_camera=False)

    def update_scene(self, reset_camera=False):
        """Update mappers, actors, and UI for the current dataset."""
        data = self.current_data
        if data is None:
            return

        # Dataset info
        num_points = data.GetNumberOfPoints()
        num_cells = data.GetNumberOfCells()
        point_scalars = data.GetPointData().GetScalars()
        cell_scalars = data.GetCellData().GetScalars()
        cell_type_summary = format_cell_type_summary(get_cell_type_summary(data))

        dims_info = ""
        if isinstance(
            data, (vtk.vtkImageData, vtk.vtkRectilinearGrid, vtk.vtkStructuredGrid)
        ):
            dims = get_structured_dimensions(data)
            if dims:
                dims_info = f"Dimensions: {dims[0]}x{dims[1]}x{dims[2]}"
        if isinstance(data, vtk.vtkImageData):
            spacing = data.GetSpacing()
            spacing_info = (
                f"Spacing: {spacing[0]:.2f}, {spacing[1]:.2f}, {spacing[2]:.2f}"
            )
        elif isinstance(data, vtk.vtkRectilinearGrid):
            spacing_info = "Spacing: non-uniform"
        else:
            spacing_info = ""

        color_mode = self.color_by_combo.currentText()
        use_cell_data = (
            color_mode == COLOR_MODE_CELL
            and cell_scalars is not None
            and cell_scalars.GetNumberOfTuples() > 0
        )
        cell_type_ids = []
        cell_lut = None
        cell_range = (0, 0)
        unique_cell_types = []

        if use_cell_data:
            cell_type_ids = [
                int(cell_scalars.GetTuple1(i))
                for i in range(cell_scalars.GetNumberOfTuples())
            ]
            cell_lut, cell_range, unique_cell_types = build_cell_type_lookup(
                cell_type_ids
            )

        if use_cell_data:
            active_field = CELL_TYPE_FIELD_NAME
            active_range = cell_range
        else:
            active_field = point_scalars.GetName() if point_scalars else "None"
            active_range = point_scalars.GetRange() if point_scalars else (0.0, 0.0)

        extra_info_parts = []
        if dims_info:
            extra_info_parts.append(dims_info)
        if spacing_info:
            extra_info_parts.append(spacing_info)
        extra_info = " | ".join(extra_info_parts)

        info_text = (
            f"<b>{self.current_display_name}</b><br/>"
            f"{self.current_description}<br/>"
            f"Points: {num_points} | Cells: {num_cells} | Color: {color_mode} | "
            f"Field: {active_field} ({active_range[0]:.2f} to {active_range[1]:.2f})"
        )
        if extra_info:
            info_text += f"<br/>{extra_info}"
        info_text += f"<br/>Cell types: {cell_type_summary}"

        # Select input for surface/edge filters
        self.shrink_filter = None
        if isinstance(data, vtk.vtkUnstructuredGrid):
            self.shrink_filter = vtk.vtkShrinkFilter()
            self.shrink_filter.SetInputData(data)
            self.shrink_filter.SetShrinkFactor(0.85)

        # Surface mapper
        if isinstance(data, vtk.vtkPolyData):
            surface_mapper = vtk.vtkPolyDataMapper()
        else:
            surface_mapper = vtk.vtkDataSetMapper()

        if self.shrink_filter:
            surface_mapper.SetInputConnection(self.shrink_filter.GetOutputPort())
        else:
            surface_mapper.SetInputData(data)

        if use_cell_data and cell_lut:
            surface_mapper.SetScalarModeToUseCellData()
            surface_mapper.SetLookupTable(cell_lut)
            surface_mapper.SetScalarRange(cell_range)
            surface_mapper.ScalarVisibilityOn()
        elif point_scalars:
            surface_mapper.SetScalarModeToUsePointData()
            surface_mapper.SetLookupTable(self.lookup_table)
            surface_mapper.SetScalarRange(active_range)
            surface_mapper.ScalarVisibilityOn()
        else:
            surface_mapper.ScalarVisibilityOff()

        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetOpacity(0.85)
        self.surface_actor.SetVisibility(self.show_surface_cb.isChecked())

        # Point glyph visualization
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(VERTEX_SPHERE_RADIUS)
        sphere_source.SetThetaResolution(16)
        sphere_source.SetPhiResolution(16)

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(data)
        glyph.SetSourceConnection(sphere_source.GetOutputPort())
        glyph.SetScaleModeToDataScalingOff()

        point_mapper = vtk.vtkPolyDataMapper()
        point_mapper.SetInputConnection(glyph.GetOutputPort())
        if use_cell_data:
            point_mapper.ScalarVisibilityOff()
            self.point_actor.GetProperty().SetColor(0.9, 0.2, 0.2)
        elif point_scalars:
            point_mapper.SetLookupTable(self.lookup_table)
            point_mapper.SetScalarRange(active_range)
            point_mapper.ScalarVisibilityOn()
        else:
            point_mapper.ScalarVisibilityOff()

        self.point_actor.SetMapper(point_mapper)
        self.point_actor.GetProperty().SetOpacity(0.9)
        self.point_actor.SetVisibility(self.show_points_cb.isChecked())

        # Connectivity graph between adjacent cells
        if self.show_connectivity_cb.isChecked():
            connectivity_data, edge_count = build_cell_adjacency_polydata(data)
            if connectivity_data and edge_count > 0:
                connectivity_mapper = vtk.vtkPolyDataMapper()
                connectivity_mapper.SetInputData(connectivity_data)
                connectivity_mapper.ScalarVisibilityOff()

                self.connectivity_actor.SetMapper(connectivity_mapper)
                self.connectivity_actor.GetProperty().SetColor(0.95, 0.7, 0.1)
                self.connectivity_actor.GetProperty().SetLineWidth(2.2)
                self.connectivity_actor.SetVisibility(True)
                info_text += (
                    f"<br/>Connectivity links: {edge_count} (shared-point adjacency)"
                )
            else:
                self.connectivity_actor.SetVisibility(False)
                info_text += "<br/>Connectivity links: none"
        else:
            self.connectivity_actor.SetVisibility(False)

        # Edges to highlight topology
        self.edge_filter = vtk.vtkExtractEdges()
        if self.shrink_filter:
            self.edge_filter.SetInputConnection(self.shrink_filter.GetOutputPort())
        else:
            self.edge_filter.SetInputData(data)

        edge_mapper = vtk.vtkPolyDataMapper()
        edge_mapper.SetInputConnection(self.edge_filter.GetOutputPort())
        edge_mapper.ScalarVisibilityOff()

        self.edge_actor.SetMapper(edge_mapper)
        self.edge_actor.GetProperty().SetColor(0.95, 0.95, 0.95)
        self.edge_actor.GetProperty().SetLineWidth(2.0)
        self.edge_actor.SetVisibility(self.show_edges_cb.isChecked())

        # Cell ID labels for explicit connectivity
        if self.show_labels_cb.isChecked():
            new_label_actor = build_cell_id_label_actor(data)
            if self.label_actor is not None:
                self.renderer.RemoveViewProp(self.label_actor)
            self.label_actor = new_label_actor
            self.renderer.AddViewProp(self.label_actor)
        elif self.label_actor is not None:
            self.renderer.RemoveViewProp(self.label_actor)
            self.label_actor = None

        # Scalar bar configuration
        if use_cell_data and cell_lut:
            self.scalar_bar.SetLookupTable(cell_lut)
            self.scalar_bar.SetTitle(CELL_TYPE_FIELD_NAME)
            self.scalar_bar.SetDrawAnnotations(1)
            self.scalar_bar.SetNumberOfLabels(max(len(unique_cell_types), 1))
        else:
            self.scalar_bar.SetLookupTable(self.lookup_table)
            self.scalar_bar.SetTitle(active_field)
            self.scalar_bar.SetDrawAnnotations(0)
            self.scalar_bar.SetNumberOfLabels(4)
        self.scalar_bar.SetVisibility(self.show_scalar_bar_cb.isChecked())

        self.info_label.setText(info_text)

        if reset_camera:
            self.renderer.ResetCamera()
            camera = self.renderer.GetActiveCamera()
            camera.Azimuth(CAMERA_AZIMUTH)
            camera.Elevation(CAMERA_ELEVATION)
            self.renderer.ResetCameraClippingRange()

        self.vtk_widget.GetRenderWindow().Render()

    def closeEvent(self, event):
        """Clean up VTK resources on close."""
        self.vtk_widget.Finalize()
        super().closeEvent(event)


def main():
    """Launch the VTK Data Types Demo application."""
    app = QApplication(sys.argv)
    window = DataTypesDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
