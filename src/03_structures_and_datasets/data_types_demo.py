"""
VTK Data Types Demo: Interactive Visualization of VTK Dataset Types

This module provides an interactive demonstration of all VTK dataset types using
a PyQt6-based interface with a combo box for selecting different data types.

Each dataset type is displayed clearly in the center of the view with:
- Visible vertices (as spheres)
- Visible edges for better understanding of topology
- Data type information displayed as text overlay

Dataset Types Covered:
---------------------
1. vtkImageData: Handles regular, rectilinear grid data
   - Uniform spacing in all directions
   - Most memory-efficient grid type
   - Used for volumetric data, images, uniform grids

2. vtkRectilinearGrid: Manages grid data with varying spacing
   - Axis-aligned cells with non-uniform spacing
   - Defined by 1D coordinate arrays
   - Used for stretched grids, boundary layer meshes

3. vtkStructuredGrid: Deals with structured data points in 3D space
   - Regular i-j-k topology with curvilinear coordinates
   - Grid lines can be curved
   - Used for body-fitted meshes, O-grids, C-grids

4. vtkPolyData: Specializes in representing polygonal data
   - Vertices, lines, polygons, and triangle strips
   - Primary format for surface meshes
   - Used for STL, OBJ, and surface visualization

5. vtkUnstructuredGrid: Ideal for representing complex, irregular grid data
   - Explicit connectivity, any cell type
   - Maximum geometric flexibility
   - Used for FEM/CFD with complex geometries

Usage:
------
Run this script directly to launch the interactive viewer:
    python data_types_demo.py

Select different data types from the combo box to visualize each one.
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


# Visual constants for rendering
VERTEX_SPHERE_RADIUS = 0.08  # Radius of spheres showing vertices
CAMERA_AZIMUTH = 30  # Initial camera azimuth angle in degrees
CAMERA_ELEVATION = 20  # Initial camera elevation angle in degrees


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
    image_data.SetDimensions(2, 2, 2)  # 2x2x2 points = 1 voxel
    image_data.SetSpacing(1.0, 1.0, 1.0)
    image_data.SetOrigin(-0.5, -0.5, -0.5)  # Center at origin

    return (
        image_data,
        "vtkImageData",
        "Regular grid with uniform spacing. "
        "Most memory-efficient for volumetric data and images.",
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
    x_coords.InsertNextValue(-0.6)
    x_coords.InsertNextValue(0.4)  # Non-uniform in X

    y_coords = vtk.vtkFloatArray()
    y_coords.InsertNextValue(-0.5)
    y_coords.InsertNextValue(0.5)

    z_coords = vtk.vtkFloatArray()
    z_coords.InsertNextValue(-0.4)
    z_coords.InsertNextValue(0.6)  # Non-uniform in Z

    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(2, 2, 2)
    rgrid.SetXCoordinates(x_coords)
    rgrid.SetYCoordinates(y_coords)
    rgrid.SetZCoordinates(z_coords)

    return (
        rgrid,
        "vtkRectilinearGrid",
        "Axis-aligned grid with non-uniform spacing. "
        "Defined by 1D coordinate arrays for each axis.",
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
    grid.SetDimensions(2, 2, 2)

    points = vtk.vtkPoints()
    # Create a slightly warped hexahedron to show curvilinear nature
    # Bottom face (z=0)
    points.InsertNextPoint(-0.6, -0.5, -0.5)  # (0,0,0)
    points.InsertNextPoint(0.6, -0.4, -0.5)  # (1,0,0) - slightly warped
    points.InsertNextPoint(-0.5, 0.5, -0.4)  # (0,1,0) - slightly warped
    points.InsertNextPoint(0.5, 0.6, -0.5)  # (1,1,0) - slightly warped
    # Top face (z=1)
    points.InsertNextPoint(-0.5, -0.6, 0.5)  # (0,0,1) - slightly warped
    points.InsertNextPoint(0.5, -0.5, 0.6)  # (1,0,1) - slightly warped
    points.InsertNextPoint(-0.6, 0.5, 0.5)  # (0,1,1) - slightly warped
    points.InsertNextPoint(0.6, 0.5, 0.5)  # (1,1,1)

    grid.SetPoints(points)

    return (
        grid,
        "vtkStructuredGrid",
        "Curvilinear structured grid. Regular i-j-k topology "
        "but points can be at any position (curved grid lines).",
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
    points.InsertNextPoint(-0.8, -0.6, 0.0)
    points.InsertNextPoint(0.8, -0.5, 0.0)
    points.InsertNextPoint(0.7, 0.7, 0.0)
    points.InsertNextPoint(-0.6, 0.6, 0.0)

    # Create a quad cell
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, 0)
    quad.GetPointIds().SetId(1, 1)
    quad.GetPointIds().SetId(2, 2)
    quad.GetPointIds().SetId(3, 3)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(quad)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    return (
        polydata,
        "vtkPolyData",
        "Polygonal surface mesh. Primary format for vertices, lines, "
        "and polygons. Used for STL/OBJ files and surface visualization.",
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
    # Hexahedron vertices (VTK ordering)
    points.InsertNextPoint(-0.5, -0.5, -0.5)  # 0
    points.InsertNextPoint(0.5, -0.5, -0.5)  # 1
    points.InsertNextPoint(0.5, 0.5, -0.5)  # 2
    points.InsertNextPoint(-0.5, 0.5, -0.5)  # 3
    points.InsertNextPoint(-0.5, -0.5, 0.5)  # 4
    points.InsertNextPoint(0.5, -0.5, 0.5)  # 5
    points.InsertNextPoint(0.5, 0.5, 0.5)  # 6
    points.InsertNextPoint(-0.5, 0.5, 0.5)  # 7

    hexa = vtk.vtkHexahedron()
    for i in range(8):
        hexa.GetPointIds().SetId(i, i)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.InsertNextCell(hexa.GetCellType(), hexa.GetPointIds())

    return (
        ugrid,
        "vtkUnstructuredGrid",
        "Flexible mesh with explicit connectivity. Supports any cell type "
        "and can mix different cell types. Used for complex FEM/CFD geometries.",
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

        header_layout.addStretch()
        layout.addLayout(header_layout)

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

        # Display first data type
        self.on_data_type_changed(self.combo_box.currentText())

    def on_data_type_changed(self, data_type_name):
        """Handle combo box selection change."""
        if data_type_name not in DATA_TYPE_CREATORS:
            return

        # Create the selected data type
        data, display_name, description = DATA_TYPE_CREATORS[data_type_name]()

        # Get dataset info
        num_points = data.GetNumberOfPoints()
        num_cells = data.GetNumberOfCells()

        # Update info label
        info_text = (
            f"<b>{display_name}</b><br/>"
            f"{description}<br/>"
            f"Points: {num_points} | Cells: {num_cells}"
        )
        self.info_label.setText(info_text)

        # Create surface mapper - handle both polydata and other dataset types
        if isinstance(data, vtk.vtkPolyData):
            surface_mapper = vtk.vtkPolyDataMapper()
            surface_mapper.SetInputData(data)
        else:
            surface_mapper = vtk.vtkDataSetMapper()
            surface_mapper.SetInputData(data)

        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)
        self.surface_actor.GetProperty().SetOpacity(0.8)
        self.surface_actor.GetProperty().SetEdgeVisibility(1)
        self.surface_actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
        self.surface_actor.GetProperty().SetLineWidth(2)

        # Create point glyph visualization
        # For all data types, we need to get points for glyph display
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

        self.point_actor.SetMapper(point_mapper)
        self.point_actor.GetProperty().SetColor(0.9, 0.2, 0.2)

        # Reset camera to show the data clearly
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
    """Launch the VTK Data Types Demo application."""
    app = QApplication(sys.argv)
    window = DataTypesDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
