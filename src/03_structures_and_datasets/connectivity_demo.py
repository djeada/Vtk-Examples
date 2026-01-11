"""
VTK Connectivity Filter Demo: Interactive Visualization of Connected Regions

This module provides an interactive demonstration of VTK's connectivity filter
using a PyQt6-based interface with a combo box for selecting different extraction modes.

The connectivity filter (vtkConnectivityFilter) is used to extract connected regions
from datasets. A region is a collection of cells that are connected through shared
points or edges.

Connectivity Concepts:
---------------------
- **Connected Region**: A group of cells where any cell can be reached from any
  other cell in the group by traversing shared points or edges.
- **Disconnected Regions**: Multiple separate groups that have no shared points
  or edges between them.

Extraction Modes:
----------------
1. **All Regions**: Extracts all connected regions and colors them differently.
   Useful for identifying how many separate objects exist in the data.

2. **Largest Region**: Extracts only the largest connected region (by cell count).
   Useful for filtering out noise or small artifacts.

3. **Closest Region**: Extracts the region closest to a specified point.
   Useful for selecting a specific region based on spatial location.

4. **Specified Region**: Extracts a specific region by its ID number.
   Useful when you know which region you want to work with.

5. **Point Seeded Region**: Extracts the region containing a specified seed point.
   Useful for interactive selection or known point-based queries.

6. **Cell Seeded Region**: Extracts the region containing a specified seed cell.
   Useful for selection based on cell ID.

Real-World Applications:
-----------------------
- Medical imaging: Separating organs or tumors from background
- CAD: Identifying separate parts in an assembly
- Mesh processing: Removing small disconnected fragments
- Scientific visualization: Isolating specific structures

Usage:
------
Run this script directly to launch the interactive viewer:
    python connectivity_demo.py

Select different extraction modes from the combo box to see how
the connectivity filter operates on a multi-object dataset.
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
    QSpinBox,
)
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# Visual constants
CAMERA_AZIMUTH = 30
CAMERA_ELEVATION = 20


def create_multi_object_dataset():
    """
    Create a dataset with multiple disconnected objects.

    This creates several shapes (spheres, cube, cylinder, cone) at
    different positions, all combined into a single polydata object.
    The connectivity filter will be able to identify each as a
    separate region.

    Returns:
        vtkPolyData containing multiple disconnected objects
    """
    append_filter = vtk.vtkAppendPolyData()

    # Create first sphere (Region 0)
    sphere1 = vtk.vtkSphereSource()
    sphere1.SetCenter(-3, 0, 0)
    sphere1.SetRadius(1.0)
    sphere1.SetThetaResolution(20)
    sphere1.SetPhiResolution(20)
    sphere1.Update()
    append_filter.AddInputData(sphere1.GetOutput())

    # Create second sphere (Region 1)
    sphere2 = vtk.vtkSphereSource()
    sphere2.SetCenter(0, 3, 0)
    sphere2.SetRadius(0.7)
    sphere2.SetThetaResolution(20)
    sphere2.SetPhiResolution(20)
    sphere2.Update()
    append_filter.AddInputData(sphere2.GetOutput())

    # Create third sphere (Region 2) - smaller
    sphere3 = vtk.vtkSphereSource()
    sphere3.SetCenter(0, -2, 2)
    sphere3.SetRadius(0.5)
    sphere3.SetThetaResolution(16)
    sphere3.SetPhiResolution(16)
    sphere3.Update()
    append_filter.AddInputData(sphere3.GetOutput())

    # Create a cube (Region 3)
    cube = vtk.vtkCubeSource()
    cube.SetCenter(3, 0, 0)
    cube.SetXLength(1.5)
    cube.SetYLength(1.5)
    cube.SetZLength(1.5)
    cube.Update()

    # Triangulate the cube for consistency
    cube_triangles = vtk.vtkTriangleFilter()
    cube_triangles.SetInputData(cube.GetOutput())
    cube_triangles.Update()
    append_filter.AddInputData(cube_triangles.GetOutput())

    # Create a cylinder (Region 4)
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetCenter(0, -3, -1)
    cylinder.SetRadius(0.6)
    cylinder.SetHeight(2.0)
    cylinder.SetResolution(20)
    cylinder.Update()

    # Triangulate the cylinder
    cylinder_triangles = vtk.vtkTriangleFilter()
    cylinder_triangles.SetInputData(cylinder.GetOutput())
    cylinder_triangles.Update()
    append_filter.AddInputData(cylinder_triangles.GetOutput())

    # Create a cone (Region 5)
    cone = vtk.vtkConeSource()
    cone.SetCenter(-2, 2, 1.5)
    cone.SetRadius(0.8)
    cone.SetHeight(1.5)
    cone.SetResolution(20)
    cone.Update()

    cone_triangles = vtk.vtkTriangleFilter()
    cone_triangles.SetInputData(cone.GetOutput())
    cone_triangles.Update()
    append_filter.AddInputData(cone_triangles.GetOutput())

    # Merge all objects
    append_filter.Update()

    # Clean up to ensure proper point sharing within each object
    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(append_filter.GetOutputPort())
    clean.Update()

    return clean.GetOutput()


def apply_connectivity_filter(polydata, mode, region_id=0, seed_point=None):
    """
    Apply connectivity filter with the specified extraction mode.

    Args:
        polydata: Input vtkPolyData with multiple disconnected regions
        mode: Extraction mode string
        region_id: Region ID for "Specified Region" mode
        seed_point: Tuple (x, y, z) for "Closest Region" and "Seeded Region" modes

    Returns:
        Tuple of (output polydata, number of regions, description)
    """
    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputData(polydata)

    if mode == "All Regions":
        connectivity.SetExtractionModeToAllRegions()
        connectivity.ColorRegionsOn()
        description = (
            "Extracts ALL connected regions and assigns each a unique RegionId.\n"
            "ColorRegionsOn() enables automatic region coloring."
        )

    elif mode == "Largest Region":
        connectivity.SetExtractionModeToLargestRegion()
        description = (
            "Extracts only the LARGEST connected region (by cell count).\n"
            "Useful for removing small noise fragments."
        )

    elif mode == "Closest Region":
        connectivity.SetExtractionModeToClosestPointRegion()
        if seed_point:
            connectivity.SetClosestPoint(*seed_point)
        else:
            connectivity.SetClosestPoint(-3, 0, 0)  # Default: first sphere
        description = (
            f"Extracts the region closest to point {seed_point or '(-3, 0, 0)'}.\n"
            "Useful for spatial selection of regions."
        )

    elif mode == "Specified Region":
        connectivity.SetExtractionModeToSpecifiedRegions()
        connectivity.InitializeSpecifiedRegionList()
        connectivity.AddSpecifiedRegion(region_id)
        description = (
            f"Extracts the specified region with ID {region_id}.\n"
            "Region IDs are assigned based on cell count (largest = 0)."
        )

    elif mode == "Point Seeded Region":
        connectivity.SetExtractionModeToPointSeededRegions()
        # Use the first point of the input as seed
        connectivity.InitializeSeedList()
        connectivity.AddSeed(0)  # First point
        description = (
            "Extracts the region containing the seed point (point ID 0).\n"
            "The seed point is the first point in the dataset."
        )

    elif mode == "Cell Seeded Region":
        connectivity.SetExtractionModeToCellSeededRegions()
        connectivity.InitializeSeedList()
        connectivity.AddSeed(0)  # First cell
        description = (
            "Extracts the region containing the seed cell (cell ID 0).\n"
            "The seed cell is the first cell in the dataset."
        )

    else:
        connectivity.SetExtractionModeToAllRegions()
        connectivity.ColorRegionsOn()
        description = "Default: Extracting all regions."

    connectivity.Update()

    # Get number of regions
    num_regions = connectivity.GetNumberOfExtractedRegions()

    return connectivity.GetOutput(), num_regions, description


# Dictionary of extraction modes
EXTRACTION_MODES = [
    "All Regions",
    "Largest Region",
    "Closest Region",
    "Specified Region",
    "Point Seeded Region",
    "Cell Seeded Region",
]


class ConnectivityDemo(QMainWindow):
    """
    Interactive VTK Connectivity Filter Demo using PyQt6.

    This widget displays a combo box for selecting different extraction modes
    and renders the result of the connectivity filter on a multi-object dataset.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("VTK Connectivity Filter Demo")
        self.resize(1000, 750)

        # Create the multi-object dataset
        self.input_data = create_multi_object_dataset()

        # Create central widget and layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Create header with combo box
        header_layout = QHBoxLayout()
        header_layout.addWidget(QLabel("Extraction Mode:"))

        self.combo_box = QComboBox()
        self.combo_box.addItems(EXTRACTION_MODES)
        self.combo_box.setMinimumWidth(200)
        self.combo_box.currentTextChanged.connect(self.on_mode_changed)
        header_layout.addWidget(self.combo_box)

        # Add region ID spinner for "Specified Region" mode
        header_layout.addWidget(QLabel("Region ID:"))
        self.region_spinner = QSpinBox()
        self.region_spinner.setMinimum(0)
        self.region_spinner.setMaximum(10)
        self.region_spinner.setValue(0)
        self.region_spinner.valueChanged.connect(self.on_region_changed)
        header_layout.addWidget(self.region_spinner)

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
        self.info_label.setMinimumHeight(80)
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
        self.renderer.AddActor(self.surface_actor)

        # Add axes widget
        self.axes = vtk.vtkAxesActor()
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(self.axes)
        self.axes_widget.SetInteractor(self.vtk_widget)
        self.axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
        self.axes_widget.SetEnabled(1)
        self.axes_widget.InteractiveOff()

        # Create a lookup table for coloring regions
        self.lut = vtk.vtkLookupTable()
        self.lut.SetNumberOfTableValues(10)
        self.lut.SetTableValue(0, 0.9, 0.2, 0.2, 1.0)  # Red
        self.lut.SetTableValue(1, 0.2, 0.9, 0.2, 1.0)  # Green
        self.lut.SetTableValue(2, 0.2, 0.2, 0.9, 1.0)  # Blue
        self.lut.SetTableValue(3, 0.9, 0.9, 0.2, 1.0)  # Yellow
        self.lut.SetTableValue(4, 0.9, 0.2, 0.9, 1.0)  # Magenta
        self.lut.SetTableValue(5, 0.2, 0.9, 0.9, 1.0)  # Cyan
        self.lut.SetTableValue(6, 0.9, 0.5, 0.2, 1.0)  # Orange
        self.lut.SetTableValue(7, 0.5, 0.2, 0.9, 1.0)  # Purple
        self.lut.SetTableValue(8, 0.2, 0.9, 0.5, 1.0)  # Teal
        self.lut.SetTableValue(9, 0.9, 0.5, 0.5, 1.0)  # Pink
        self.lut.Build()

        # Initialize interactor
        self.vtk_widget.Initialize()

        # Display first mode
        self.on_mode_changed(self.combo_box.currentText())

    def on_mode_changed(self, mode):
        """Handle combo box selection change."""
        self.update_visualization(mode)

    def on_region_changed(self, value):
        """Handle region ID spinner change."""
        if self.combo_box.currentText() == "Specified Region":
            self.update_visualization("Specified Region")

    def update_visualization(self, mode):
        """Update the visualization based on the selected mode."""
        # Get region ID for specified region mode
        region_id = self.region_spinner.value()

        # Apply connectivity filter
        output_data, num_regions, description = apply_connectivity_filter(
            self.input_data, mode, region_id=region_id
        )

        # Update info label
        input_cells = self.input_data.GetNumberOfCells()
        input_points = self.input_data.GetNumberOfPoints()
        output_cells = output_data.GetNumberOfCells()
        output_points = output_data.GetNumberOfPoints()

        info_text = (
            f"<b>Mode: {mode}</b><br/>"
            f"{description}<br/><br/>"
            f"<b>Input:</b> {input_points} points, {input_cells} cells | "
            f"<b>Output:</b> {output_points} points, {output_cells} cells | "
            f"<b>Regions found:</b> {num_regions}"
        )
        self.info_label.setText(info_text)

        # Create mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(output_data)

        # For "All Regions" mode, color by region ID
        if mode == "All Regions" and output_data.GetPointData().GetScalars():
            mapper.SetScalarModeToUsePointData()
            mapper.SetScalarRange(0, num_regions - 1)
            mapper.SetLookupTable(self.lut)
            mapper.ScalarVisibilityOn()
        else:
            mapper.ScalarVisibilityOff()

        self.surface_actor.SetMapper(mapper)

        # Set appearance for non-colored modes
        if mode != "All Regions":
            self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)

        self.surface_actor.GetProperty().SetOpacity(1.0)
        self.surface_actor.GetProperty().SetEdgeVisibility(0)

        # Reset camera for first render, otherwise just render
        if not hasattr(self, "_camera_initialized"):
            self.renderer.ResetCamera()
            camera = self.renderer.GetActiveCamera()
            camera.Azimuth(CAMERA_AZIMUTH)
            camera.Elevation(CAMERA_ELEVATION)
            self.renderer.ResetCameraClippingRange()
            self._camera_initialized = True

        # Render the scene
        self.vtk_widget.GetRenderWindow().Render()

    def closeEvent(self, event):
        """Clean up VTK resources on close."""
        self.vtk_widget.Finalize()
        super().closeEvent(event)


def print_educational_summary():
    """Print educational information about connectivity in VTK."""
    print("\n" + "=" * 70)
    print("VTK CONNECTIVITY FILTER: Educational Summary")
    print("=" * 70)
    print("""
┌─────────────────────────────────────────────────────────────────────┐
│                    CONNECTIVITY IN VTK                              │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  What is Connectivity?                                              │
│  ─────────────────────                                              │
│  Connectivity refers to how cells are connected through shared      │
│  points. A connected region is a group of cells where any cell      │
│  can be reached from any other through shared points.               │
│                                                                     │
│  The vtkConnectivityFilter identifies and extracts these regions.   │
│                                                                     │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  EXTRACTION MODES:                                                  │
│                                                                     │
│  1. AllRegions       - Extract all regions, color each uniquely     │
│  2. LargestRegion    - Extract only the largest region              │
│  3. ClosestPoint     - Extract region nearest to a point            │
│  4. SpecifiedRegions - Extract specific region(s) by ID             │
│  5. PointSeededRegions - Extract region containing seed point(s)    │
│  6. CellSeededRegions  - Extract region containing seed cell(s)     │
│                                                                     │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  REAL-WORLD APPLICATIONS:                                           │
│                                                                     │
│  • Medical: Isolate organs/tumors from CT/MRI scans                 │
│  • CAD: Identify separate parts in an assembly                      │
│  • Mesh Processing: Remove small disconnected fragments             │
│  • Scientific Viz: Extract specific structures of interest          │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

In this demonstration:
- 6 separate objects (spheres, cube, cylinder, cone) are combined
- Use the combo box to select different extraction modes
- Watch how each mode extracts different portions of the data
""")


def main():
    """Launch the VTK Connectivity Filter Demo application."""
    print_educational_summary()

    print("\nStarting interactive demonstration...")
    print("Use the combo box to select different extraction modes.")
    print("-" * 70)

    app = QApplication(sys.argv)
    window = ConnectivityDemo()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
