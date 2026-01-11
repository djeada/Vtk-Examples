"""
Topology vs Geometry: Understanding the Fundamental Distinction in VTK

This module provides an interactive demonstration of the difference between
topology and geometry in VTK and computational mesh representations.

What is Topology?
-----------------
Topology describes the CONNECTIVITY and STRUCTURE of a mesh:
- Which points form a cell
- Which cells are neighbors
- How elements are logically connected
- The "skeleton" of relationships between elements

Topology answers: "Who is connected to whom?"

Topology does NOT care about:
- Distances between points
- Angles between edges
- Shape or size of cells
- Position in physical space

What is Geometry?
-----------------
Geometry describes the SPATIAL POSITION of points:
- Coordinates (x, y, z) of each vertex
- Cell shape, size, and orientation
- Curvature and deformation
- Physical location in space

Geometry answers: "Where are things located?"

Geometry does NOT care about:
- Connectivity rules
- Which points form which cells
- Neighbor relationships

Key Insight:
------------
You can:
- Change GEOMETRY without changing TOPOLOGY (warp, stretch, twist a mesh)
- Change TOPOLOGY without changing GEOMETRY (re-mesh the same shape)

This example demonstrates a hexahedral cell (8 vertices, 12 edges, 6 faces)
where the TOPOLOGY remains constant while the GEOMETRY is interactively
modified through various transformations:
- Scaling (stretch/compress)
- Shearing (skew)
- Twisting (torsion)
- Bending (curvature)

All transformations preserve the topological relationships:
- Still 8 vertices
- Still 12 edges
- Still 6 faces
- Same connectivity between vertices

But the GEOMETRY changes dramatically, creating different physical shapes
from the same topological structure.
"""

import numpy as np
import vtk
from vtkmodules.vtkInteractionWidgets import vtkSliderWidget, vtkSliderRepresentation2D


# Default hexahedron vertices (unit cube centered at origin)
DEFAULT_VERTICES = np.array([
    [-0.5, -0.5, -0.5],  # 0: bottom-front-left
    [ 0.5, -0.5, -0.5],  # 1: bottom-front-right
    [ 0.5,  0.5, -0.5],  # 2: bottom-back-right
    [-0.5,  0.5, -0.5],  # 3: bottom-back-left
    [-0.5, -0.5,  0.5],  # 4: top-front-left
    [ 0.5, -0.5,  0.5],  # 5: top-front-right
    [ 0.5,  0.5,  0.5],  # 6: top-back-right
    [-0.5,  0.5,  0.5],  # 7: top-back-left
], dtype=np.float64)


class TopologyGeometryDemo:
    """
    Interactive demonstration of topology vs geometry concepts.

    This class creates a VTK visualization that allows users to interactively
    modify the geometry of a hexahedral cell while the topology remains constant.
    """

    def __init__(self):
        """Initialize the demo with default values."""
        # Transformation parameters
        self.scale_x = 1.0
        self.scale_y = 1.0
        self.scale_z = 1.0
        self.shear_xy = 0.0
        self.twist_z = 0.0
        self.bend_factor = 0.0

        # Create VTK objects
        self.points = vtk.vtkPoints()
        self.ugrid = vtk.vtkUnstructuredGrid()
        self.setup_topology()

        # Renderer and window setup
        self.renderer = vtk.vtkRenderer()
        self.render_window = vtk.vtkRenderWindow()
        self.interactor = vtk.vtkRenderWindowInteractor()

        # Actors
        self.surface_actor = None
        self.edge_actor = None
        self.point_actor = None

        # Sliders
        self.sliders = []

    def setup_topology(self):
        """
        Set up the hexahedral cell topology.

        This defines the CONNECTIVITY - which points form the cell.
        The topology is set once and never changes, demonstrating
        that we can modify geometry while preserving topology.
        """
        # Initialize points with default geometry
        self.update_geometry()

        # Create hexahedron cell - THIS IS THE TOPOLOGY
        # It defines that points 0-7 form a hexahedron
        hexa = vtk.vtkHexahedron()
        for i in range(8):
            hexa.GetPointIds().SetId(i, i)

        # Add cell to unstructured grid
        self.ugrid.SetPoints(self.points)
        self.ugrid.InsertNextCell(vtk.VTK_HEXAHEDRON, hexa.GetPointIds())

    def apply_transformations(self, vertices):
        """
        Apply geometric transformations to vertices.

        These transformations change the GEOMETRY (positions)
        while the TOPOLOGY (connectivity) remains unchanged.

        Args:
            vertices: numpy array of shape (8, 3) with vertex coordinates

        Returns:
            Transformed vertices as numpy array
        """
        result = vertices.copy()

        # Apply scaling (changes size, not connectivity)
        result[:, 0] *= self.scale_x
        result[:, 1] *= self.scale_y
        result[:, 2] *= self.scale_z

        # Apply shearing (skews the shape, not connectivity)
        result[:, 0] += self.shear_xy * result[:, 1]

        # Apply twist around Z axis (rotates based on height)
        if abs(self.twist_z) > 0.01:
            for i in range(len(result)):
                z = result[i, 2]
                angle = np.radians(self.twist_z * z)
                x, y = result[i, 0], result[i, 1]
                result[i, 0] = x * np.cos(angle) - y * np.sin(angle)
                result[i, 1] = x * np.sin(angle) + y * np.cos(angle)

        # Apply bending (curves the shape based on height)
        if abs(self.bend_factor) > 0.01:
            for i in range(len(result)):
                z = result[i, 2]
                result[i, 0] += self.bend_factor * z * z

        return result

    def update_geometry(self):
        """
        Update the geometry of the hexahedron.

        This modifies the GEOMETRY (point positions) while
        the TOPOLOGY (cell connectivity) remains unchanged.
        """
        # Apply transformations to get new vertex positions
        transformed = self.apply_transformations(DEFAULT_VERTICES)

        # Update VTK points
        self.points.Reset()
        for i, (x, y, z) in enumerate(transformed):
            self.points.InsertNextPoint(x, y, z)

        self.points.Modified()
        if self.ugrid.GetPoints():
            self.ugrid.Modified()

    def create_visualization(self):
        """
        Create the VTK visualization pipeline.

        Sets up actors for:
        - Transparent surface rendering
        - Edge visualization (shows topology)
        - Vertex spheres (shows points)
        """
        # Surface mapper and actor (transparent to see structure)
        surface_mapper = vtk.vtkDataSetMapper()
        surface_mapper.SetInputData(self.ugrid)

        self.surface_actor = vtk.vtkActor()
        self.surface_actor.SetMapper(surface_mapper)
        self.surface_actor.GetProperty().SetColor(0.3, 0.6, 0.9)  # Light blue
        self.surface_actor.GetProperty().SetOpacity(0.3)
        self.surface_actor.GetProperty().SetEdgeVisibility(1)
        self.surface_actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.6)
        self.surface_actor.GetProperty().SetLineWidth(3)

        # Extract edges for clearer topology visualization
        edges = vtk.vtkExtractEdges()
        edges.SetInputData(self.ugrid)

        edge_mapper = vtk.vtkPolyDataMapper()
        edge_mapper.SetInputConnection(edges.GetOutputPort())

        self.edge_actor = vtk.vtkActor()
        self.edge_actor.SetMapper(edge_mapper)
        self.edge_actor.GetProperty().SetColor(0.1, 0.1, 0.4)
        self.edge_actor.GetProperty().SetLineWidth(4)

        # Point visualization using spheres
        sphere_source = vtk.vtkSphereSource()
        sphere_source.SetRadius(0.08)
        sphere_source.SetThetaResolution(16)
        sphere_source.SetPhiResolution(16)

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(self.ugrid)
        glyph.SetSourceConnection(sphere_source.GetOutputPort())
        glyph.SetScaleModeToDataScalingOff()

        point_mapper = vtk.vtkPolyDataMapper()
        point_mapper.SetInputConnection(glyph.GetOutputPort())

        self.point_actor = vtk.vtkActor()
        self.point_actor.SetMapper(point_mapper)
        self.point_actor.GetProperty().SetColor(0.9, 0.2, 0.2)  # Red

        # Add actors to renderer
        self.renderer.AddActor(self.surface_actor)
        self.renderer.AddActor(self.edge_actor)
        self.renderer.AddActor(self.point_actor)

    def create_info_text(self):
        """
        Create text annotations explaining topology vs geometry.
        """
        # Title
        title = vtk.vtkTextActor()
        title.SetInput("Topology vs Geometry Demonstration")
        title.GetTextProperty().SetFontSize(24)
        title.GetTextProperty().SetColor(1.0, 1.0, 1.0)
        title.GetTextProperty().SetBold(True)
        title.SetPosition(10, 560)
        self.renderer.AddActor2D(title)

        # Topology info (constant)
        topology_text = vtk.vtkTextActor()
        topology_text.SetInput(
            "TOPOLOGY (Constant):\n"
            "  • 8 vertices\n"
            "  • 12 edges\n"
            "  • 6 faces\n"
            "  • 1 hexahedral cell\n"
            "  → Connectivity unchanged!"
        )
        topology_text.GetTextProperty().SetFontSize(14)
        topology_text.GetTextProperty().SetColor(0.4, 0.9, 0.4)  # Green
        topology_text.SetPosition(10, 420)
        self.renderer.AddActor2D(topology_text)

        # Geometry info (variable)
        geometry_text = vtk.vtkTextActor()
        geometry_text.SetInput(
            "GEOMETRY (Variable):\n"
            "  • Point positions change\n"
            "  • Shape/size changes\n"
            "  • Angles change\n"
            "  → Use sliders to modify!"
        )
        geometry_text.GetTextProperty().SetFontSize(14)
        geometry_text.GetTextProperty().SetColor(0.9, 0.6, 0.3)  # Orange
        geometry_text.SetPosition(10, 320)
        self.renderer.AddActor2D(geometry_text)

        # Instructions
        instructions = vtk.vtkTextActor()
        instructions.SetInput(
            "Drag sliders to change GEOMETRY while TOPOLOGY stays constant"
        )
        instructions.GetTextProperty().SetFontSize(12)
        instructions.GetTextProperty().SetColor(0.8, 0.8, 0.8)
        instructions.SetPosition(10, 10)
        self.renderer.AddActor2D(instructions)

    def create_slider(self, title, min_val, max_val, initial_val, y_pos, callback):
        """
        Create a slider widget for parameter control.

        Args:
            title: Slider title
            min_val: Minimum value
            max_val: Maximum value
            initial_val: Initial value
            y_pos: Vertical position (normalized 0-1)
            callback: Function to call on value change

        Returns:
            vtkSliderWidget
        """
        slider_rep = vtkSliderRepresentation2D()
        slider_rep.SetMinimumValue(min_val)
        slider_rep.SetMaximumValue(max_val)
        slider_rep.SetValue(initial_val)
        slider_rep.SetTitleText(title)

        # Position slider on right side of window
        slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint1Coordinate().SetValue(0.65, y_pos)
        slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint2Coordinate().SetValue(0.95, y_pos)

        # Style
        slider_rep.GetTitleProperty().SetColor(1, 1, 1)
        slider_rep.GetSliderProperty().SetColor(0.3, 0.6, 0.9)
        slider_rep.GetCapProperty().SetColor(0.2, 0.4, 0.6)
        slider_rep.GetSelectedProperty().SetColor(0.5, 0.8, 1.0)

        slider_widget = vtkSliderWidget()
        slider_widget.SetInteractor(self.interactor)
        slider_widget.SetRepresentation(slider_rep)
        slider_widget.EnabledOn()
        slider_widget.AddObserver("InteractionEvent", callback)

        return slider_widget

    def on_scale_x_changed(self, obj, event):
        """Callback for X scale slider."""
        self.scale_x = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def on_scale_y_changed(self, obj, event):
        """Callback for Y scale slider."""
        self.scale_y = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def on_scale_z_changed(self, obj, event):
        """Callback for Z scale slider."""
        self.scale_z = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def on_shear_changed(self, obj, event):
        """Callback for shear slider."""
        self.shear_xy = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def on_twist_changed(self, obj, event):
        """Callback for twist slider."""
        self.twist_z = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def on_bend_changed(self, obj, event):
        """Callback for bend slider."""
        self.bend_factor = obj.GetRepresentation().GetValue()
        self.update_geometry()
        self.render_window.Render()

    def setup_sliders(self):
        """Set up all parameter sliders."""
        slider_configs = [
            ("Scale X", 0.2, 3.0, 1.0, 0.90, self.on_scale_x_changed),
            ("Scale Y", 0.2, 3.0, 1.0, 0.78, self.on_scale_y_changed),
            ("Scale Z", 0.2, 3.0, 1.0, 0.66, self.on_scale_z_changed),
            ("Shear XY", -1.0, 1.0, 0.0, 0.54, self.on_shear_changed),
            ("Twist (deg)", -90.0, 90.0, 0.0, 0.42, self.on_twist_changed),
            ("Bend", -0.5, 0.5, 0.0, 0.30, self.on_bend_changed),
        ]

        for title, min_val, max_val, initial, y_pos, callback in slider_configs:
            slider = self.create_slider(title, min_val, max_val, initial, y_pos, callback)
            self.sliders.append(slider)

    def setup_axes(self):
        """Add orientation axes widget."""
        axes = vtk.vtkAxesActor()
        axes_widget = vtk.vtkOrientationMarkerWidget()
        axes_widget.SetOrientationMarker(axes)
        axes_widget.SetInteractor(self.interactor)
        axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
        axes_widget.SetEnabled(1)
        axes_widget.InteractiveOff()

    def run(self):
        """Run the interactive demonstration."""
        # Set up renderer
        self.renderer.SetBackground(0.1, 0.1, 0.15)

        # Set up render window
        self.render_window.AddRenderer(self.renderer)
        self.render_window.SetSize(1200, 600)
        self.render_window.SetWindowName("Topology vs Geometry - VTK Demonstration")

        # Set up interactor
        self.interactor.SetRenderWindow(self.render_window)

        # Create visualization components
        self.create_visualization()
        self.create_info_text()
        self.setup_sliders()
        self.setup_axes()

        # Set up camera for good initial view
        camera = self.renderer.GetActiveCamera()
        camera.SetPosition(3, -3, 2)
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 0, 1)
        self.renderer.ResetCamera()

        # Start interaction
        self.interactor.Initialize()
        self.render_window.Render()
        self.interactor.Start()


def print_educational_summary():
    """Print educational information about topology vs geometry."""
    print("\n" + "=" * 70)
    print("TOPOLOGY vs GEOMETRY: Educational Summary")
    print("=" * 70)
    print("""
┌─────────────────────────────────────────────────────────────────────┐
│                    TOPOLOGY vs GEOMETRY                             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  TOPOLOGY                          GEOMETRY                         │
│  ────────                          ────────                         │
│  • Connectivity                    • Position                       │
│  • "Who connects to whom"          • "Where things are"             │
│  • Structure                       • Shape                          │
│  • Relationships                   • Size                           │
│  • Cell definitions                • Coordinates                    │
│                                                                     │
│  Examples:                         Examples:                        │
│  • 8 vertices in a hex cell        • Cube: vertices at ±0.5        │
│  • Edge between points 0-1         • Stretched: x scaled by 2      │
│  • Face from points 0,1,2,3        • Twisted: rotated with height  │
│                                                                     │
│  ─────────────────────────────────────────────────────────────────  │
│                                                                     │
│  KEY INSIGHT:                                                       │
│  Same TOPOLOGY + Different GEOMETRY = Different SHAPES              │
│                                                                     │
│  A cube, a stretched box, and a twisted prism can all have:        │
│  • 8 vertices                                                       │
│  • 12 edges                                                         │
│  • 6 faces                                                          │
│  • Same connectivity pattern                                        │
│                                                                     │
│  But look completely different because their GEOMETRY differs!      │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘

In this demonstration:
- The hexahedral cell's TOPOLOGY is fixed (8 vertices, 12 edges, 6 faces)
- Use the sliders to modify GEOMETRY (positions, shape, orientation)
- Observe how the shape changes while connectivity stays the same!

CFD/FEA Relevance:
- Mesh deformation changes geometry but preserves topology
- Structured grids have implicit topology but variable geometry
- Understanding this distinction helps with mesh quality assessment
""")


def main():
    """Main function to run the topology vs geometry demonstration."""
    print_educational_summary()

    print("\nStarting interactive demonstration...")
    print("Use the sliders on the right to modify geometry.")
    print("Notice how the TOPOLOGY (vertices, edges, faces) remains constant!")
    print("-" * 70)

    demo = TopologyGeometryDemo()
    demo.run()


if __name__ == "__main__":
    main()
