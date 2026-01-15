"""
NACA 0012 Airfoil Generation and Visualization Demo

This script demonstrates generating and visualizing the NACA 0012 airfoil,
a symmetric airfoil profile commonly used in aerodynamic research.

For comprehensive NACA airfoil generation with support for:
- NACA 4-digit airfoils (symmetric and cambered)
- NACA 5-digit airfoils
- 3D wing generation (tapered, swept, twisted)
- Multiple export formats (VTP, STL, OBJ, CSV, DAT)
- Comparison visualizations

See the naca_airfoil_generator.py module in this directory.

Usage:
    python generate_naca0012.py

    # Or use the comprehensive generator:
    python naca_airfoil_generator.py 0012 --visualize
    python naca_airfoil_generator.py 0012 --export output.stl
    python naca_airfoil_generator.py 0012 --wing --span 10
"""

import vtk

from naca_airfoil_generator import NACAGenerator, WingParameters


def main():
    """Generate and visualize NACA 0012 airfoil using the comprehensive generator."""

    print("=" * 70)
    print("NACA 0012 Airfoil Generation Demo")
    print("=" * 70)

    # Create NACA 0012 generator
    generator = NACAGenerator(
        naca_code="0012",
        chord=1.0,
        num_points=100,
        cosine_spacing=True,  # Better resolution at LE and TE
    )

    # Print airfoil information
    print(generator)
    print()

    # Create 2D polydata
    polydata_2d = generator.create_polydata_2d()

    # Export to VTP file
    output_file = "naca0012.vtp"
    generator.export(output_file)
    print(f"Exported 2D profile to: {output_file}")

    # Also export to other formats for demonstration
    generator.export("naca0012.csv")
    print("Exported coordinates to: naca0012.csv")

    generator.export("naca0012.dat")
    print("Exported XFOIL format to: naca0012.dat")

    print()
    print("-" * 70)
    print("Visualizing 2D airfoil profile...")
    print("-" * 70)

    # Visualize 2D profile
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata_2d)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.2, 0.4, 0.8)  # Blue color
    actor.GetProperty().SetEdgeVisibility(1)
    actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
    actor.GetProperty().SetLineWidth(2)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1, 1, 1)
    renderer.AddActor(actor)

    # Add text annotation
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("NACA 0012 Airfoil\nSymmetric, 12% thickness")
    text_actor.GetTextProperty().SetFontSize(18)
    text_actor.GetTextProperty().SetColor(0.1, 0.1, 0.1)
    text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    text_actor.GetPositionCoordinate().SetValue(0.02, 0.92)
    renderer.AddActor2D(text_actor)

    # Add axes
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1000, 500)
    render_window.SetWindowName("NACA 0012 Airfoil - 2D Profile")

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Set camera for 2D view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0.5, 0, 2)
    camera.SetFocalPoint(0.5, 0, 0)
    camera.SetViewUp(0, 1, 0)
    renderer.ResetCamera()

    interactor.Initialize()
    render_window.Render()
    interactor.Start()


if __name__ == "__main__":
    main()
