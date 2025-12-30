"""
VTK Points: The Foundation of Computational Geometry

This module provides an educational demonstration of VTK points, which are the
fundamental building blocks for all geometric representations in Computational
Fluid Dynamics (CFD) and Finite Element Analysis (FEA).

What are Points?
----------------
In VTK (and numerical simulations), a point is a location in 3D space defined
by (x, y, z) coordinates. Points are the vertices that define the geometry of
all graphical primitives - lines, polygons, and volumes.

The vtkPoints class stores an array of 3D points and provides methods for:
- Inserting and retrieving points
- Coordinate transformations
- Memory-efficient storage

Points vs Nodes:
----------------
In CFD/FEA terminology:
- "Points" (VTK term) = "Nodes" (FEA term) = "Vertices" (geometry term)
- These are locations where field values (velocity, pressure, temperature) are
  defined or interpolated
- The quality of point placement affects simulation accuracy

CFD/FEA Context:
----------------
Points are crucial in computational simulations:

1. MESH GENERATION:
   - Points define the vertices of mesh cells
   - Point distribution affects mesh quality and solution accuracy
   - Clustering points near walls captures boundary layers
   - Uniform spacing in bulk flow regions

2. SOLUTION STORAGE:
   - Field values (velocity, pressure, temperature) stored at points
   - Interpolation between points for visualization
   - Gradient calculations use point coordinates

3. POINT CLOUDS:
   - Experimental data (PIV, LiDAR) often comes as point clouds
   - VTK can import, process, and visualize point cloud data
   - Surface reconstruction from scattered points

4. COORDINATE SYSTEMS:
   - Cartesian (x, y, z) - most common
   - Cylindrical (r, θ, z) - pipe flows, rotating machinery
   - Spherical (r, θ, φ) - atmospheric simulations

This example demonstrates creating, manipulating, and converting VTK points,
including interoperability with NumPy for scientific computing workflows.
"""

import math

import numpy as np
import vtk
import vtk.util.numpy_support as vtk_np


# Default point size for visualization
DEFAULT_POINT_SIZE = 15

# Colors for different point sets
POINT_COLORS = {
    "original": (1.0, 0.0, 0.0),    # Red
    "transformed": (0.0, 1.0, 0.0),  # Green
    "numpy": (0.0, 0.0, 1.0),        # Blue
}


def create_vtk_points_simple():
    """
    Create a basic vtkPoints object with manually inserted points.

    This demonstrates the fundamental way to create points in VTK:
    - Points are added one at a time using InsertNextPoint()
    - Each point is automatically assigned an index (0, 1, 2, ...)
    - These indices are used by cells to reference points

    In CFD, this approach is used for:
    - Defining specific probe locations
    - Creating custom geometries programmatically
    - Setting up boundary condition points

    Returns:
        vtkPoints: Object containing three points forming a triangle
    """
    points = vtk.vtkPoints()
    # Define three points forming an equilateral triangle in the XY plane
    points.InsertNextPoint(0.0, 0.0, 0.0)       # Point 0: Origin
    points.InsertNextPoint(1.0, 0.0, 0.0)       # Point 1: X-axis
    points.InsertNextPoint(0.5, 0.866, 0.0)    # Point 2: Apex (height = sqrt(3)/2)
    return points


def create_structured_point_grid(nx, ny, nz, spacing=1.0):
    """
    Create a structured grid of points (common in CFD preprocessing).

    Structured grids are the foundation of many CFD simulations:
    - Regular spacing simplifies finite difference calculations
    - Easy to implement high-order numerical schemes
    - Efficient memory access patterns

    Args:
        nx, ny, nz: Number of points in each direction
        spacing: Distance between adjacent points

    Returns:
        vtkPoints: Structured grid of points
    """
    points = vtk.vtkPoints()

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = i * spacing
                y = j * spacing
                z = k * spacing
                points.InsertNextPoint(x, y, z)

    return points


def create_boundary_layer_points(wall_y=0.0, first_cell_height=0.01,
                                  growth_rate=1.2, num_layers=10, x_points=5,
                                  streamwise_spacing=0.5):
    """
    Create points with boundary layer clustering near a wall.

    Boundary layer meshes are critical in CFD because:
    - Velocity gradients are highest near walls
    - Accurate wall shear stress requires fine mesh at walls
    - y+ (dimensionless wall distance) determines mesh requirements

    The growth rate determines how quickly cells expand away from the wall:
    - Growth rate 1.0: Uniform spacing
    - Growth rate 1.1-1.2: Typical for well-resolved boundary layers
    - Growth rate > 1.5: May cause numerical errors

    Args:
        wall_y: Y-coordinate of the wall
        first_cell_height: Height of first cell (determines y+)
        growth_rate: Expansion ratio between successive layers
        num_layers: Number of boundary layer cells
        x_points: Number of points in streamwise direction
        streamwise_spacing: Distance between points in x-direction

    Returns:
        vtkPoints: Points clustered near the wall
    """
    points = vtk.vtkPoints()

    # Calculate y-coordinates using geometric progression
    y_coords = [wall_y]
    current_height = first_cell_height
    current_y = wall_y

    for _ in range(num_layers):
        current_y += current_height
        y_coords.append(current_y)
        current_height *= growth_rate

    # Create 2D grid of points
    for y in y_coords:
        for i in range(x_points):
            x = i * streamwise_spacing
            points.InsertNextPoint(x, y, 0.0)

    return points


def create_cylindrical_points(radius=1.0, num_radial=8, num_axial=5, length=2.0):
    """
    Create points arranged in a cylindrical pattern.

    Cylindrical coordinates are natural for:
    - Pipe flow simulations
    - Rotating machinery (turbines, pumps)
    - Axisymmetric geometries

    Args:
        radius: Cylinder radius
        num_radial: Number of points around circumference
        num_axial: Number of points along axis
        length: Cylinder length

    Returns:
        vtkPoints: Points arranged cylindrically
    """
    points = vtk.vtkPoints()

    for k in range(num_axial):
        z = k * length / max(num_axial - 1, 1)
        for i in range(num_radial):
            theta = 2.0 * math.pi * i / num_radial
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            points.InsertNextPoint(x, y, z)

    return points


def print_points_info(points, name="Points"):
    """
    Print detailed information about a vtkPoints object.

    Understanding point data is essential for debugging CFD meshes:
    - Point count indicates mesh resolution
    - Coordinate ranges show domain extents
    - Individual point inspection for quality checks

    Args:
        points: vtkPoints object to analyze
        name: Descriptive name for output
    """
    n_points = points.GetNumberOfPoints()
    print(f"\n{'='*60}")
    print(f"{name}: {n_points} points")
    print('='*60)

    if n_points == 0:
        print("  No points defined")
        return

    # Calculate bounds
    bounds = [float('inf'), float('-inf'),
              float('inf'), float('-inf'),
              float('inf'), float('-inf')]

    for i in range(n_points):
        x, y, z = points.GetPoint(i)
        bounds[0] = min(bounds[0], x)
        bounds[1] = max(bounds[1], x)
        bounds[2] = min(bounds[2], y)
        bounds[3] = max(bounds[3], y)
        bounds[4] = min(bounds[4], z)
        bounds[5] = max(bounds[5], z)

    print(f"  X range: [{bounds[0]:.4f}, {bounds[1]:.4f}]")
    print(f"  Y range: [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"  Z range: [{bounds[4]:.4f}, {bounds[5]:.4f}]")

    # Print first few points as examples
    print(f"\n  First {min(5, n_points)} points:")
    for i in range(min(5, n_points)):
        x, y, z = points.GetPoint(i)
        print(f"    Point {i}: ({x:.4f}, {y:.4f}, {z:.4f})")

    if n_points > 5:
        print(f"    ... ({n_points - 5} more points)")


def vtk_points_to_numpy(points):
    """
    Convert vtkPoints to a NumPy array.

    NumPy integration is essential for CFD workflows:
    - Apply mathematical operations (scaling, rotation)
    - Interface with scientific Python libraries (SciPy, scikit-learn)
    - Perform statistical analysis on point distributions
    - Export to other formats

    Args:
        points: vtkPoints object

    Returns:
        numpy.ndarray: Shape (N, 3) array of coordinates
    """
    return vtk_np.vtk_to_numpy(points.GetData())


def numpy_to_vtk_points(numpy_array):
    """
    Convert a NumPy array to vtkPoints.

    This enables importing data from:
    - Experimental measurements (PIV, pressure taps)
    - Other simulation tools (OpenFOAM, ANSYS)
    - Mathematical point generation (parametric surfaces)

    Args:
        numpy_array: Shape (N, 3) array of coordinates

    Returns:
        vtkPoints: VTK points object
    """
    vtk_points = vtk.vtkPoints()
    vtk_data = vtk_np.numpy_to_vtk(numpy_array.astype(np.float64))
    vtk_points.SetData(vtk_data)
    return vtk_points


def transform_points(points, translation=(0, 0, 0), scale=1.0, rotation_deg=0.0):
    """
    Apply geometric transformation to points.

    Transformations are common in CFD preprocessing:
    - Scaling: Convert between unit systems (mm to m)
    - Translation: Position geometries in the domain
    - Rotation: Align with flow direction or coordinate system

    Args:
        points: vtkPoints to transform
        translation: (dx, dy, dz) translation vector
        scale: Uniform scaling factor
        rotation_deg: Rotation angle around Z-axis (degrees)

    Returns:
        vtkPoints: Transformed points
    """
    # Convert to numpy for easy manipulation
    np_points = vtk_points_to_numpy(points)

    # Apply rotation around Z-axis
    angle_rad = math.radians(rotation_deg)
    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
    rotation_matrix = np.array([
        [cos_a, -sin_a, 0],
        [sin_a, cos_a, 0],
        [0, 0, 1]
    ])
    np_points = np_points @ rotation_matrix.T

    # Apply scaling
    np_points *= scale

    # Apply translation
    np_points += np.array(translation)

    return numpy_to_vtk_points(np_points)


def visualize_points(points_list, names=None, colors=None):
    """
    Visualize multiple point sets with different colors.

    Point cloud visualization is used in CFD for:
    - Mesh quality inspection
    - Experimental data overlay
    - Solution sampling locations
    - Particle tracking results

    Args:
        points_list: List of vtkPoints objects
        names: List of names for each point set
        colors: List of (r, g, b) color tuples
    """
    if names is None:
        names = [f"Points {i}" for i in range(len(points_list))]
    if colors is None:
        colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1)]

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.15, 0.15, 0.2)

    for i, points in enumerate(points_list):
        color = colors[i % len(colors)]

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        # Use vertex glyph filter to render points
        glyph = vtk.vtkVertexGlyphFilter()
        glyph.SetInputData(polydata)
        glyph.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetPointSize(DEFAULT_POINT_SIZE)

        renderer.AddActor(actor)

    # Create legend
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(len(points_list))
    legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPositionCoordinate().SetValue(0.02, 0.85)
    legend.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPosition2Coordinate().SetValue(0.25, 0.98)
    legend.UseBackgroundOn()
    legend.SetBackgroundColor(0.1, 0.1, 0.1)

    sphere = vtk.vtkSphereSource()
    sphere.Update()

    for i, (name, color) in enumerate(zip(names, colors)):
        legend.SetEntry(i, sphere.GetOutput(), name, color)

    renderer.AddActor(legend)

    # Add axes
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Points: Foundation of CFD Meshes")

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    renderer.ResetCamera()
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


def print_educational_summary():
    """
    Print educational summary about VTK points in CFD context.
    """
    print("\n" + "=" * 70)
    print("VTK Points: Educational Summary for CFD/FEA")
    print("=" * 70)
    print("""
1. WHAT ARE POINTS?
   - Locations in 3D space defined by (x, y, z) coordinates
   - Foundation for all geometric representations in VTK
   - Equivalent to "nodes" in FEA or "vertices" in geometry

2. POINT STORAGE:
   - vtkPoints stores coordinates as contiguous array
   - Efficient memory layout for large meshes
   - Supports float or double precision

3. POINT OPERATIONS:
   - InsertNextPoint(): Add points sequentially
   - GetPoint(i): Retrieve coordinates by index
   - SetPoint(i, x, y, z): Modify existing point
   - GetNumberOfPoints(): Count total points

4. CFD APPLICATIONS:
   - Mesh vertices for finite volume/element methods
   - Probe locations for solution sampling
   - Particle positions in Lagrangian tracking
   - Experimental data import (PIV, LiDAR)

5. NUMPY INTEROPERABILITY:
   - vtk_to_numpy(): Convert VTK to NumPy
   - numpy_to_vtk(): Convert NumPy to VTK
   - Enables scientific Python workflows

6. BEST PRACTICES:
   - Use structured grids when geometry allows
   - Cluster points near walls for boundary layers
   - Check point distribution quality before simulation
   - Use appropriate coordinate system for geometry
""")


def main():
    """
    Main function demonstrating VTK points for CFD applications.

    Demonstrates:
    - Creating points manually and programmatically
    - Structured grid generation
    - Boundary layer point clustering
    - NumPy conversion and transformation
    - Multi-set visualization
    """
    print_educational_summary()

    # Create different point sets
    print("\n" + "-" * 70)
    print("Creating various point configurations...")
    print("-" * 70)

    # Simple manual points
    simple_points = create_vtk_points_simple()
    print_points_info(simple_points, "Simple Triangle Points")

    # Structured grid
    grid_points = create_structured_point_grid(4, 3, 2, spacing=0.5)
    print_points_info(grid_points, "Structured Grid Points")

    # Boundary layer points
    bl_points = create_boundary_layer_points(
        wall_y=0.0, first_cell_height=0.02,
        growth_rate=1.3, num_layers=8, x_points=6
    )
    print_points_info(bl_points, "Boundary Layer Points")

    # Cylindrical points
    cyl_points = create_cylindrical_points(radius=0.5, num_radial=12, num_axial=4)
    print_points_info(cyl_points, "Cylindrical Points")

    # NumPy conversion demonstration
    print("\n" + "-" * 70)
    print("NumPy Interoperability:")
    print("-" * 70)

    np_array = vtk_points_to_numpy(simple_points)
    print(f"  Converted to NumPy array shape: {np_array.shape}")
    print(f"  Array contents:\n{np_array}")

    # Create new points from NumPy
    new_np_points = np.array([
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [2.5, 0.866, 0.0]
    ])
    from_numpy = numpy_to_vtk_points(new_np_points)
    print_points_info(from_numpy, "Points from NumPy")

    # Transform points
    transformed = transform_points(simple_points, translation=(0, 0, 1),
                                   scale=1.5, rotation_deg=45)
    print_points_info(transformed, "Transformed Points")

    # Visualize multiple point sets
    visualize_points(
        [simple_points, from_numpy, transformed],
        names=["Original", "From NumPy", "Transformed"],
        colors=[(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    )


if __name__ == "__main__":
    main()
