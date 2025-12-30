"""
VTK PolyData: The Versatile Surface Representation for CFD

This module provides an educational demonstration of vtkPolyData, which is one
of the most important data structures in VTK for representing and manipulating
surface geometry in Computational Fluid Dynamics (CFD) and visualization.

What is PolyData?
-----------------
vtkPolyData is a concrete dataset that represents a geometric structure
consisting of vertices, lines, polygons, and triangle strips. It's the
primary format for surface meshes in VTK.

PolyData can contain:
- Vertices (verts): 0D point cells
- Lines: 1D edge cells
- Polygons (polys): 2D surface cells (triangles, quads, etc.)
- Triangle strips: Efficient representation of connected triangles

PolyData Components:
--------------------
1. POINTS (vtkPoints):
   - The geometric coordinates of all vertices
   - Shared by all cell types in the polydata

2. CELLS (vtkCellArray):
   - Define connectivity between points
   - Stored separately for verts, lines, polys, and strips

3. POINT DATA:
   - Scalar/vector/tensor fields at each vertex
   - Used for coloring, contouring, glyph sizing

4. CELL DATA:
   - Fields defined per cell (not per point)
   - Useful for cell-centered quantities

CFD/FEA Applications:
--------------------
PolyData is essential for CFD workflows:

1. SURFACE MESHES:
   - Vehicle body surfaces for aerodynamics
   - Wing surfaces for aerospace CFD
   - Pipe walls for internal flow

2. BOUNDARY REPRESENTATION:
   - Wall boundaries for no-slip conditions
   - Inlet/outlet surfaces for boundary conditions
   - Symmetry planes

3. RESULT VISUALIZATION:
   - Isosurfaces of pressure, temperature, velocity magnitude
   - Streamline ribbons and tubes
   - Cutting planes through 3D data

4. GEOMETRY IMPORT/EXPORT:
   - STL files (3D printing, CAD exchange)
   - OBJ files (graphics, visualization)
   - PLY files (point cloud data)

This example demonstrates creating, manipulating, and visualizing polydata
representing various CFD-relevant geometries with associated field data.
"""

import math

import vtk


# Colors for different geometry types
GEOMETRY_COLORS = {
    "triangle": (1.0, 0.3, 0.3),     # Light red
    "quad": (0.3, 1.0, 0.3),         # Light green
    "polygon": (0.3, 0.3, 1.0),      # Light blue
    "airfoil": (1.0, 0.8, 0.2),      # Gold
    "cylinder": (0.8, 0.3, 0.8),     # Purple
}


def create_triangle_polydata():
    """
    Create a simple triangle as polydata.

    Triangles are the most fundamental polygon in CFD surface meshes:
    - Always planar (3 points define a plane)
    - Can approximate any surface geometry
    - Used in STL files, the standard for CAD/CFD exchange
    - Robust for intersection and Boolean operations

    Returns:
        vtkPolyData: Triangle geometry with normal vector
    """
    points = vtk.vtkPoints()
    # Equilateral triangle with unit side length
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(0.5, 0.866, 0.0)

    # Create triangle cell
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 0)
    triangle.GetPointIds().SetId(1, 1)
    triangle.GetPointIds().SetId(2, 2)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(triangle)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    # Compute normals (important for lighting and visualization)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.ComputePointNormalsOn()
    normals.Update()

    return normals.GetOutput()


def create_quad_polydata():
    """
    Create a quadrilateral (quad) as polydata.

    Quads are preferred in many CFD applications:
    - More efficient than triangles (fewer cells for same area)
    - Better aligned with structured flow patterns
    - Required for some mesh generation tools
    - Provide higher-order interpolation accuracy

    Returns:
        vtkPolyData: Quad geometry
    """
    points = vtk.vtkPoints()
    # Unit square
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 0.0, 0.0)
    points.InsertNextPoint(1.0, 1.0, 0.0)
    points.InsertNextPoint(0.0, 1.0, 0.0)

    # Create quad cell (vertices in counter-clockwise order)
    quad = vtk.vtkQuad()
    for i in range(4):
        quad.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(quad)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    return polydata


def create_airfoil_polydata(chord=1.0, num_points=50, thickness=0.12):
    """
    Create a NACA 4-digit airfoil profile as polydata.

    Airfoils are fundamental to aerospace CFD:
    - Wings, turbine blades, propellers
    - NACA profiles: standardized shapes with known properties
    - Used in 2D section analysis and 3D wing construction

    The thickness distribution follows the NACA 00xx formula:
    y = 5t * c * (0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
    where t is the maximum thickness ratio and c is the chord length.

    Args:
        chord: Airfoil chord length
        num_points: Number of points along each surface
        thickness: Maximum thickness as fraction of chord (e.g., 0.12 = 12%)

    Returns:
        vtkPolyData: Airfoil profile (closed polygon)
    """
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    # Generate upper and lower surface points
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(2 * num_points)

    point_id = 0

    # Upper surface (trailing edge to leading edge)
    for i in range(num_points):
        x = chord * (1 - i / (num_points - 1))  # From TE to LE
        x_norm = x / chord

        # NACA thickness distribution
        yt = thickness / 0.2 * chord * (
            0.2969 * math.sqrt(x_norm) -
            0.1260 * x_norm -
            0.3516 * x_norm**2 +
            0.2843 * x_norm**3 -
            0.1015 * x_norm**4
        )

        points.InsertNextPoint(x, yt, 0)
        polygon.GetPointIds().SetId(point_id, point_id)
        point_id += 1

    # Lower surface (leading edge to trailing edge)
    for i in range(num_points):
        x = chord * (i / (num_points - 1))  # From LE to TE
        x_norm = x / chord

        # NACA thickness distribution (negative for lower surface)
        yt = -thickness / 0.2 * chord * (
            0.2969 * math.sqrt(max(x_norm, 1e-10)) -
            0.1260 * x_norm -
            0.3516 * x_norm**2 +
            0.2843 * x_norm**3 -
            0.1015 * x_norm**4
        )

        points.InsertNextPoint(x, yt, 0)
        polygon.GetPointIds().SetId(point_id, point_id)
        point_id += 1

    cells.InsertNextCell(polygon)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

    return polydata


def create_cylinder_surface(radius=0.5, height=2.0, num_radial=20, num_axial=10):
    """
    Create a cylinder surface as polydata.

    Cylinders are common in CFD:
    - Pipe flow simulations
    - Heat exchangers (tube bundles)
    - Vortex shedding studies (flow past cylinder)
    - Pressure vessel analysis

    Args:
        radius: Cylinder radius
        height: Cylinder height
        num_radial: Points around circumference
        num_axial: Points along axis

    Returns:
        vtkPolyData: Cylinder surface mesh
    """
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(height)
    cylinder.SetResolution(num_radial)
    cylinder.CappingOn()
    cylinder.Update()

    return cylinder.GetOutput()


def create_sphere_surface(radius=0.5, theta_resolution=20, phi_resolution=20):
    """
    Create a sphere surface as polydata.

    Spheres are used in CFD for:
    - Flow around bluff bodies (drag coefficient studies)
    - Particle simulations (droplets, bubbles)
    - Validation cases (Stokes flow, potential flow)

    Args:
        radius: Sphere radius
        theta_resolution: Divisions around pole axis
        phi_resolution: Divisions from pole to pole

    Returns:
        vtkPolyData: Sphere surface mesh
    """
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(theta_resolution)
    sphere.SetPhiResolution(phi_resolution)
    sphere.Update()

    return sphere.GetOutput()


def add_scalar_field(polydata, field_name="Pressure", mode="distance"):
    """
    Add a scalar field to polydata points.

    Scalar fields on surfaces are fundamental for CFD visualization:
    - Pressure distribution on airfoils, vehicles
    - Temperature on heat transfer surfaces
    - Wall shear stress for boundary layer analysis
    - y+ values for turbulence modeling

    Args:
        polydata: vtkPolyData to modify
        field_name: Name of the scalar field
        mode: How to compute values ("distance", "height", "random")

    Returns:
        vtkPolyData: Modified polydata with scalar field
    """
    scalars = vtk.vtkFloatArray()
    scalars.SetName(field_name)
    scalars.SetNumberOfComponents(1)

    n_points = polydata.GetNumberOfPoints()

    for i in range(n_points):
        x, y, z = polydata.GetPoint(i)

        if mode == "distance":
            # Distance from origin (good for sphere/cylinder)
            value = math.sqrt(x**2 + y**2 + z**2)
        elif mode == "height":
            # Y-coordinate (good for airfoils, vertical surfaces)
            value = y
        elif mode == "x_position":
            # X-coordinate (good for streamwise variation)
            value = x
        else:
            # Random noise
            value = math.sin(x * 10) * math.cos(y * 10)

        scalars.InsertNextValue(value)

    polydata.GetPointData().AddArray(scalars)
    polydata.GetPointData().SetActiveScalars(field_name)

    return polydata


def add_vector_field(polydata, field_name="SurfaceNormal"):
    """
    Add a vector field (surface normals) to polydata.

    Surface normals are essential in CFD:
    - Define wall boundary orientation
    - Calculate surface forces (pressure * normal)
    - Determine flux directions for boundary conditions
    - Used in rendering for proper lighting

    Args:
        polydata: vtkPolyData to modify
        field_name: Name of the vector field

    Returns:
        vtkPolyData: Modified polydata with normals
    """
    normals_filter = vtk.vtkPolyDataNormals()
    normals_filter.SetInputData(polydata)
    normals_filter.ComputePointNormalsOn()
    normals_filter.ComputeCellNormalsOn()
    normals_filter.SplittingOff()
    normals_filter.Update()

    return normals_filter.GetOutput()


def compute_surface_properties(polydata):
    """
    Compute and print surface properties.

    Surface properties are important for CFD validation:
    - Surface area for force coefficient calculations
    - Number of cells for mesh resolution assessment
    - Bounds for domain sizing

    Args:
        polydata: vtkPolyData to analyze
    """
    # Compute mass properties (area, volume for closed surfaces)
    mass_props = vtk.vtkMassProperties()
    mass_props.SetInputData(polydata)

    # Get basic properties
    n_points = polydata.GetNumberOfPoints()
    n_cells = polydata.GetNumberOfCells()
    bounds = polydata.GetBounds()

    print(f"\n{'='*60}")
    print("Surface Properties")
    print('='*60)
    print(f"  Points: {n_points}")
    print(f"  Cells: {n_cells}")
    print(f"  Bounds:")
    print(f"    X: [{bounds[0]:.4f}, {bounds[1]:.4f}]")
    print(f"    Y: [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"    Z: [{bounds[4]:.4f}, {bounds[5]:.4f}]")

    try:
        surface_area = mass_props.GetSurfaceArea()
        print(f"  Surface Area: {surface_area:.4f}")
    except Exception:
        print("  Surface Area: (computation requires closed surface)")


def visualize_polydata_collection(polydata_list, names=None, show_edges=True):
    """
    Visualize multiple polydata objects with proper CFD visualization.

    This demonstrates CFD-style surface visualization:
    - Color mapping by scalar field
    - Edge visibility for mesh inspection
    - Multiple surfaces with legend
    - Interactive rotation/zoom

    Args:
        polydata_list: List of vtkPolyData objects
        names: List of names for legend
        show_edges: Whether to show mesh edges
    """
    if names is None:
        names = [f"Surface {i}" for i in range(len(polydata_list))]

    colors = list(GEOMETRY_COLORS.values())

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Offset for positioning multiple objects
    x_offset = 0

    for i, polydata in enumerate(polydata_list):
        color = colors[i % len(colors)]

        # Check if polydata has scalars for color mapping
        has_scalars = polydata.GetPointData().GetScalars() is not None

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)

        if has_scalars:
            scalar_range = polydata.GetPointData().GetScalars().GetRange()
            mapper.SetScalarRange(scalar_range)
            mapper.SetScalarModeToUsePointData()

            # Create color lookup table
            lut = vtk.vtkLookupTable()
            lut.SetHueRange(0.667, 0.0)  # Blue to red
            lut.Build()
            mapper.SetLookupTable(lut)
        else:
            mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        if not has_scalars:
            actor.GetProperty().SetColor(color)

        if show_edges:
            actor.GetProperty().SetEdgeVisibility(1)
            actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
            actor.GetProperty().SetLineWidth(1)

        # Position objects side by side
        bounds = polydata.GetBounds()
        width = bounds[1] - bounds[0]
        actor.SetPosition(x_offset, 0, 0)
        x_offset += width + 0.5

        renderer.AddActor(actor)

    # Create legend
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(len(polydata_list))
    legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPositionCoordinate().SetValue(0.02, 0.75)
    legend.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPosition2Coordinate().SetValue(0.2, 0.95)
    legend.UseBackgroundOn()
    legend.SetBackgroundColor(0.1, 0.1, 0.1)

    sphere = vtk.vtkSphereSource()
    sphere.Update()

    for i, (name, color) in enumerate(zip(names, colors)):
        legend.SetEntry(i, sphere.GetOutput(), name, color)

    renderer.AddActor(legend)

    # Axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK PolyData: CFD Surface Meshes")

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
    Print educational summary about VTK PolyData in CFD context.
    """
    print("\n" + "=" * 70)
    print("VTK PolyData: Educational Summary for CFD/FEA")
    print("=" * 70)
    print("""
1. WHAT IS POLYDATA?
   - VTK's primary format for surface geometry
   - Contains points, cells, and associated data
   - Supports vertices, lines, polygons, and triangle strips

2. POLYDATA STRUCTURE:
   - Points: Vertex coordinates (shared by all cells)
   - Verts: Point cells (0D)
   - Lines: Edge cells (1D)
   - Polys: Surface cells (2D triangles, quads, polygons)
   - Strips: Efficient triangle strip representation

3. ASSOCIATED DATA:
   - Point Data: Fields defined at vertices (interpolated)
   - Cell Data: Fields defined per cell (constant per cell)
   - Field Data: Global data (metadata, parameters)

4. CFD APPLICATIONS:
   - Surface meshes for boundary representation
   - Wall boundaries for CFD simulations
   - Isosurfaces from volumetric data
   - Streamline tubes and ribbons
   - Cutting planes through 3D domains

5. FILE FORMATS:
   - VTK/VTP: Native VTK format
   - STL: Standard for CAD/CFD exchange
   - OBJ: Wavefront format with materials
   - PLY: Point cloud with colors

6. COMMON OPERATIONS:
   - Triangulation: Convert to triangles only
   - Decimation: Reduce cell count
   - Smoothing: Reduce surface noise
   - Normals: Compute surface normals
   - Subdivision: Increase mesh resolution
""")


def main():
    """
    Main function demonstrating VTK PolyData for CFD applications.

    Demonstrates:
    - Creating various polydata geometries
    - Adding scalar fields for visualization
    - Computing surface properties
    - Multi-object visualization
    """
    print_educational_summary()

    # Create different polydata objects
    print("\n" + "-" * 70)
    print("Creating CFD-relevant surface geometries...")
    print("-" * 70)

    # Simple shapes
    triangle = create_triangle_polydata()
    compute_surface_properties(triangle)
    print("  Created: Triangle (basic CFD mesh element)")

    quad = create_quad_polydata()
    compute_surface_properties(quad)
    print("  Created: Quadrilateral (structured mesh element)")

    # Airfoil profile
    airfoil = create_airfoil_polydata(chord=1.0, num_points=40, thickness=0.12)
    add_scalar_field(airfoil, "Pressure", mode="x_position")
    compute_surface_properties(airfoil)
    print("  Created: NACA 0012 Airfoil (aerospace CFD)")

    # 3D shapes
    cylinder = create_cylinder_surface(radius=0.3, height=1.0)
    add_scalar_field(cylinder, "Temperature", mode="height")
    compute_surface_properties(cylinder)
    print("  Created: Cylinder (pipe flow, heat exchanger)")

    sphere = create_sphere_surface(radius=0.3)
    add_scalar_field(sphere, "Pressure", mode="distance")
    compute_surface_properties(sphere)
    print("  Created: Sphere (bluff body, particle)")

    # Visualize all surfaces
    visualize_polydata_collection(
        [triangle, quad, airfoil, cylinder, sphere],
        names=["Triangle", "Quad", "Airfoil", "Cylinder", "Sphere"],
        show_edges=True
    )


if __name__ == "__main__":
    main()
