"""
VTK Structured Grid: Regular Mesh Topology for CFD

This module provides an educational demonstration of vtkStructuredGrid, which is
one of the most important data structures for Computational Fluid Dynamics (CFD)
and Finite Difference/Finite Volume methods.

What is a Structured Grid?
--------------------------
A structured grid is a mesh where:
- Points are arranged in a regular i-j-k (or i-j) topology
- Each interior point has the same number of neighbors
- Connectivity is IMPLICIT (determined by i,j,k indices, not stored explicitly)
- Points can be non-uniformly spaced (unlike vtkImageData)
- Grid lines can be curved (unlike vtkRectilinearGrid)

Structured vs Unstructured:
---------------------------
STRUCTURED GRIDS:
- Implicit connectivity (memory efficient)
- Fast neighbor lookup (O(1) using indices)
- Limited to simple topologies (boxes, O-grids, C-grids)
- Higher accuracy for aligned flows
- Simpler numerical schemes

UNSTRUCTURED GRIDS:
- Explicit connectivity (stored cell-by-cell)
- Slower neighbor lookup
- Can mesh any geometry
- More flexible but more complex

CFD Applications:
-----------------
Structured grids are preferred for:

1. CHANNEL AND PIPE FLOWS:
   - Flow is aligned with grid lines
   - High accuracy with fewer cells
   - Easy boundary layer resolution

2. AERODYNAMIC BODIES:
   - O-grids wrap around airfoils/wings
   - C-grids for wake resolution
   - Body-fitted coordinates

3. TURBOMACHINERY:
   - Blade passages with H-O-H topology
   - Periodic boundaries
   - Tip gap resolution

4. HEAT EXCHANGERS:
   - Regular tube arrangements
   - Structured block approach

Grid Quality Considerations:
----------------------------
- Orthogonality: Angle between grid lines (ideal = 90°)
- Aspect ratio: Ratio of cell dimensions (affects numerical diffusion)
- Smoothness: Gradual cell size changes (affects accuracy)
- Skewness: Deviation from ideal cell shape

This example demonstrates creating various structured grid configurations
commonly used in CFD, including uniform grids, stretched grids for boundary
layers, and body-fitted grids for curved geometries.
"""

import math

import vtk


# Grid quality thresholds for CFD
MAX_RECOMMENDED_ASPECT_RATIO = 100.0
MAX_RECOMMENDED_EXPANSION_RATIO = 1.3
MIN_ORTHOGONALITY_ANGLE = 30.0  # degrees


def create_uniform_structured_grid(nx, ny, nz, spacing=1.0):
    """
    Create a uniform structured grid with constant spacing.

    Uniform grids are the simplest structured mesh:
    - Equal spacing in all directions
    - Ideal for problems with uniform gradients
    - Common for DNS (Direct Numerical Simulation)

    In CFD, uniform grids are used when:
    - Resolution requirements are similar everywhere
    - Studying homogeneous turbulence
    - Initial mesh before refinement

    Args:
        nx, ny, nz: Number of points in each direction
        spacing: Uniform distance between points

    Returns:
        vtkStructuredGrid: Uniform structured grid
    """
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(nx, ny, nz)

    points = vtk.vtkPoints()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                x = i * spacing
                y = j * spacing
                z = k * spacing
                points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)
    return grid


def create_stretched_grid(nx, ny, nz, length_x, length_y, length_z,
                          stretch_y=1.2, wall_clustering=True):
    """
    Create a stretched structured grid with boundary layer clustering.

    Stretched grids concentrate cells near boundaries where gradients are high:
    - First cell height controls y+ (wall units)
    - Growth rate determines boundary layer resolution
    - Typically 10-30 cells across boundary layer

    The stretching function uses geometric progression:
    dy[i+1] = dy[i] * stretch_factor

    Args:
        nx, ny, nz: Number of points in each direction
        length_x, length_y, length_z: Domain dimensions
        stretch_y: Stretching ratio in y-direction
        wall_clustering: If True, cluster at both y=0 and y=length_y

    Returns:
        vtkStructuredGrid: Stretched grid with boundary layer clustering
    """
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(nx, ny, nz)

    # Compute stretched y-coordinates
    if wall_clustering and ny > 2:
        # Hyperbolic tangent stretching for two-sided clustering
        y_coords = []
        for j in range(ny):
            eta = j / (ny - 1)  # 0 to 1
            # Hyperbolic tangent stretching concentrates points at both ends
            beta = 1.5  # Stretching parameter
            y_normalized = 0.5 * (1 + math.tanh(beta * (2 * eta - 1)) / math.tanh(beta))
            y_coords.append(y_normalized * length_y)
    else:
        # Geometric stretching from one wall
        y_coords = []
        if ny > 1:
            # Calculate first cell size for geometric progression
            first_dy = length_y * (stretch_y - 1) / (stretch_y ** (ny - 1) - 1)
            y = 0.0
            for j in range(ny):
                y_coords.append(y)
                if j < ny - 1:
                    dy = first_dy * (stretch_y ** j)
                    y += dy
        else:
            y_coords = [0.0]

    points = vtk.vtkPoints()
    for k in range(nz):
        z = k * length_z / max(nz - 1, 1)
        for j in range(ny):
            y = y_coords[j]
            for i in range(nx):
                x = i * length_x / max(nx - 1, 1)
                points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)
    return grid


def create_ogrid_around_cylinder(radius_inner, radius_outer, nradial, ntheta, nz,
                                  length_z=1.0):
    """
    Create an O-grid topology around a cylinder.

    O-grids are essential for body-fitted CFD meshes:
    - Grid wraps completely around the body
    - No singularities at surface
    - Natural boundary layer resolution
    - Used for airfoils, cylinders, spheres

    The O-grid topology provides:
    - Orthogonal cells at the surface (important for accuracy)
    - Gradual cell growth away from body
    - Smooth grid line distribution

    Args:
        radius_inner: Inner radius (body surface)
        radius_outer: Outer radius (far-field boundary)
        nradial: Number of points in radial direction
        ntheta: Number of points around circumference
        nz: Number of points in axial direction
        length_z: Length in z-direction

    Returns:
        vtkStructuredGrid: O-grid topology around cylinder
    """
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(nradial, ntheta, nz)

    # Geometric stretching for radial direction
    stretch_ratio = 1.15
    radii = []
    dr_first = (radius_outer - radius_inner) * (stretch_ratio - 1) / (stretch_ratio ** (nradial - 1) - 1)
    r = radius_inner
    for i in range(nradial):
        radii.append(r)
        if i < nradial - 1:
            dr = dr_first * (stretch_ratio ** i)
            r += dr

    points = vtk.vtkPoints()
    for k in range(nz):
        z = k * length_z / max(nz - 1, 1)
        for j in range(ntheta):
            theta = 2 * math.pi * j / ntheta
            for i in range(nradial):
                r = radii[i]
                x = r * math.cos(theta)
                y = r * math.sin(theta)
                points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)
    return grid


def create_channel_flow_grid(length, height, width, nx, ny, nz,
                              wall_stretch=1.2):
    """
    Create a structured grid for channel flow CFD simulation.

    Channel flow is a fundamental CFD validation case:
    - Fully developed flow has analytical solution (Poiseuille)
    - Tests boundary layer resolution
    - Used for turbulence model validation

    Grid requirements for channel flow:
    - Fine mesh near walls (y+ ~ 1 for wall-resolved LES)
    - Coarser mesh in core region
    - Uniform in streamwise direction (or stretched for development)

    Args:
        length: Channel length (x-direction, streamwise)
        height: Channel height (y-direction, wall-normal)
        width: Channel width (z-direction, spanwise)
        nx, ny, nz: Grid resolution
        wall_stretch: Stretching factor near walls

    Returns:
        vtkStructuredGrid: Channel flow grid with wall clustering
    """
    return create_stretched_grid(nx, ny, nz, length, height, width,
                                  stretch_y=wall_stretch, wall_clustering=True)


def add_velocity_field(grid, flow_type="channel"):
    """
    Add a velocity field to the structured grid.

    Demonstrates field data on structured grids:
    - Velocity stored at grid points (vertex-centered)
    - Vector field with 3 components (u, v, w)
    - Can represent analytical solutions or CFD results

    Args:
        grid: vtkStructuredGrid to add field to
        flow_type: Type of velocity profile ("channel", "uniform", "vortex")

    Returns:
        vtkStructuredGrid: Grid with velocity field attached
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    nx, ny, nz = dims

    velocity = vtk.vtkFloatArray()
    velocity.SetName("Velocity")
    velocity.SetNumberOfComponents(3)

    # Get grid bounds for normalization
    bounds = grid.GetBounds()
    y_min, y_max = bounds[2], bounds[3]
    height = y_max - y_min

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                point = grid.GetPoint(i + j * nx + k * nx * ny)
                x, y, z = point

                if flow_type == "channel":
                    # Parabolic velocity profile (Poiseuille flow)
                    y_norm = (y - y_min) / height
                    u = 4.0 * y_norm * (1.0 - y_norm)  # Max at center
                    v = 0.0
                    w = 0.0
                elif flow_type == "vortex":
                    # Taylor-Green vortex
                    u = math.sin(x) * math.cos(y)
                    v = -math.cos(x) * math.sin(y)
                    w = 0.0
                else:  # uniform
                    u = 1.0
                    v = 0.0
                    w = 0.0

                velocity.InsertNextTuple3(u, v, w)

    grid.GetPointData().AddArray(velocity)
    grid.GetPointData().SetActiveVectors("Velocity")

    return grid


def add_pressure_field(grid, dp_dx=-0.1):
    """
    Add a pressure field to the structured grid.

    Pressure is a fundamental CFD variable:
    - Scalar field (single value per point)
    - Drives the flow (pressure gradient)
    - Boundary condition for incompressible flow

    Args:
        grid: vtkStructuredGrid to add field to
        dp_dx: Pressure gradient in x-direction

    Returns:
        vtkStructuredGrid: Grid with pressure field attached
    """
    bounds = grid.GetBounds()
    x_min = bounds[0]

    pressure = vtk.vtkFloatArray()
    pressure.SetName("Pressure")
    pressure.SetNumberOfComponents(1)

    n_points = grid.GetNumberOfPoints()
    for i in range(n_points):
        x, y, z = grid.GetPoint(i)
        p = 1.0 + dp_dx * (x - x_min)  # Linear pressure drop
        pressure.InsertNextValue(p)

    grid.GetPointData().AddArray(pressure)
    grid.GetPointData().SetActiveScalars("Pressure")

    return grid


def compute_grid_quality(grid):
    """
    Compute and print grid quality metrics.

    Grid quality is critical for CFD accuracy:
    - Aspect ratio affects numerical diffusion
    - Orthogonality affects gradient accuracy
    - Smoothness affects solution convergence

    Args:
        grid: vtkStructuredGrid to analyze

    Returns:
        dict: Grid quality metrics
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    nx, ny, nz = dims

    print("\n" + "=" * 60)
    print("Structured Grid Quality Analysis")
    print("=" * 60)

    print(f"\nGrid dimensions: {nx} x {ny} x {nz}")
    print(f"Total points: {grid.GetNumberOfPoints()}")
    print(f"Total cells: {grid.GetNumberOfCells()}")

    bounds = grid.GetBounds()
    print(f"\nDomain bounds:")
    print(f"  X: [{bounds[0]:.4f}, {bounds[1]:.4f}]")
    print(f"  Y: [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"  Z: [{bounds[4]:.4f}, {bounds[5]:.4f}]")

    # Compute cell size statistics along j-direction (for boundary layer grids)
    if ny > 1:
        dy_values = []
        for j in range(ny - 1):
            p1 = grid.GetPoint(j * nx)
            p2 = grid.GetPoint((j + 1) * nx)
            dy = abs(p2[1] - p1[1])
            dy_values.append(dy)

        min_dy = min(dy_values)
        max_dy = max(dy_values)
        expansion_ratio = max(dy_values[i + 1] / dy_values[i]
                              for i in range(len(dy_values) - 1)) if len(dy_values) > 1 else 1.0

        print(f"\nY-direction cell sizes:")
        print(f"  First cell height: {dy_values[0]:.6f}")
        print(f"  Min cell height: {min_dy:.6f}")
        print(f"  Max cell height: {max_dy:.6f}")
        print(f"  Max expansion ratio: {expansion_ratio:.3f}")

        if expansion_ratio > MAX_RECOMMENDED_EXPANSION_RATIO:
            print(f"  WARNING: Expansion ratio > {MAX_RECOMMENDED_EXPANSION_RATIO} may affect accuracy")

    return {"dims": dims, "bounds": bounds}


def visualize_structured_grid(grid, show_edges=True, color_by_field=None):
    """
    Visualize the structured grid with CFD-style rendering.

    Demonstrates visualization techniques for structured grids:
    - Wireframe to show grid topology
    - Color mapping for scalar fields
    - Vector glyphs for velocity fields

    Args:
        grid: vtkStructuredGrid to visualize
        show_edges: Whether to show grid lines
        color_by_field: Name of scalar field for color mapping
    """
    # Create mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)

    if color_by_field and grid.GetPointData().GetArray(color_by_field):
        mapper.SetScalarModeToUsePointFieldData()
        mapper.SelectColorArray(color_by_field)
        scalar_range = grid.GetPointData().GetArray(color_by_field).GetRange()
        mapper.SetScalarRange(scalar_range)

        # Create lookup table
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)  # Blue to red
        lut.Build()
        mapper.SetLookupTable(lut)

    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    if show_edges:
        actor.GetProperty().SetEdgeVisibility(1)
        actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
        actor.GetProperty().SetLineWidth(1)

    actor.GetProperty().SetOpacity(0.8)

    # Create velocity glyphs if velocity field exists
    actors = [actor]

    if grid.GetPointData().GetArray("Velocity"):
        arrow = vtk.vtkArrowSource()

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(grid)
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetVectorModeToUseVector()
        glyph.SetInputArrayToProcess(1, 0, 0, 0, "Velocity")
        glyph.SetScaleModeToScaleByVector()
        glyph.SetScaleFactor(0.3)
        glyph.OrientOn()
        glyph.Update()

        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyph.GetOutputPort())
        glyph_mapper.ScalarVisibilityOff()

        glyph_actor = vtk.vtkActor()
        glyph_actor.SetMapper(glyph_mapper)
        glyph_actor.GetProperty().SetColor(0.2, 0.2, 0.8)

        actors.append(glyph_actor)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    for a in actors:
        renderer.AddActor(a)
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Add scalar bar if coloring by field
    if color_by_field and grid.GetPointData().GetArray(color_by_field):
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(mapper.GetLookupTable())
        scalar_bar.SetTitle(color_by_field)
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetPosition(0.85, 0.1)
        scalar_bar.SetWidth(0.1)
        scalar_bar.SetHeight(0.8)
        renderer.AddActor(scalar_bar)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Structured Grid: CFD Mesh Topology")

    # Interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Camera positioning
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Elevation(20)
    camera.Azimuth(30)

    # Add axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    interactor.Initialize()
    render_window.Render()
    interactor.Start()


def print_educational_summary():
    """
    Print educational summary about structured grids in CFD.
    """
    print("\n" + "=" * 70)
    print("VTK Structured Grid: Educational Summary for CFD")
    print("=" * 70)
    print("""
1. WHAT IS A STRUCTURED GRID?
   - Regular i-j-k topology with implicit connectivity
   - Each point identified by (i,j,k) indices
   - Grid lines are continuous across the domain
   - Points can be non-uniformly spaced

2. ADVANTAGES FOR CFD:
   - Memory efficient (no connectivity storage)
   - Fast neighbor access using indices
   - Higher order accuracy for aligned flows
   - Simple implementation of finite difference schemes

3. GRID TOPOLOGIES:
   - H-grid: Simple rectangular domain
   - O-grid: Wraps around body (airfoils, cylinders)
   - C-grid: Open at downstream (wake capture)
   - H-O-H: Combined for complex geometries

4. GRID STRETCHING:
   - Geometric: dy[i+1] = ratio * dy[i]
   - Hyperbolic tangent: clusters at both walls
   - Exponential: gradual near-wall refinement
   - Polynomial: smooth variation

5. QUALITY METRICS:
   - Aspect ratio: dx/dy (should be < 100)
   - Expansion ratio: dy[i+1]/dy[i] (should be < 1.3)
   - Orthogonality: angle between grid lines (> 30°)
   - Smoothness: gradual size changes

6. CFD APPLICATIONS:
   - Channel flow, pipe flow, boundary layers
   - External aerodynamics (with O/C grids)
   - Turbomachinery passages
   - Heat exchanger analysis
""")


def main():
    """
    Main function demonstrating VTK structured grids for CFD.

    Demonstrates:
    - Uniform grid creation
    - Boundary layer stretching
    - O-grid topology
    - Field data attachment
    - Grid quality analysis
    """
    print_educational_summary()

    # Create and analyze different grid types
    print("\n" + "-" * 70)
    print("Creating structured grid examples...")
    print("-" * 70)

    # 1. Uniform grid
    print("\n1. UNIFORM STRUCTURED GRID")
    uniform_grid = create_uniform_structured_grid(10, 8, 3, spacing=0.5)
    compute_grid_quality(uniform_grid)

    # 2. Stretched grid for channel flow
    print("\n2. CHANNEL FLOW GRID (with wall stretching)")
    channel_grid = create_channel_flow_grid(
        length=5.0, height=2.0, width=1.0,
        nx=15, ny=20, nz=5,
        wall_stretch=1.15
    )
    compute_grid_quality(channel_grid)
    add_velocity_field(channel_grid, flow_type="channel")
    add_pressure_field(channel_grid)

    # 3. O-grid around cylinder
    print("\n3. O-GRID AROUND CYLINDER")
    ogrid = create_ogrid_around_cylinder(
        radius_inner=0.5, radius_outer=3.0,
        nradial=15, ntheta=32, nz=3
    )
    compute_grid_quality(ogrid)

    # Visualize the channel flow grid with velocity field
    print("\n" + "-" * 70)
    print("Visualizing channel flow grid with velocity profile...")
    print("-" * 70)

    visualize_structured_grid(channel_grid, show_edges=True,
                               color_by_field="Pressure")


if __name__ == "__main__":
    main()
