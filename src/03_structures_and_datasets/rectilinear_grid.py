"""
VTK Rectilinear Grid: Axis-Aligned Non-Uniform Mesh for CFD

This module provides an educational demonstration of vtkRectilinearGrid, which
is an efficient data structure for axis-aligned meshes with non-uniform spacing
in Computational Fluid Dynamics (CFD).

What is a Rectilinear Grid?
---------------------------
A rectilinear grid is a mesh where:
- Grid lines are parallel to coordinate axes (X, Y, Z)
- Spacing can be non-uniform in each direction
- Points are defined by three 1D coordinate arrays
- Topology is implicit (like structured grid)
- Cells are always axis-aligned rectangular hexahedra

Rectilinear vs Other Grid Types:
--------------------------------
RECTILINEAR GRID (vtkRectilinearGrid):
- Axis-aligned cells only
- Non-uniform spacing allowed
- Defined by 1D coordinate arrays (memory efficient)
- Fast interpolation and lookup

IMAGE DATA (vtkImageData):
- Uniform spacing required
- Most memory efficient
- Simplest grid type

STRUCTURED GRID (vtkStructuredGrid):
- Can have curved grid lines
- Most flexible structured type
- Higher memory (stores all 3D points)

CFD Applications:
-----------------
Rectilinear grids are commonly used for:

1. GRID STRETCHING:
   - Clustering near walls with simple specification
   - Each direction stretched independently
   - No grid skewness (always orthogonal)

2. CARTESIAN METHODS:
   - Immersed boundary methods
   - Cut-cell methods
   - Adaptive mesh refinement (AMR)

3. DNS/LES SIMULATIONS:
   - Spectral methods require structured grids
   - Box turbulence simulations
   - Channel flow (periodic in streamwise/spanwise)

4. HEAT TRANSFER:
   - Rectangular enclosures
   - Heat exchangers (regular geometries)
   - Building energy simulation

Advantages:
-----------
- Memory efficient (only store 1D arrays)
- Simple indexing and fast interpolation
- No grid skewness (always orthogonal)
- Easy to implement finite difference methods

This example demonstrates creating rectilinear grids with various stretching
functions commonly used in CFD simulations.
"""

import math

import vtk


# Default stretching parameters
DEFAULT_STRETCH_RATIO = 1.15
MIN_CELL_SIZE = 0.001


def create_uniform_rectilinear_grid(nx, ny, nz, spacing=1.0):
    """
    Create a uniform rectilinear grid (equivalent to ImageData but more flexible).

    Uniform grids are the baseline for many CFD methods:
    - Equal spacing in all directions
    - Simple finite difference implementation
    - Good for DNS with spectral accuracy

    Args:
        nx, ny, nz: Number of points in each direction
        spacing: Uniform spacing between points

    Returns:
        vtkRectilinearGrid: Uniform rectilinear grid
    """
    # Create 1D coordinate arrays
    x_coords = vtk.vtkFloatArray()
    for i in range(nx):
        x_coords.InsertNextValue(i * spacing)

    y_coords = vtk.vtkFloatArray()
    for j in range(ny):
        y_coords.InsertNextValue(j * spacing)

    z_coords = vtk.vtkFloatArray()
    for k in range(nz):
        z_coords.InsertNextValue(k * spacing)

    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(nx, ny, nz)
    rgrid.SetXCoordinates(x_coords)
    rgrid.SetYCoordinates(y_coords)
    rgrid.SetZCoordinates(z_coords)

    return rgrid


def create_geometric_stretch_coords(n_points, length, stretch_ratio, reverse=False):
    """
    Create 1D coordinate array with geometric stretching.

    Geometric stretching is the most common CFD mesh stretching:
    - Each cell is stretch_ratio times larger than previous
    - Good for boundary layer resolution
    - Simple to specify and understand

    The cell sizes follow: dx[i] = dx[0] * ratio^i

    Args:
        n_points: Number of grid points
        length: Total length to cover
        stretch_ratio: Ratio between successive cell sizes
        reverse: If True, cluster at end instead of beginning

    Returns:
        vtkFloatArray: 1D coordinate array
    """
    coords = vtk.vtkFloatArray()

    if n_points <= 1:
        coords.InsertNextValue(0.0)
        return coords

    # Calculate first cell size for geometric series
    # Sum of geometric series: S = dx[0] * (ratio^n - 1) / (ratio - 1)
    if abs(stretch_ratio - 1.0) < 1e-10:
        # Uniform spacing
        dx0 = length / (n_points - 1)
    else:
        n_cells = n_points - 1
        dx0 = length * (stretch_ratio - 1) / (stretch_ratio**n_cells - 1)

    # Generate coordinates
    positions = [0.0]
    x = 0.0
    for i in range(n_points - 1):
        dx = dx0 * (stretch_ratio**i)
        x += dx
        positions.append(x)

    # Reverse if clustering at end
    if reverse:
        positions = [length - p for p in reversed(positions)]

    for p in positions:
        coords.InsertNextValue(p)

    return coords


def create_tanh_stretch_coords(n_points, length, beta=1.5):
    """
    Create 1D coordinate array with hyperbolic tangent stretching.

    Tanh stretching clusters points at BOTH ends:
    - Natural for channel flow (walls at both ends)
    - Symmetric distribution
    - Smooth variation of cell sizes

    The stretching function: x = L/2 * (1 + tanh(beta*(eta-0.5))/tanh(beta/2))
    where eta goes from 0 to 1

    Args:
        n_points: Number of grid points
        length: Total length
        beta: Stretching parameter (higher = more clustering)

    Returns:
        vtkFloatArray: 1D coordinate array with two-sided clustering
    """
    coords = vtk.vtkFloatArray()

    for i in range(n_points):
        eta = i / max(n_points - 1, 1)
        x = length / 2 * (1 + math.tanh(beta * (eta - 0.5)) / math.tanh(beta / 2))
        coords.InsertNextValue(x)

    return coords


def create_boundary_layer_grid(
    length_x, length_y, length_z, nx, ny, nz, y_stretch=1.2, cluster_both_walls=True
):
    """
    Create a rectilinear grid optimized for boundary layer resolution.

    Boundary layer grids require:
    - Fine mesh near walls to resolve velocity gradients
    - First cell height controls y+ (wall units)
    - Growth rate affects boundary layer resolution

    For wall-resolved LES: y+ ~ 1 (first cell in viscous sublayer)
    For wall-modeled LES: y+ ~ 30-50 (first cell in log layer)
    For RANS: depends on wall treatment

    Args:
        length_x, length_y, length_z: Domain dimensions
        nx, ny, nz: Number of points in each direction
        y_stretch: Stretching ratio in y-direction
        cluster_both_walls: If True, cluster at y=0 and y=length_y

    Returns:
        vtkRectilinearGrid: Grid with boundary layer clustering
    """
    # X-direction: uniform (streamwise, often periodic)
    x_coords = vtk.vtkFloatArray()
    for i in range(nx):
        x_coords.InsertNextValue(i * length_x / max(nx - 1, 1))

    # Y-direction: stretched for boundary layers
    if cluster_both_walls:
        y_coords = create_tanh_stretch_coords(ny, length_y, beta=1.5)
    else:
        y_coords = create_geometric_stretch_coords(ny, length_y, y_stretch)

    # Z-direction: uniform (spanwise, often periodic)
    z_coords = vtk.vtkFloatArray()
    for k in range(nz):
        z_coords.InsertNextValue(k * length_z / max(nz - 1, 1))

    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(nx, ny, nz)
    rgrid.SetXCoordinates(x_coords)
    rgrid.SetYCoordinates(y_coords)
    rgrid.SetZCoordinates(z_coords)

    return rgrid


def create_refinement_region_grid(length, n_base, n_refined, refine_start, refine_end):
    """
    Create coordinate array with local refinement region.

    Local refinement is useful for:
    - Wake regions behind bodies
    - Shear layers and mixing zones
    - Regions with expected high gradients

    Args:
        length: Total length
        n_base: Base number of cells
        n_refined: Number of cells in refined region
        refine_start: Start of refined region (0 to 1)
        refine_end: End of refined region (0 to 1)

    Returns:
        vtkFloatArray: Coordinate array with refinement
    """
    coords = vtk.vtkFloatArray()

    # Calculate regions
    start_pos = refine_start * length
    end_pos = refine_end * length

    # Cells before refinement
    n_before = int(n_base * refine_start)
    dx_before = start_pos / max(n_before, 1)

    # Cells in refined region
    dx_refined = (end_pos - start_pos) / n_refined

    # Cells after refinement
    n_after = n_base - n_before
    dx_after = (length - end_pos) / max(n_after, 1)

    # Generate coordinates
    x = 0.0
    coords.InsertNextValue(x)

    # Before region
    for _ in range(n_before):
        x += dx_before
        coords.InsertNextValue(x)

    # Refined region
    x = start_pos
    for _ in range(n_refined):
        x += dx_refined
        coords.InsertNextValue(x)

    # After region
    x = end_pos
    for _ in range(n_after):
        x += dx_after
        coords.InsertNextValue(min(x, length))

    return coords


def add_scalar_field(rgrid, field_name="Temperature", mode="gradient"):
    """
    Add a scalar field to the rectilinear grid.

    Demonstrates field data on rectilinear grids:
    - Scalar stored at each grid point
    - Easy indexing using i,j,k structure
    - Interpolation uses grid structure

    Args:
        rgrid: vtkRectilinearGrid to modify
        field_name: Name of the scalar field
        mode: Type of field ("gradient", "sinusoidal", "gaussian")

    Returns:
        vtkRectilinearGrid: Grid with scalar field
    """
    dims = rgrid.GetDimensions()
    bounds = rgrid.GetBounds()

    scalars = vtk.vtkFloatArray()
    scalars.SetName(field_name)
    scalars.SetNumberOfComponents(1)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                idx = i + j * dims[0] + k * dims[0] * dims[1]
                x, y, z = rgrid.GetPoint(idx)

                # Normalize coordinates
                x_norm = (x - bounds[0]) / max(bounds[1] - bounds[0], 1e-10)
                y_norm = (y - bounds[2]) / max(bounds[3] - bounds[2], 1e-10)

                if mode == "gradient":
                    # Linear temperature gradient
                    value = 100 + 200 * x_norm
                elif mode == "sinusoidal":
                    # Sinusoidal variation
                    value = math.sin(2 * math.pi * x_norm) * math.cos(math.pi * y_norm)
                elif mode == "gaussian":
                    # Gaussian hot spot
                    r2 = (x_norm - 0.5) ** 2 + (y_norm - 0.5) ** 2
                    value = math.exp(-10 * r2)
                else:
                    value = x_norm + y_norm

                scalars.InsertNextValue(value)

    rgrid.GetPointData().AddArray(scalars)
    rgrid.GetPointData().SetActiveScalars(field_name)

    return rgrid


def add_velocity_field(rgrid):
    """
    Add a velocity field representing channel flow.

    Channel flow velocity profile (Poiseuille flow):
    - Parabolic profile in wall-normal direction
    - Zero at walls, maximum at center
    - Exact solution for laminar fully-developed flow

    Args:
        rgrid: vtkRectilinearGrid to modify

    Returns:
        vtkRectilinearGrid: Grid with velocity field
    """
    dims = rgrid.GetDimensions()
    bounds = rgrid.GetBounds()
    y_min, y_max = bounds[2], bounds[3]
    height = y_max - y_min

    velocity = vtk.vtkFloatArray()
    velocity.SetName("Velocity")
    velocity.SetNumberOfComponents(3)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                idx = i + j * dims[0] + k * dims[0] * dims[1]
                x, y, z = rgrid.GetPoint(idx)

                # Parabolic profile: u = u_max * (1 - (2y/H - 1)^2)
                y_norm = (y - y_min) / height
                u = 4.0 * y_norm * (1.0 - y_norm)  # Max = 1 at center
                v = 0.0
                w = 0.0

                velocity.InsertNextTuple3(u, v, w)

    rgrid.GetPointData().AddArray(velocity)
    rgrid.GetPointData().SetActiveVectors("Velocity")

    return rgrid


def compute_grid_statistics(rgrid):
    """
    Compute and print grid statistics.

    Args:
        rgrid: vtkRectilinearGrid to analyze
    """
    dims = rgrid.GetDimensions()

    print("\n" + "=" * 60)
    print("Rectilinear Grid Statistics")
    print("=" * 60)

    print(f"\nDimensions: {dims[0]} x {dims[1]} x {dims[2]}")
    print(f"Total points: {rgrid.GetNumberOfPoints()}")
    print(f"Total cells: {rgrid.GetNumberOfCells()}")

    bounds = rgrid.GetBounds()
    print(f"\nDomain bounds:")
    print(f"  X: [{bounds[0]:.4f}, {bounds[1]:.4f}]")
    print(f"  Y: [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"  Z: [{bounds[4]:.4f}, {bounds[5]:.4f}]")

    # Analyze cell sizes in each direction
    for direction, name in [(0, "X"), (1, "Y"), (2, "Z")]:
        if direction == 0:
            coords = rgrid.GetXCoordinates()
        elif direction == 1:
            coords = rgrid.GetYCoordinates()
        else:
            coords = rgrid.GetZCoordinates()

        n = coords.GetNumberOfTuples()
        if n > 1:
            sizes = []
            for i in range(n - 1):
                d = coords.GetValue(i + 1) - coords.GetValue(i)
                sizes.append(d)

            min_size = min(sizes)
            max_size = max(sizes)
            ratios = (
                [sizes[i + 1] / sizes[i] for i in range(len(sizes) - 1)]
                if len(sizes) > 1
                else [1.0]
            )
            max_ratio = max(ratios)

            print(f"\n{name}-direction cell sizes:")
            print(f"  Min: {min_size:.6f}")
            print(f"  Max: {max_size:.6f}")
            print(f"  Max expansion ratio: {max_ratio:.4f}")


def visualize_rectilinear_grid(rgrid, show_edges=True, color_by_field=None):
    """
    Visualize the rectilinear grid with CFD-style rendering.

    Args:
        rgrid: vtkRectilinearGrid to visualize
        show_edges: Whether to show grid lines
        color_by_field: Name of field for color mapping
    """
    # Convert to polydata for rendering
    geometry = vtk.vtkRectilinearGridGeometryFilter()
    geometry.SetInputData(rgrid)
    geometry.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometry.GetOutputPort())

    if color_by_field and rgrid.GetPointData().GetArray(color_by_field):
        mapper.SetScalarModeToUsePointFieldData()
        mapper.SelectColorArray(color_by_field)
        scalar_range = rgrid.GetPointData().GetArray(color_by_field).GetRange()
        mapper.SetScalarRange(scalar_range)

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)
        lut.Build()
        mapper.SetLookupTable(lut)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    if show_edges:
        actor.GetProperty().SetEdgeVisibility(1)
        actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
        actor.GetProperty().SetLineWidth(1)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.15, 0.15, 0.2)

    if color_by_field and rgrid.GetPointData().GetArray(color_by_field):
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(mapper.GetLookupTable())
        scalar_bar.SetTitle(color_by_field)
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetPosition(0.85, 0.1)
        scalar_bar.SetWidth(0.1)
        scalar_bar.SetHeight(0.8)
        renderer.AddActor(scalar_bar)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Rectilinear Grid: Non-Uniform CFD Mesh")

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Elevation(30)
    camera.Azimuth(30)

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
    Print educational summary about rectilinear grids in CFD.
    """
    print("\n" + "=" * 70)
    print("VTK Rectilinear Grid: Educational Summary for CFD")
    print("=" * 70)
    print(
        """
1. WHAT IS A RECTILINEAR GRID?
   - Axis-aligned grid with non-uniform spacing
   - Defined by three 1D coordinate arrays (x, y, z)
   - Cells are always rectangular hexahedra
   - Implicit connectivity (like structured grid)

2. ADVANTAGES:
   - Memory efficient (only 1D arrays stored)
   - Always orthogonal (no grid skewness)
   - Simple finite difference implementation
   - Fast interpolation and lookup

3. STRETCHING FUNCTIONS:
   - Geometric: dx[i+1] = ratio * dx[i]
   - Hyperbolic tangent: clusters at both ends
   - Exponential: smooth near-wall clustering
   - Piecewise: different regions with different spacing

4. CFD APPLICATIONS:
   - Channel and duct flows
   - Heat transfer in enclosures
   - DNS/LES with spectral methods
   - Cartesian immersed boundary methods

5. GRID QUALITY:
   - Expansion ratio should be < 1.2-1.3
   - Aspect ratio affects numerical diffusion
   - Smooth transitions between regions
   - First cell size controls y+

6. COMPARISON WITH OTHER GRIDS:
   - vs ImageData: allows non-uniform spacing
   - vs StructuredGrid: restricted to axis-aligned
   - vs UnstructuredGrid: more memory efficient
"""
    )


def main():
    """
    Main function demonstrating VTK rectilinear grids for CFD.

    Demonstrates:
    - Uniform grid creation
    - Geometric stretching for boundary layers
    - Hyperbolic tangent stretching for channel flow
    - Field data attachment
    - Grid statistics and visualization
    """
    print_educational_summary()

    # Create different grid types
    print("\n" + "-" * 70)
    print("Creating rectilinear grid examples...")
    print("-" * 70)

    # 1. Uniform grid
    print("\n1. UNIFORM RECTILINEAR GRID")
    uniform = create_uniform_rectilinear_grid(10, 8, 4, spacing=0.5)
    compute_grid_statistics(uniform)

    # 2. Boundary layer grid
    print("\n2. BOUNDARY LAYER GRID (channel flow)")
    bl_grid = create_boundary_layer_grid(
        length_x=10.0,
        length_y=2.0,
        length_z=3.0,
        nx=20,
        ny=30,
        nz=8,
        y_stretch=1.15,
        cluster_both_walls=True,
    )
    compute_grid_statistics(bl_grid)

    # Add fields to boundary layer grid
    add_scalar_field(bl_grid, "Temperature", "gradient")
    add_velocity_field(bl_grid)

    # 3. Visualize
    print("\n" + "-" * 70)
    print("Visualizing boundary layer grid...")
    print("-" * 70)

    visualize_rectilinear_grid(bl_grid, show_edges=True, color_by_field="Temperature")


if __name__ == "__main__":
    main()
