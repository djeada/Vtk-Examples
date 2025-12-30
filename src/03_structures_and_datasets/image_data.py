"""
VTK Image Data: Uniform Grid for CFD and Image Processing

This module provides an educational demonstration of vtkImageData, which is
the simplest and most memory-efficient grid type in VTK, commonly used for
Computational Fluid Dynamics (CFD) and image-based analysis.

What is Image Data?
-------------------
vtkImageData (also known as uniform rectilinear grid) is a mesh where:
- Points are arranged in a regular 3D lattice
- Spacing is UNIFORM in each direction (dx, dy, dz are constant)
- Origin and spacing completely define point positions
- No explicit point coordinates stored (maximally memory efficient)

Image Data Properties:
----------------------
- Origin: Starting corner point (x0, y0, z0)
- Spacing: Cell size in each direction (dx, dy, dz)
- Dimensions: Number of points in each direction (nx, ny, nz)
- Extent: Index ranges [x_min, x_max, y_min, y_max, z_min, z_max]

Point position is computed: P(i,j,k) = Origin + (i*dx, j*dy, k*dz)

CFD Applications:
-----------------
Image data is used extensively in CFD for:

1. CARTESIAN GRID METHODS:
   - Immersed boundary method (IBM)
   - Cut-cell methods for complex geometries
   - Volume of Fluid (VOF) for multiphase

2. DNS/LES SIMULATIONS:
   - Homogeneous isotropic turbulence
   - Channel flow with periodic BCs
   - Spectral methods (FFT-based solvers)

3. POST-PROCESSING:
   - Sampling CFD results on regular grids
   - Interpolation between different meshes
   - Data export for visualization tools

4. IMAGE-BASED CFD:
   - Medical imaging (CT/MRI) to flow simulation
   - Porous media flow from microCT
   - Cardiovascular flow analysis

Advantages:
-----------
- Minimal memory footprint (only store origin, spacing, dimensions)
- Fast point/cell location (O(1) using indices)
- Simple implementation of stencil operations
- Natural for FFT-based methods

This example demonstrates creating image data for CFD applications,
including volumetric data and field visualization.
"""

import math

import vtk


def create_uniform_image_data(nx, ny, nz, spacing=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0)):
    """
    Create a uniform vtkImageData grid.

    This is the most memory-efficient VTK grid type:
    - Only stores dimensions, spacing, and origin
    - Point positions computed on-the-fly
    - Ideal for large uniform grids

    Args:
        nx, ny, nz: Number of points in each direction
        spacing: (dx, dy, dz) tuple of cell sizes
        origin: (x0, y0, z0) tuple for grid origin

    Returns:
        vtkImageData: Uniform grid
    """
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(nx, ny, nz)
    image_data.SetSpacing(*spacing)
    image_data.SetOrigin(*origin)

    return image_data


def create_channel_flow_domain(length, height, width, resolution):
    """
    Create image data representing a channel flow domain.

    Channel flow is a fundamental CFD validation case:
    - Rectangular domain with walls at y=0 and y=height
    - Periodic in x (streamwise) and z (spanwise)
    - Uniform grid suitable for spectral methods

    Args:
        length: Domain length (x-direction)
        height: Domain height (y-direction)
        width: Domain width (z-direction)
        resolution: Grid points per unit length

    Returns:
        vtkImageData: Channel flow domain
    """
    nx = int(length * resolution) + 1
    ny = int(height * resolution) + 1
    nz = int(width * resolution) + 1

    dx = length / max(nx - 1, 1)
    dy = height / max(ny - 1, 1)
    dz = width / max(nz - 1, 1)

    return create_uniform_image_data(nx, ny, nz, spacing=(dx, dy, dz))


def create_box_turbulence_domain(size, n_points):
    """
    Create image data for homogeneous isotropic turbulence simulation.

    Box turbulence (HIT) simulations:
    - Cubic domain with periodic boundaries
    - Uniform grid for spectral accuracy
    - Common DNS configuration

    Args:
        size: Cubic domain side length
        n_points: Number of points in each direction

    Returns:
        vtkImageData: Cubic domain for turbulence simulation
    """
    spacing = size / n_points

    return create_uniform_image_data(
        n_points, n_points, n_points,
        spacing=(spacing, spacing, spacing)
    )


def add_velocity_field(image_data, flow_type="channel"):
    """
    Add a velocity field to the image data.

    Demonstrates vector field storage on uniform grids:
    - 3-component vector at each point
    - Common for CFD velocity data
    - Can represent analytical or computed solutions

    Args:
        image_data: vtkImageData to modify
        flow_type: Type of flow ("channel", "vortex", "uniform")

    Returns:
        vtkImageData: Grid with velocity field
    """
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()
    origin = image_data.GetOrigin()

    height = (dims[1] - 1) * spacing[1]

    velocity = vtk.vtkFloatArray()
    velocity.SetName("Velocity")
    velocity.SetNumberOfComponents(3)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                x = origin[0] + i * spacing[0]
                y = origin[1] + j * spacing[1]
                z = origin[2] + k * spacing[2]

                if flow_type == "channel":
                    # Poiseuille flow profile
                    y_norm = y / height
                    u = 4.0 * y_norm * (1.0 - y_norm)
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

    image_data.GetPointData().AddArray(velocity)
    image_data.GetPointData().SetActiveVectors("Velocity")

    return image_data


def add_scalar_field(image_data, field_name="Temperature", mode="gradient"):
    """
    Add a scalar field to the image data.

    Scalar fields are common in CFD:
    - Temperature for heat transfer
    - Pressure (though often at cell centers)
    - Concentration for species transport

    Args:
        image_data: vtkImageData to modify
        field_name: Name of the scalar field
        mode: Field distribution type

    Returns:
        vtkImageData: Grid with scalar field
    """
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()
    origin = image_data.GetOrigin()

    # Calculate domain size
    length = (dims[0] - 1) * spacing[0]
    height = (dims[1] - 1) * spacing[1]

    scalars = vtk.vtkFloatArray()
    scalars.SetName(field_name)
    scalars.SetNumberOfComponents(1)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                x = origin[0] + i * spacing[0]
                y = origin[1] + j * spacing[1]
                z = origin[2] + k * spacing[2]

                x_norm = x / max(length, 1e-10)
                y_norm = y / max(height, 1e-10)

                if mode == "gradient":
                    # Linear gradient
                    value = 300 + 100 * x_norm
                elif mode == "hot_spot":
                    # Gaussian hot spot
                    r2 = (x_norm - 0.5)**2 + (y_norm - 0.5)**2
                    value = 300 + 200 * math.exp(-20 * r2)
                elif mode == "layers":
                    # Horizontal temperature layers
                    value = 300 + 100 * y_norm
                else:
                    value = 350.0

                scalars.InsertNextValue(value)

    image_data.GetPointData().AddArray(scalars)
    image_data.GetPointData().SetActiveScalars(field_name)

    return image_data


def add_vorticity_magnitude(image_data):
    """
    Compute and add vorticity magnitude from velocity field.

    Vorticity is a key quantity in CFD:
    - Measures local fluid rotation
    - ω = ∇ × u (curl of velocity)
    - Used to identify vortex structures

    This computes a simplified 2D vorticity for demonstration.

    Args:
        image_data: vtkImageData with velocity field

    Returns:
        vtkImageData: Grid with vorticity magnitude
    """
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()

    velocity = image_data.GetPointData().GetArray("Velocity")
    if not velocity:
        return image_data

    vorticity = vtk.vtkFloatArray()
    vorticity.SetName("VorticityMagnitude")
    vorticity.SetNumberOfComponents(1)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                idx = i + j * dims[0] + k * dims[0] * dims[1]

                # Compute gradients using central differences
                if 0 < i < dims[0] - 1 and 0 < j < dims[1] - 1:
                    idx_ip = (i + 1) + j * dims[0] + k * dims[0] * dims[1]
                    idx_im = (i - 1) + j * dims[0] + k * dims[0] * dims[1]
                    idx_jp = i + (j + 1) * dims[0] + k * dims[0] * dims[1]
                    idx_jm = i + (j - 1) * dims[0] + k * dims[0] * dims[1]

                    u_ip = velocity.GetTuple3(idx_ip)
                    u_im = velocity.GetTuple3(idx_im)
                    u_jp = velocity.GetTuple3(idx_jp)
                    u_jm = velocity.GetTuple3(idx_jm)

                    # ω_z = dv/dx - du/dy
                    dvdx = (u_ip[1] - u_im[1]) / (2 * spacing[0])
                    dudy = (u_jp[0] - u_jm[0]) / (2 * spacing[1])

                    omega_z = dvdx - dudy
                else:
                    omega_z = 0.0

                vorticity.InsertNextValue(abs(omega_z))

    image_data.GetPointData().AddArray(vorticity)

    return image_data


def compute_grid_info(image_data):
    """
    Compute and print image data grid information.

    Args:
        image_data: vtkImageData to analyze
    """
    dims = image_data.GetDimensions()
    spacing = image_data.GetSpacing()
    origin = image_data.GetOrigin()
    bounds = image_data.GetBounds()

    print("\n" + "=" * 60)
    print("Image Data Grid Information")
    print("=" * 60)

    print(f"\nDimensions: {dims[0]} x {dims[1]} x {dims[2]}")
    print(f"Total points: {image_data.GetNumberOfPoints()}")
    print(f"Total cells: {image_data.GetNumberOfCells()}")

    print(f"\nOrigin: ({origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f})")
    print(f"Spacing: ({spacing[0]:.4f}, {spacing[1]:.4f}, {spacing[2]:.4f})")

    print(f"\nBounds:")
    print(f"  X: [{bounds[0]:.4f}, {bounds[1]:.4f}]")
    print(f"  Y: [{bounds[2]:.4f}, {bounds[3]:.4f}]")
    print(f"  Z: [{bounds[4]:.4f}, {bounds[5]:.4f}]")

    # Memory estimation
    bytes_per_point = 0
    pd = image_data.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(i)
        bytes_per_point += arr.GetNumberOfComponents() * 4  # Assuming float32

    total_memory_mb = (image_data.GetNumberOfPoints() * bytes_per_point) / (1024 * 1024)
    print(f"\nEstimated field data memory: {total_memory_mb:.2f} MB")

    print("\nPoint data arrays:")
    if pd.GetNumberOfArrays() == 0:
        print("  (none)")
    for i in range(pd.GetNumberOfArrays()):
        arr = pd.GetArray(i)
        print(f"  {arr.GetName()}: {arr.GetNumberOfTuples()} values, "
              f"{arr.GetNumberOfComponents()} components")


def visualize_image_data(image_data, color_by_field=None, show_outline=True):
    """
    Visualize image data with CFD-style rendering.

    Args:
        image_data: vtkImageData to visualize
        color_by_field: Scalar field for color mapping
        show_outline: Whether to show domain outline
    """
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Show domain outline
    if show_outline:
        outline = vtk.vtkOutlineFilter()
        outline.SetInputData(image_data)

        outline_mapper = vtk.vtkPolyDataMapper()
        outline_mapper.SetInputConnection(outline.GetOutputPort())

        outline_actor = vtk.vtkActor()
        outline_actor.SetMapper(outline_mapper)
        outline_actor.GetProperty().SetColor(1, 1, 1)
        outline_actor.GetProperty().SetLineWidth(2)

        renderer.AddActor(outline_actor)

    # Create slice through the middle
    dims = image_data.GetDimensions()

    # XY slice
    plane = vtk.vtkPlane()
    center = image_data.GetCenter()
    plane.SetOrigin(center)
    plane.SetNormal(0, 0, 1)

    cutter = vtk.vtkCutter()
    cutter.SetInputData(image_data)
    cutter.SetCutFunction(plane)
    cutter.Update()

    slice_mapper = vtk.vtkPolyDataMapper()
    slice_mapper.SetInputConnection(cutter.GetOutputPort())

    if color_by_field and image_data.GetPointData().GetArray(color_by_field):
        slice_mapper.SetScalarModeToUsePointFieldData()
        slice_mapper.SelectColorArray(color_by_field)
        scalar_range = image_data.GetPointData().GetArray(color_by_field).GetRange()
        slice_mapper.SetScalarRange(scalar_range)

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)
        lut.Build()
        slice_mapper.SetLookupTable(lut)

    slice_actor = vtk.vtkActor()
    slice_actor.SetMapper(slice_mapper)
    renderer.AddActor(slice_actor)

    # Add scalar bar
    if color_by_field and image_data.GetPointData().GetArray(color_by_field):
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(slice_mapper.GetLookupTable())
        scalar_bar.SetTitle(color_by_field)
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetPosition(0.85, 0.1)
        scalar_bar.SetWidth(0.1)
        scalar_bar.SetHeight(0.8)
        renderer.AddActor(scalar_bar)

    # Add velocity glyphs if available
    if image_data.GetPointData().GetArray("Velocity"):
        # Subsample for glyphs
        arrow = vtk.vtkArrowSource()

        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(image_data)
        glyph.SetSourceConnection(arrow.GetOutputPort())
        glyph.SetVectorModeToUseVector()
        glyph.SetInputArrayToProcess(1, 0, 0, 0, "Velocity")
        glyph.SetScaleModeToScaleByVector()
        glyph.SetScaleFactor(0.2)
        glyph.OrientOn()
        glyph.Update()

        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyph.GetOutputPort())
        glyph_mapper.ScalarVisibilityOff()

        glyph_actor = vtk.vtkActor()
        glyph_actor.SetMapper(glyph_mapper)
        glyph_actor.GetProperty().SetColor(0.2, 0.8, 0.2)

        renderer.AddActor(glyph_actor)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Image Data: Uniform CFD Grid")

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
    Print educational summary about image data in CFD.
    """
    print("\n" + "=" * 70)
    print("VTK Image Data: Educational Summary for CFD")
    print("=" * 70)
    print("""
1. WHAT IS IMAGE DATA?
   - Uniform rectilinear grid (constant spacing)
   - Defined by origin, spacing, and dimensions
   - No explicit point storage (computed on-the-fly)
   - Most memory-efficient VTK grid type

2. KEY PROPERTIES:
   - Origin: (x0, y0, z0) starting corner
   - Spacing: (dx, dy, dz) cell sizes
   - Dimensions: (nx, ny, nz) point counts
   - Point(i,j,k) = Origin + (i*dx, j*dy, k*dz)

3. CFD APPLICATIONS:
   - Cartesian grid methods (IBM, cut-cell)
   - DNS/LES with spectral methods
   - Turbulence box simulations
   - Post-processing and data exchange

4. ADVANTAGES:
   - Minimal memory for grid structure
   - O(1) point/cell location
   - Simple stencil operations
   - Natural for FFT-based solvers

5. LIMITATIONS:
   - Uniform spacing only (use RectilinearGrid for stretching)
   - Axis-aligned cells only
   - Not suitable for complex geometries directly

6. COMMON OPERATIONS:
   - Slicing and probing
   - Contour extraction
   - Volume rendering
   - Resampling between grids
""")


def main():
    """
    Main function demonstrating VTK image data for CFD.

    Demonstrates:
    - Creating uniform grids for CFD
    - Adding velocity and scalar fields
    - Computing derived quantities (vorticity)
    - Visualization with slices and glyphs
    """
    print_educational_summary()

    # Create channel flow domain
    print("\n" + "-" * 70)
    print("Creating channel flow domain...")
    print("-" * 70)

    channel = create_channel_flow_domain(
        length=4.0, height=2.0, width=2.0,
        resolution=10
    )

    # Add fields
    add_velocity_field(channel, flow_type="channel")
    add_scalar_field(channel, "Temperature", mode="gradient")
    add_vorticity_magnitude(channel)

    compute_grid_info(channel)

    # Create box turbulence domain
    print("\n" + "-" * 70)
    print("Creating box turbulence domain...")
    print("-" * 70)

    box = create_box_turbulence_domain(size=2.0 * math.pi, n_points=32)
    add_velocity_field(box, flow_type="vortex")

    compute_grid_info(box)

    # Visualize channel flow
    print("\n" + "-" * 70)
    print("Visualizing channel flow domain...")
    print("-" * 70)

    visualize_image_data(channel, color_by_field="Temperature")


if __name__ == "__main__":
    main()
