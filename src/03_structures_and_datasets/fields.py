"""
VTK Fields: Data Attributes in Computational Meshes

This module provides an educational demonstration of VTK fields, which represent
physical quantities in Computational Fluid Dynamics (CFD) and Finite Element
Analysis (FEA) simulations.

What are Fields?
----------------
In VTK (and numerical simulations), fields are data arrays that store physical
quantities associated with the geometric elements of a mesh. Fields can be:
- Scalar fields: Single values (temperature, pressure, density)
- Vector fields: 3-component values (velocity, force, displacement)
- Tensor fields: 9-component values (stress, strain tensors)

Field Location (Point Data vs Cell Data):
-----------------------------------------
Fields can be defined at different locations in the mesh:

1. POINT DATA (Vertex-centered):
   - Values stored at mesh vertices/nodes
   - Interpolated within cells for visualization
   - Common in FEA where quantities are solved at nodes
   - Examples: nodal displacements, temperature at nodes

2. CELL DATA (Cell-centered):
   - Values stored at cell centers
   - Constant within each cell (piecewise constant)
   - Common in Finite Volume CFD methods
   - Examples: cell pressure, cell-averaged velocity

CFD/FEA Context:
----------------
Understanding fields is crucial for CFD/FEA because:

1. Conservation equations are discretized as:
   - FVM (Finite Volume Method): Cell-centered values, fluxes at faces
   - FEM (Finite Element Method): Node-centered values, shape functions

2. Post-processing requires mapping between:
   - Simulation results (solver-native location)
   - Visualization requirements (may need interpolation)

3. Physical interpretation:
   - Scalar fields: Contour plots, isosurfaces
   - Vector fields: Streamlines, glyphs (arrows)
   - Tensor fields: Principal stress directions

This example demonstrates creating a structured grid with multiple physical
fields commonly found in CFD simulations (velocity, pressure, temperature)
and visualizing them using appropriate techniques.
"""

import math

import numpy as np
import vtk


def create_structured_grid(nx, ny, nz, spacing=1.0):
    """
    Create a structured grid representing a computational domain.

    Structured grids are commonly used in CFD for:
    - Channel flows and pipe flows (simple geometries)
    - Boundary layer meshes (aligned with flow)
    - DNS/LES simulations (regular spacing for spectral accuracy)

    Args:
        nx, ny, nz: Number of points in each direction
        spacing: Distance between grid points

    Returns:
        vtkStructuredGrid: The computational mesh
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


def create_velocity_field(grid, u_bulk=1.0):
    """
    Create a velocity field representing Poiseuille flow (channel flow).

    Poiseuille flow is the exact solution for fully-developed laminar flow
    between parallel plates or in a pipe. The velocity profile is parabolic:
    - Maximum velocity at the center
    - Zero velocity at the walls (no-slip boundary condition)

    This is a fundamental CFD validation case because:
    - Analytical solution exists for comparison
    - Tests boundary condition implementation
    - Demonstrates viscous effects

    The velocity profile is: u(y) = u_max * (1 - (2y/H)^2)
    where H is the channel height and u_max = 1.5 * u_bulk

    Args:
        grid: vtkStructuredGrid to attach the field to
        u_bulk: Bulk (average) velocity

    Returns:
        vtkFloatArray: Velocity vectors at each point
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    ny = dims[1]
    h = (ny - 1)  # Channel height in grid units
    u_max = 1.5 * u_bulk  # Maximum velocity for parabolic profile

    velocity = vtk.vtkFloatArray()
    velocity.SetName("Velocity")
    velocity.SetNumberOfComponents(3)

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                # Parabolic profile: u = u_max * (1 - (2*y/H - 1)^2)
                y_normalized = 2.0 * j / h - 1.0  # Range [-1, 1]
                u_x = u_max * (1.0 - y_normalized * y_normalized)
                u_y = 0.0
                u_z = 0.0
                velocity.InsertNextTuple3(u_x, u_y, u_z)

    return velocity


def create_pressure_field(grid, dp_dx=-0.1):
    """
    Create a pressure field with linear gradient (driving force for flow).

    In Poiseuille flow, the pressure gradient is the driving force that
    balances viscous forces. The pressure decreases linearly in the flow
    direction:
    - Constant pressure gradient dp/dx
    - Related to viscosity and velocity by: dp/dx = -12 * mu * u_bulk / H^2

    Understanding pressure fields is essential for CFD because:
    - Pressure is the primary variable in incompressible flow
    - Pressure gradients drive the flow
    - Pressure boundary conditions are critical for convergence

    Args:
        grid: vtkStructuredGrid to attach the field to
        dp_dx: Pressure gradient in x-direction (negative = flow in +x)

    Returns:
        vtkFloatArray: Pressure scalar at each point
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    pressure = vtk.vtkFloatArray()
    pressure.SetName("Pressure")
    pressure.SetNumberOfComponents(1)

    # Reference pressure at inlet
    p_inlet = 1.0

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                # Linear pressure drop: p = p_inlet + dp/dx * x
                x = i
                p = p_inlet + dp_dx * x
                pressure.InsertNextValue(p)

    return pressure


def create_temperature_field(grid, t_wall=100.0, t_inlet=20.0):
    """
    Create a temperature field for heated channel flow.

    This represents thermal development in a channel with heated walls.
    The temperature profile develops from uniform inlet temperature to
    the fully-developed thermal profile.

    Key thermal concepts for CFD:
    - Thermal boundary layer development
    - Nusselt number: ratio of convective to conductive heat transfer
    - Prandtl number: ratio of momentum to thermal diffusivity
    - Energy equation coupling with momentum equation

    The simplified profile used here assumes:
    - Constant wall temperature (Dirichlet BC)
    - Linear thermal development in flow direction

    Args:
        grid: vtkStructuredGrid to attach the field to
        t_wall: Wall temperature (constant)
        t_inlet: Inlet fluid temperature

    Returns:
        vtkFloatArray: Temperature scalar at each point
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    nx, ny, nz = dims[0], dims[1], dims[2]

    temperature = vtk.vtkFloatArray()
    temperature.SetName("Temperature")
    temperature.SetNumberOfComponents(1)

    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Development factor: 0 at inlet, 1 at outlet
                development = i / max(nx - 1, 1)

                # Distance from center (normalized 0-1)
                y_center = abs(2.0 * j / max(ny - 1, 1) - 1.0)

                # Blend from uniform inlet to developed profile
                # Developed profile: higher temp near walls
                t_developed = t_inlet + (t_wall - t_inlet) * y_center

                # Actual temperature: blend inlet to developed
                t = t_inlet + development * (t_developed - t_inlet)
                temperature.InsertNextValue(t)

    return temperature


def create_vorticity_field(grid, velocity_field):
    """
    Compute vorticity (curl of velocity) from the velocity field.

    Vorticity is a key quantity in fluid dynamics:
    - Measures local rotation of fluid elements
    - Zero in irrotational flow (potential flow)
    - Generated at boundaries and by body forces
    - Important for turbulence characterization

    In CFD, vorticity is used for:
    - Identifying vortex structures (Q-criterion, lambda2)
    - Vortex methods (alternative to solving NS equations)
    - Turbulence modeling (omega in k-omega models)

    For Poiseuille flow, vorticity is constant and in z-direction:
    omega_z = -du/dy = 2 * u_max * (2y/H) / H

    Args:
        grid: vtkStructuredGrid containing the mesh
        velocity_field: vtkFloatArray with velocity vectors

    Returns:
        vtkFloatArray: Vorticity vectors at each point
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    ny = dims[1]
    h = max(ny - 1, 1)

    vorticity = vtk.vtkFloatArray()
    vorticity.SetName("Vorticity")
    vorticity.SetNumberOfComponents(3)

    # For Poiseuille flow, du/dy = -2 * u_max * (2y/H) / H
    u_max = 1.5  # From velocity field creation

    for k in range(dims[2]):
        for j in range(dims[1]):
            for i in range(dims[0]):
                y_normalized = 2.0 * j / h - 1.0
                # Vorticity omega_z = -du/dy
                omega_z = 2.0 * u_max * 2.0 * y_normalized / h
                vorticity.InsertNextTuple3(0.0, 0.0, omega_z)

    return vorticity


def create_color_lookup_table(min_val, max_val, preset="rainbow"):
    """
    Create a color lookup table for scalar field visualization.

    Color maps are essential for scientific visualization:
    - Rainbow: Traditional but can obscure features
    - Diverging: Good for data with meaningful center (e.g., velocity deviation)
    - Sequential: Good for monotonic data (e.g., temperature)
    - Perceptually uniform: Better for quantitative reading

    Args:
        min_val: Minimum scalar value
        max_val: Maximum scalar value
        preset: Color map type ("rainbow", "cool_to_warm", "blue_to_red")

    Returns:
        vtkLookupTable configured for the scalar range
    """
    lut = vtk.vtkLookupTable()
    lut.SetRange(min_val, max_val)
    lut.SetNumberOfTableValues(256)

    if preset == "cool_to_warm":
        lut.SetHueRange(0.667, 0.0)  # Blue to Red
    elif preset == "blue_to_red":
        lut.SetHueRange(0.667, 0.0)
        lut.SetSaturationRange(1.0, 1.0)
    else:  # rainbow
        lut.SetHueRange(0.667, 0.0)

    lut.Build()
    return lut


def visualize_fields(grid):
    """
    Visualize the computational grid with multiple field representations.

    This demonstrates common CFD visualization techniques:
    1. Scalar color mapping (temperature, pressure)
    2. Vector glyphs (velocity arrows)
    3. Contour lines
    4. Streamlines (not shown here but important for vector fields)

    Args:
        grid: vtkStructuredGrid with attached fields
    """
    # Get field ranges for color mapping
    temp_range = grid.GetPointData().GetArray("Temperature").GetRange()
    vel_array = grid.GetPointData().GetArray("Velocity")

    # Compute velocity magnitude for coloring
    vel_magnitude = vtk.vtkFloatArray()
    vel_magnitude.SetName("VelocityMagnitude")
    vel_magnitude.SetNumberOfComponents(1)
    for i in range(vel_array.GetNumberOfTuples()):
        v = vel_array.GetTuple3(i)
        mag = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        vel_magnitude.InsertNextValue(mag)
    grid.GetPointData().AddArray(vel_magnitude)

    # Create surface representation colored by temperature
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(grid)
    surface_mapper.SetScalarModeToUsePointFieldData()
    surface_mapper.SelectColorArray("Temperature")
    surface_mapper.SetScalarRange(temp_range)

    lut = create_color_lookup_table(temp_range[0], temp_range[1], "cool_to_warm")
    surface_mapper.SetLookupTable(lut)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetEdgeVisibility(1)
    surface_actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
    surface_actor.GetProperty().SetOpacity(0.7)

    # Create velocity glyphs (arrows showing flow direction)
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetShaftResolution(16)

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
    glyph_mapper.SetScalarModeToUsePointFieldData()
    glyph_mapper.SelectColorArray("VelocityMagnitude")
    vel_range = vel_magnitude.GetRange()
    glyph_mapper.SetScalarRange(vel_range)

    glyph_actor = vtk.vtkActor()
    glyph_actor.SetMapper(glyph_mapper)

    # Create scalar bar (legend) for temperature
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("Temperature")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetPosition(0.85, 0.1)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.8)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(glyph_actor)
    renderer.AddActor(scalar_bar)
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Fields: CFD Flow Visualization")

    # Interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Camera positioning
    camera = renderer.GetActiveCamera()
    camera.SetPosition(10, -15, 10)
    camera.SetFocalPoint(5, 2, 1)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()

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


def print_field_info(grid):
    """
    Print educational information about the fields in the grid.

    Demonstrates how to query field properties in VTK, which is
    essential for debugging and validating CFD results.

    Args:
        grid: vtkStructuredGrid with attached fields
    """
    print("\n" + "=" * 70)
    print("VTK Fields: Educational Overview for CFD/FEA")
    print("=" * 70)

    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    print(f"\nGrid dimensions: {dims[0]} x {dims[1]} x {dims[2]}")
    print(f"Total points: {grid.GetNumberOfPoints()}")
    print(f"Total cells: {grid.GetNumberOfCells()}")

    print("\n" + "-" * 70)
    print("Point Data Fields (defined at mesh vertices):")
    print("-" * 70)

    point_data = grid.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(i)
        name = array.GetName()
        num_components = array.GetNumberOfComponents()
        num_tuples = array.GetNumberOfTuples()
        range_vals = array.GetRange()

        field_type = "Scalar" if num_components == 1 else f"Vector ({num_components}D)"
        print(f"\n  {name}:")
        print(f"    Type: {field_type}")
        print(f"    Number of values: {num_tuples}")
        if num_components == 1:
            print(f"    Range: [{range_vals[0]:.4f}, {range_vals[1]:.4f}]")

    print("\n" + "-" * 70)
    print("Cell Data Fields (defined at cell centers):")
    print("-" * 70)

    cell_data = grid.GetCellData()
    if cell_data.GetNumberOfArrays() == 0:
        print("\n  No cell data defined (point data only in this example)")

    print("\n" + "=" * 70)
    print("Key Concepts for CFD Field Data:")
    print("=" * 70)
    print("""
1. FIELD LOCATION: Choose point vs cell data based on:
   - Discretization scheme (FVM typically uses cell-centered)
   - Interpolation requirements (point data interpolates smoothly)
   - Solver requirements (many FEM solvers use node-centered)

2. FIELD TYPES for CFD:
   - Scalar: Pressure, temperature, density, turbulent kinetic energy
   - Vector: Velocity, vorticity, force, heat flux
   - Tensor: Stress tensor, strain rate tensor, Reynolds stress

3. VISUALIZATION TECHNIQUES:
   - Scalars: Color mapping, contour lines, isosurfaces
   - Vectors: Glyphs (arrows), streamlines, streaklines
   - Derived quantities: Vorticity, Q-criterion, lambda2

4. DATA INTERPOLATION:
   - VTK interpolates point data within cells for smooth visualization
   - Cell data appears as constant per cell (piecewise constant)
   - Higher-order elements can represent curved variations
""")


def main():
    """
    Main function demonstrating VTK fields for CFD applications.

    Creates a structured grid representing a channel flow domain and
    attaches multiple physical fields commonly found in CFD:
    - Velocity (vector field with Poiseuille profile)
    - Pressure (scalar field with linear gradient)
    - Temperature (scalar field with thermal development)
    - Vorticity (derived vector field)
    """
    # Create computational domain (channel flow)
    nx, ny, nz = 12, 6, 3  # Grid resolution
    grid = create_structured_grid(nx, ny, nz, spacing=1.0)

    # Attach physical fields
    velocity = create_velocity_field(grid, u_bulk=1.0)
    grid.GetPointData().AddArray(velocity)
    grid.GetPointData().SetActiveVectors("Velocity")

    pressure = create_pressure_field(grid, dp_dx=-0.1)
    grid.GetPointData().AddArray(pressure)

    temperature = create_temperature_field(grid, t_wall=100.0, t_inlet=20.0)
    grid.GetPointData().AddArray(temperature)

    vorticity = create_vorticity_field(grid, velocity)
    grid.GetPointData().AddArray(vorticity)

    # Print educational information
    print_field_info(grid)

    # Visualize the fields
    visualize_fields(grid)


if __name__ == "__main__":
    main()
