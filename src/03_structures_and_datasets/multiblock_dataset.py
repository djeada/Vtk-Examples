"""
VTK Multiblock Datasets: Composite Data Structures for Complex CFD Domains

This module provides an educational demonstration of VTK multiblock datasets,
which are essential for managing complex computational domains in CFD and FEA.

What are Multiblock Datasets?
-----------------------------
A vtkMultiBlockDataSet is a composite data structure that can contain multiple
datasets of different types. Think of it as a container that groups related
datasets together while maintaining their individual identities.

Why Use Multiblock Datasets?
----------------------------
In CFD/FEA, real-world simulations often require:

1. DOMAIN DECOMPOSITION:
   - Large domains split into smaller blocks for parallel processing
   - Each processor handles one or more blocks
   - Multiblock structure maintains the global view

2. MULTI-REGION MESHES:
   - Different mesh types for different regions (structured + unstructured)
   - Boundary layer meshes near walls (structured, fine)
   - Core flow region (unstructured, coarser)
   - Overset/chimera meshes for moving components

3. MULTI-PHYSICS COUPLING:
   - Solid regions (heat conduction in walls)
   - Fluid regions (convective heat transfer)
   - Interface between domains

4. ASSEMBLY/COMPONENT MODELING:
   - Complex geometries broken into components
   - Engine with pistons, valves, cylinders as separate blocks
   - Aircraft with fuselage, wings, engines

CFD/FEA Context:
----------------
Multiblock datasets mirror how real CFD software organizes data:

- OpenFOAM: Uses regions for multi-physics
- ANSYS Fluent: Cell zones and partitions
- Star-CCM+: Regions and interfaces
- Paraview: Native support for multiblock visualization

This example demonstrates creating a multiblock dataset representing
a simplified pipe flow domain with different mesh zones commonly
found in industrial CFD simulations.
"""

import vtk


# Default base temperatures for each zone (Kelvin)
ZONE_BASE_TEMPERATURES = [300.0, 350.0, 400.0, 500.0]

# Temperature gradient along flow direction (K per unit length)
TEMP_GRADIENT_ALONG_FLOW = 5.0


def create_pipe_inlet_zone(radius=1.0, length=2.0, resolution=20):
    """
    Create the inlet zone of a pipe flow domain.

    The inlet zone in CFD simulations typically requires:
    - Boundary condition specification (velocity/pressure inlet)
    - Possible turbulence specification (k, omega, or intensity)
    - May have refined mesh for inlet profile development

    This uses a cylinder to represent the inlet section of pipe flow.

    Args:
        radius: Pipe radius
        length: Zone length in flow direction
        resolution: Mesh resolution (number of facets)

    Returns:
        vtkPolyData: Inlet zone geometry with metadata
    """
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(length)
    cylinder.SetResolution(resolution)
    cylinder.SetCenter(0, length / 2, 0)  # Position at origin
    cylinder.Update()

    # Add zone identifier as field data
    zone_id = vtk.vtkIntArray()
    zone_id.SetName("ZoneID")
    zone_id.SetNumberOfTuples(1)
    zone_id.SetValue(0, 1)  # Zone 1 = Inlet

    zone_name = vtk.vtkStringArray()
    zone_name.SetName("ZoneName")
    zone_name.SetNumberOfTuples(1)
    zone_name.SetValue(0, "Inlet")

    output = vtk.vtkPolyData()
    output.ShallowCopy(cylinder.GetOutput())
    output.GetFieldData().AddArray(zone_id)
    output.GetFieldData().AddArray(zone_name)

    return output


def create_pipe_core_zone(radius=1.0, length=6.0, resolution=20, x_offset=2.0):
    """
    Create the core flow zone of a pipe flow domain.

    The core/main zone is where:
    - Flow develops from inlet conditions
    - Most of the computational cells reside
    - Key physics (turbulence, heat transfer) are resolved
    - Mesh may be coarser than near-wall regions

    In real CFD, this zone often uses:
    - Hexahedral cells for structured pipe flows
    - Polyhedral cells for complex geometries
    - Appropriate refinement for capturing flow features

    Args:
        radius: Pipe radius
        length: Zone length in flow direction
        resolution: Mesh resolution
        x_offset: Position offset from inlet

    Returns:
        vtkPolyData: Core zone geometry with metadata
    """
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(length)
    cylinder.SetResolution(resolution)
    cylinder.SetCenter(0, x_offset + length / 2, 0)
    cylinder.Update()

    # Add zone identifier
    zone_id = vtk.vtkIntArray()
    zone_id.SetName("ZoneID")
    zone_id.SetNumberOfTuples(1)
    zone_id.SetValue(0, 2)  # Zone 2 = Core

    zone_name = vtk.vtkStringArray()
    zone_name.SetName("ZoneName")
    zone_name.SetNumberOfTuples(1)
    zone_name.SetValue(0, "Core")

    output = vtk.vtkPolyData()
    output.ShallowCopy(cylinder.GetOutput())
    output.GetFieldData().AddArray(zone_id)
    output.GetFieldData().AddArray(zone_name)

    return output


def create_pipe_outlet_zone(radius=1.0, length=2.0, resolution=20, x_offset=8.0):
    """
    Create the outlet zone of a pipe flow domain.

    The outlet zone in CFD simulations:
    - Specifies outflow boundary conditions
    - Common types: pressure outlet, outflow, convective outlet
    - May need careful treatment to avoid reflections
    - Often uses coarser mesh than inlet/core

    Outlet boundary conditions are critical for:
    - Preventing backflow (negative velocities)
    - Maintaining mass conservation
    - Avoiding numerical instabilities

    Args:
        radius: Pipe radius
        length: Zone length in flow direction
        resolution: Mesh resolution
        x_offset: Position offset from inlet

    Returns:
        vtkPolyData: Outlet zone geometry with metadata
    """
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(length)
    cylinder.SetResolution(resolution)
    cylinder.SetCenter(0, x_offset + length / 2, 0)
    cylinder.Update()

    # Add zone identifier
    zone_id = vtk.vtkIntArray()
    zone_id.SetName("ZoneID")
    zone_id.SetNumberOfTuples(1)
    zone_id.SetValue(0, 3)  # Zone 3 = Outlet

    zone_name = vtk.vtkStringArray()
    zone_name.SetName("ZoneName")
    zone_name.SetNumberOfTuples(1)
    zone_name.SetValue(0, "Outlet")

    output = vtk.vtkPolyData()
    output.ShallowCopy(cylinder.GetOutput())
    output.GetFieldData().AddArray(zone_id)
    output.GetFieldData().AddArray(zone_name)

    return output


def create_obstacle(radius=0.3, center=(0, 5, 0)):
    """
    Create an obstacle within the flow domain.

    Obstacles in CFD simulations:
    - Create flow separation and recirculation
    - Generate vortex shedding (von Karman street)
    - Require refined mesh in wake region
    - Test turbulence model performance

    Common obstacle shapes in CFD validation:
    - Sphere: 3D wake, axisymmetric base case
    - Cylinder: 2D wake, vortex shedding studies
    - Cube: Sharp edge separation, building aerodynamics

    Args:
        radius: Obstacle radius
        center: Position in domain

    Returns:
        vtkPolyData: Obstacle geometry with metadata
    """
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetCenter(center)
    sphere.SetThetaResolution(20)
    sphere.SetPhiResolution(20)
    sphere.Update()

    # Add zone identifier
    zone_id = vtk.vtkIntArray()
    zone_id.SetName("ZoneID")
    zone_id.SetNumberOfTuples(1)
    zone_id.SetValue(0, 4)  # Zone 4 = Obstacle

    zone_name = vtk.vtkStringArray()
    zone_name.SetName("ZoneName")
    zone_name.SetNumberOfTuples(1)
    zone_name.SetValue(0, "Obstacle")

    output = vtk.vtkPolyData()
    output.ShallowCopy(sphere.GetOutput())
    output.GetFieldData().AddArray(zone_id)
    output.GetFieldData().AddArray(zone_name)

    return output


def create_multiblock_dataset(blocks, block_names=None):
    """
    Assemble individual blocks into a multiblock dataset.

    This function demonstrates the key multiblock operations:
    1. Creating the container dataset
    2. Setting the number of blocks
    3. Assigning each block with optional metadata

    Multiblock metadata is important for:
    - Identifying zones in post-processing
    - Setting up boundary conditions in preprocessing
    - Parallel load balancing

    Args:
        blocks: List of vtkDataSet objects (blocks)
        block_names: Optional list of names for each block

    Returns:
        vtkMultiBlockDataSet containing all blocks
    """
    multiblock = vtk.vtkMultiBlockDataSet()
    multiblock.SetNumberOfBlocks(len(blocks))

    for i, block in enumerate(blocks):
        multiblock.SetBlock(i, block)
        if block_names and i < len(block_names):
            multiblock.GetMetaData(i).Set(vtk.vtkCompositeDataSet.NAME(), block_names[i])

    return multiblock


def add_block_scalar_field(multiblock, field_name="Temperature"):
    """
    Add a scalar field to each block in the multiblock dataset.

    In CFD, different zones may have different field values:
    - Inlet: prescribed boundary values
    - Core: computed solution values
    - Outlet: extrapolated or computed values
    - Solid zones: different physics (conduction vs convection)

    This demonstrates how to iterate over blocks and add fields.

    Args:
        multiblock: vtkMultiBlockDataSet to modify
        field_name: Name of the scalar field

    Returns:
        None (modifies multiblock in place)
    """
    iterator = multiblock.NewIterator()
    iterator.InitTraversal()

    block_idx = 0
    while not iterator.IsDoneWithTraversal():
        block = iterator.GetCurrentDataObject()
        if block:
            # Create temperature field based on zone
            temp_array = vtk.vtkFloatArray()
            temp_array.SetName(field_name)
            temp_array.SetNumberOfComponents(1)

            base_temp = ZONE_BASE_TEMPERATURES[block_idx % len(ZONE_BASE_TEMPERATURES)]

            # Add some variation based on position
            for i in range(block.GetNumberOfPoints()):
                point = block.GetPoint(i)
                # Temperature increases along flow direction (y)
                temp = base_temp + point[1] * TEMP_GRADIENT_ALONG_FLOW
                temp_array.InsertNextValue(temp)

            block.GetPointData().AddArray(temp_array)

        iterator.GoToNextItem()
        block_idx += 1


def visualize_multiblock(multiblock):
    """
    Visualize the multiblock dataset with color-coded zones.

    Demonstrates multiblock visualization techniques:
    1. Using vtkCompositePolyDataMapper for efficient rendering
    2. Color-coding different blocks/zones
    3. Adding block visibility controls

    In production CFD tools, similar techniques are used to:
    - Show/hide specific zones
    - Apply different rendering to different regions
    - Highlight interfaces between zones

    Args:
        multiblock: vtkMultiBlockDataSet to visualize
    """
    # Convert to polydata for rendering
    geometry_filter = vtk.vtkCompositeDataGeometryFilter()
    geometry_filter.SetInputDataObject(multiblock)
    geometry_filter.Update()

    # Create color map for temperature
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)  # Blue to Red
    lut.SetNumberOfTableValues(256)
    lut.Build()

    # Main mapper with scalar coloring
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometry_filter.GetOutputPort())
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("Temperature")
    mapper.SetScalarRange(300, 500)
    mapper.SetLookupTable(lut)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetEdgeVisibility(1)
    actor.GetProperty().SetEdgeColor(0.1, 0.1, 0.1)
    actor.GetProperty().SetLineWidth(1)

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("Temperature (K)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetPosition(0.85, 0.1)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.8)

    # Create legend for zones
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(4)
    legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPositionCoordinate().SetValue(0.02, 0.7)
    legend.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPosition2Coordinate().SetValue(0.2, 0.95)
    legend.UseBackgroundOn()
    legend.SetBackgroundColor(0.1, 0.1, 0.1)

    zone_colors = [
        (0.2, 0.4, 0.8),   # Inlet - Blue
        (0.2, 0.8, 0.2),   # Core - Green
        (0.8, 0.4, 0.2),   # Outlet - Orange
        (0.8, 0.2, 0.2),   # Obstacle - Red
    ]
    zone_names = ["Inlet Zone", "Core Zone", "Outlet Zone", "Obstacle"]

    sphere = vtk.vtkSphereSource()
    sphere.Update()
    sphere_output = sphere.GetOutput()

    for i, (name, color) in enumerate(zip(zone_names, zone_colors)):
        legend.SetEntry(i, sphere_output, name, color)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.AddActor(scalar_bar)
    renderer.AddActor(legend)
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Multiblock: CFD Domain Decomposition")

    # Interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Camera positioning for side view of pipe
    camera = renderer.GetActiveCamera()
    camera.SetPosition(15, 5, 8)
    camera.SetFocalPoint(0, 5, 0)
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


def print_multiblock_info(multiblock):
    """
    Print educational information about the multiblock structure.

    Demonstrates how to traverse and query multiblock datasets,
    which is essential for:
    - Debugging mesh issues
    - Extracting zone-specific data
    - Setting up boundary conditions

    Args:
        multiblock: vtkMultiBlockDataSet to analyze
    """
    print("\n" + "=" * 70)
    print("VTK Multiblock Dataset: Educational Overview for CFD")
    print("=" * 70)

    print(f"\nNumber of blocks: {multiblock.GetNumberOfBlocks()}")

    print("\n" + "-" * 70)
    print("Block Details:")
    print("-" * 70)

    total_points = 0
    total_cells = 0

    for i in range(multiblock.GetNumberOfBlocks()):
        block = multiblock.GetBlock(i)
        if block:
            metadata = multiblock.GetMetaData(i)
            name = "Unnamed"
            if metadata:
                name = metadata.Get(vtk.vtkCompositeDataSet.NAME()) or f"Block {i}"

            num_points = block.GetNumberOfPoints()
            num_cells = block.GetNumberOfCells()
            total_points += num_points
            total_cells += num_cells

            # Get zone info from field data if available
            zone_name_array = block.GetFieldData().GetAbstractArray("ZoneName")
            zone_name = zone_name_array.GetValue(0) if zone_name_array else "Unknown"

            print(f"\n  Block {i}: {name}")
            print(f"    Zone Name: {zone_name}")
            print(f"    Points: {num_points}")
            print(f"    Cells: {num_cells}")

            # List point data arrays
            pd = block.GetPointData()
            if pd.GetNumberOfArrays() > 0:
                print(f"    Point Data Arrays:")
                for j in range(pd.GetNumberOfArrays()):
                    arr = pd.GetArray(j)
                    if arr:
                        print(f"      - {arr.GetName()}: {arr.GetNumberOfTuples()} values")

    print(f"\n  TOTAL: {total_points} points, {total_cells} cells")

    print("\n" + "=" * 70)
    print("Key Concepts for CFD Multiblock Datasets:")
    print("=" * 70)
    print("""
1. DOMAIN DECOMPOSITION:
   - Multiblock enables parallel CFD by splitting domain across processors
   - Each block can be solved independently with interface communication
   - Load balancing ensures equal work distribution

2. ZONE TYPES in CFD:
   - Fluid zones: Where flow equations are solved
   - Solid zones: Conduction only (no flow)
   - Porous zones: Modified momentum equations
   - Interface zones: Coupling between regions

3. BOUNDARY CONDITIONS:
   - Each block boundary needs appropriate conditions
   - Internal interfaces: Continuity/conservation
   - External boundaries: Inlet, outlet, wall, symmetry

4. MESH HANDLING:
   - Non-conformal interfaces: Different mesh densities
   - Sliding interfaces: Rotating machinery
   - Overset/chimera: Moving components

5. POST-PROCESSING:
   - Extract individual zones for analysis
   - Compute zone-averaged quantities
   - Visualize interface data transfer
""")


def main():
    """
    Main function demonstrating VTK multiblock datasets for CFD.

    Creates a simplified pipe flow domain with multiple zones:
    - Inlet zone: Entry region with specified velocity
    - Core zone: Main flow development region
    - Outlet zone: Exit with pressure specification
    - Obstacle: Bluff body for wake generation
    """
    # Create individual zones
    inlet = create_pipe_inlet_zone(radius=1.0, length=2.0)
    core = create_pipe_core_zone(radius=1.0, length=6.0, x_offset=2.0)
    outlet = create_pipe_outlet_zone(radius=1.0, length=2.0, x_offset=8.0)
    obstacle = create_obstacle(radius=0.3, center=(0, 5, 0))

    # Assemble into multiblock
    blocks = [inlet, core, outlet, obstacle]
    block_names = ["Inlet", "Core", "Outlet", "Obstacle"]
    multiblock = create_multiblock_dataset(blocks, block_names)

    # Add field data to each block
    add_block_scalar_field(multiblock, "Temperature")

    # Print educational information
    print_multiblock_info(multiblock)

    # Visualize
    visualize_multiblock(multiblock)


if __name__ == "__main__":
    main()
