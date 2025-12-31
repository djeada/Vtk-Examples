"""
Potential Flow Around a Cylinder - Educational CFD Visualization

This module implements a 2D potential flow simulation around a circular cylinder
using analytical solutions from classical aerodynamics. This is a fundamental
problem in fluid mechanics that introduces key concepts used in computational
fluid dynamics (CFD) and aerodynamic analysis.

Physical Problem:
-----------------
Uniform flow past a circular cylinder is a classic problem in fluid mechanics.
The analytical solution combines:
1. Uniform freestream flow (constant velocity U_inf in x-direction)
2. A doublet at the origin (to enforce zero normal velocity at the surface)
3. Optional circulation (vortex at origin for lift generation)

The resulting flow demonstrates:
- Stagnation points where velocity is zero
- Flow acceleration around the cylinder sides
- Pressure variations according to Bernoulli's equation
- Lift generation when circulation is added (Kutta-Joukowski theorem)

Governing Equations:
--------------------
For irrotational, incompressible flow, we use potential flow theory:

1. Velocity Potential (φ):
   u = ∂φ/∂x,  v = ∂φ/∂y
   
2. Laplace's Equation:
   ∇²φ = ∂²φ/∂x² + ∂²φ/∂y² = 0

3. Stream Function (ψ) - lines of constant ψ are streamlines:
   u = ∂ψ/∂y,  v = -∂ψ/∂x

For flow past a cylinder of radius R with freestream velocity U_inf:

Velocity potential (polar coordinates):
   φ = U_inf * r * (1 + R²/r²) * cos(θ) - (Γ/2π) * θ

Stream function:
   ψ = U_inf * r * (1 - R²/r²) * sin(θ) + (Γ/2π) * ln(r/R)

where Γ is the circulation (positive counterclockwise).

Velocity Components (Cartesian):
   u = U_inf * (1 - R²*(x²-y²)/(x²+y²)² ) + Γ*y / (2π*(x²+y²))
   v = -2 * U_inf * R² * x*y / (x²+y²)² - Γ*x / (2π*(x²+y²))

Bernoulli's Equation:
--------------------
For steady, inviscid, incompressible flow along a streamline:
   p + ½ρV² = p_∞ + ½ρU_∞²

The pressure coefficient is defined as:
   Cp = (p - p_∞) / (½ρU_∞²) = 1 - (V/U_∞)²

Key CFD Concepts Demonstrated:
------------------------------
1. Potential Flow: Assumes inviscid (no viscosity) and irrotational flow
2. Superposition: Complex flows built from elementary solutions
3. Boundary Conditions: No-flow through solid surfaces
4. Pressure-Velocity Coupling: Bernoulli equation relates them
5. Circulation and Lift: Kutta-Joukowski theorem L = ρ * U_∞ * Γ

Limitations of Potential Flow:
------------------------------
- No viscous effects (no boundary layers, no drag from friction)
- No flow separation or wake behind the cylinder
- d'Alembert's paradox: zero drag for symmetric flow
- Real flows separate, creating a wake with pressure drag

This educational example helps understand:
- How CFD approximates continuous fields on discrete grids
- Velocity and pressure field computation
- Streamline visualization techniques
- The role of circulation in generating lift

References:
-----------
1. Anderson, J.D. (2016). "Fundamentals of Aerodynamics" 6th Ed., McGraw-Hill.
2. Kundu, P.K. & Cohen, I.M. (2015). "Fluid Mechanics" 6th Ed., Academic Press.
3. Katz, J. & Plotkin, A. (2001). "Low-Speed Aerodynamics" 2nd Ed., Cambridge.
"""

from typing import Tuple

import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk


# =============================================================================
# Physical and Simulation Parameters
# =============================================================================

# Cylinder parameters
CYLINDER_RADIUS: float = 1.0  # Radius of the cylinder (m)

# Flow parameters
U_INF: float = 1.0  # Freestream velocity (m/s)
RHO: float = 1.225  # Air density at sea level (kg/m³)
P_INF: float = 101325.0  # Freestream pressure (Pa)

# Circulation for lift (set to 0 for symmetric flow)
# For lift: Γ > 0 gives upward lift (counterclockwise circulation)
CIRCULATION: float = 2.0 * np.pi * CYLINDER_RADIUS * U_INF  # Gives significant lift

# Grid parameters
GRID_RESOLUTION: int = 50  # Number of grid points in each direction
DOMAIN_SIZE: float = 5.0  # Domain extends from -DOMAIN_SIZE to +DOMAIN_SIZE


# =============================================================================
# Core Potential Flow Functions
# =============================================================================

def compute_velocity_field(
    x: np.ndarray,
    y: np.ndarray,
    radius: float = CYLINDER_RADIUS,
    u_inf: float = U_INF,
    gamma: float = CIRCULATION
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the velocity field for potential flow around a cylinder.
    
    Uses the analytical solution for uniform flow + doublet + vortex.
    
    The velocity components are derived from the complex potential:
        w(z) = U_inf * (z + R²/z) + i*Γ/(2π) * ln(z)
    
    where z = x + iy is the complex coordinate.
    
    Args:
        x: 2D array of x-coordinates.
        y: 2D array of y-coordinates.
        radius: Cylinder radius (m).
        u_inf: Freestream velocity (m/s).
        gamma: Circulation strength (m²/s). Positive = counterclockwise.
    
    Returns:
        u: x-component of velocity field.
        v: y-component of velocity field.
    """
    # Distance squared from origin
    r_squared = x**2 + y**2
    
    # Avoid division by zero at the origin
    r_squared = np.maximum(r_squared, 1e-10)
    
    # Ratio R²/r² for doublet contribution
    R2_over_r2 = radius**2 / r_squared
    R2_over_r4 = radius**2 / (r_squared**2)
    
    # Velocity from uniform flow + doublet (non-lifting flow)
    # The formula u = U_inf * (1 - R²*(x²-y²)/r⁴) is equivalent to:
    # u = U_inf * (1 - R²/r² + 2*R²*y²/r⁴) since (x²-y²)/r⁴ = 1/r² - 2*y²/r⁴
    # The expanded form is used for numerical efficiency
    u = u_inf * (1.0 - R2_over_r2 + 2.0 * R2_over_r4 * y**2)
    v = -u_inf * 2.0 * R2_over_r4 * x * y
    
    # Add circulation (vortex at origin)
    # From stream function: u += Γ*y/(2π*r²), v -= Γ*x/(2π*r²)
    if abs(gamma) > 1e-10:
        u += gamma * y / (2.0 * np.pi * r_squared)
        v -= gamma * x / (2.0 * np.pi * r_squared)
    
    return u, v


def compute_stream_function(
    x: np.ndarray,
    y: np.ndarray,
    radius: float = CYLINDER_RADIUS,
    u_inf: float = U_INF,
    gamma: float = CIRCULATION
) -> np.ndarray:
    """
    Compute the stream function for potential flow around a cylinder.
    
    The stream function ψ is constant along streamlines, making it
    ideal for visualization. Contours of ψ are the streamlines.
    
    ψ = U_inf * y * (1 - R²/r²) + (Γ/2π) * ln(r/R)
    
    Args:
        x: 2D array of x-coordinates.
        y: 2D array of y-coordinates.
        radius: Cylinder radius (m).
        u_inf: Freestream velocity (m/s).
        gamma: Circulation strength (m²/s).
    
    Returns:
        psi: Stream function values.
    """
    r_squared = x**2 + y**2
    r = np.sqrt(np.maximum(r_squared, 1e-10))
    
    # Stream function for uniform flow + doublet
    psi = u_inf * y * (1.0 - radius**2 / r_squared)
    
    # Add circulation contribution
    if abs(gamma) > 1e-10:
        psi += gamma / (2.0 * np.pi) * np.log(r / radius)
    
    return psi


def compute_pressure_coefficient(
    u: np.ndarray,
    v: np.ndarray,
    u_inf: float = U_INF
) -> np.ndarray:
    """
    Compute the pressure coefficient using Bernoulli's equation.
    
    The pressure coefficient Cp relates local pressure to dynamic pressure:
        Cp = (p - p_∞) / (½ρU_∞²) = 1 - (V/U_∞)²
    
    This is a dimensionless measure of pressure:
    - Cp = 1 at stagnation points (V = 0)
    - Cp = 0 where V = U_∞
    - Cp < 0 where V > U_∞ (accelerated flow, suction)
    
    Args:
        u: x-velocity component.
        v: y-velocity component.
        u_inf: Freestream velocity (m/s).
    
    Returns:
        cp: Pressure coefficient field.
    """
    velocity_squared = u**2 + v**2
    cp = 1.0 - velocity_squared / (u_inf**2)
    return cp


def compute_pressure(
    cp: np.ndarray,
    rho: float = RHO,
    u_inf: float = U_INF,
    p_inf: float = P_INF
) -> np.ndarray:
    """
    Convert pressure coefficient to absolute pressure.
    
    From the definition of Cp:
        p = p_∞ + Cp * (½ρU_∞²)
    
    Args:
        cp: Pressure coefficient field.
        rho: Fluid density (kg/m³).
        u_inf: Freestream velocity (m/s).
        p_inf: Freestream pressure (Pa).
    
    Returns:
        p: Absolute pressure field (Pa).
    """
    dynamic_pressure = 0.5 * rho * u_inf**2
    return p_inf + cp * dynamic_pressure


def mask_inside_cylinder(
    field: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    radius: float = CYLINDER_RADIUS,
    mask_value: float = np.nan
) -> np.ndarray:
    """
    Mask field values inside the cylinder.
    
    The potential flow solution is only valid outside the cylinder.
    Points inside are masked to avoid displaying non-physical values.
    
    Args:
        field: Field to mask.
        x: x-coordinates.
        y: y-coordinates.
        radius: Cylinder radius.
        mask_value: Value to use for masked points.
    
    Returns:
        Masked field array.
    """
    inside = x**2 + y**2 < radius**2
    field_masked = field.copy()
    field_masked[inside] = mask_value
    return field_masked


# =============================================================================
# Geometry Creation
# =============================================================================

def create_cylinder_geometry(
    radius: float = CYLINDER_RADIUS,
    num_points: int = 100
) -> vtk.vtkPolyData:
    """
    Create a 2D cylinder (circle) geometry for visualization.
    
    The cylinder is represented as a closed polygon in the x-y plane.
    This serves as the solid body around which flow is computed.
    
    Args:
        radius: Cylinder radius (m).
        num_points: Number of points defining the circle.
    
    Returns:
        VTK PolyData representing the cylinder cross-section.
    """
    # Create circle points
    points = vtk.vtkPoints()
    for i in range(num_points):
        theta = 2.0 * np.pi * i / num_points
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        points.InsertNextPoint(x, y, 0.0)
    
    # Create polygon
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(num_points)
    for i in range(num_points):
        polygon.GetPointIds().SetId(i, i)
    
    # Create PolyData
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polygon)
    polydata.SetPolys(cells)
    
    return polydata


# =============================================================================
# VTK Data Structure Creation
# =============================================================================

def create_flow_field_grid(
    resolution: int = GRID_RESOLUTION,
    domain_size: float = DOMAIN_SIZE
) -> Tuple[vtk.vtkStructuredGrid, np.ndarray, np.ndarray]:
    """
    Create a VTK structured grid for the flow field visualization.
    
    The grid covers the domain around the cylinder. Points inside
    the cylinder will have their values masked in the field data.
    
    Args:
        resolution: Number of points in each direction.
        domain_size: Half-width of the square domain.
    
    Returns:
        grid: VTK structured grid.
        x: 2D array of x-coordinates.
        y: 2D array of y-coordinates.
    """
    # Create coordinate arrays
    x_1d = np.linspace(-domain_size, domain_size, resolution)
    y_1d = np.linspace(-domain_size, domain_size, resolution)
    x, y = np.meshgrid(x_1d, y_1d, indexing='ij')
    
    # Create VTK points
    points = vtk.vtkPoints()
    for i in range(resolution):
        for j in range(resolution):
            points.InsertNextPoint(x[i, j], y[i, j], 0.0)
    
    # Create structured grid
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(resolution, resolution, 1)
    grid.SetPoints(points)
    
    return grid, x, y


def add_velocity_to_grid(
    grid: vtk.vtkStructuredGrid,
    u: np.ndarray,
    v: np.ndarray
) -> None:
    """
    Add velocity vectors to the VTK grid as point data.
    
    The velocity is stored as a 3-component vector (u, v, 0) for
    compatibility with VTK's 3D vector operations.
    
    Args:
        grid: VTK structured grid.
        u: x-velocity component (2D array).
        v: y-velocity component (2D array).
    """
    # Replace NaN with zero for VTK (inside cylinder)
    u_clean = np.nan_to_num(u, nan=0.0)
    v_clean = np.nan_to_num(v, nan=0.0)
    
    # Create 3D velocity vectors (w=0 for 2D flow)
    w = np.zeros_like(u_clean)
    vectors = np.column_stack([u_clean.ravel(), v_clean.ravel(), w.ravel()])
    
    # Convert to VTK array
    vtk_vectors = numpy_to_vtk(vectors, deep=True, array_type=vtk.VTK_FLOAT)
    vtk_vectors.SetName("Velocity")
    vtk_vectors.SetNumberOfComponents(3)
    
    grid.GetPointData().SetVectors(vtk_vectors)


def add_scalar_to_grid(
    grid: vtk.vtkStructuredGrid,
    scalar: np.ndarray,
    name: str
) -> None:
    """
    Add a scalar field to the VTK grid as point data.
    
    Args:
        grid: VTK structured grid.
        scalar: 2D scalar array.
        name: Name for the scalar field.
    """
    scalar_clean = np.nan_to_num(scalar, nan=0.0)
    vtk_scalar = numpy_to_vtk(scalar_clean.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_scalar.SetName(name)
    grid.GetPointData().AddArray(vtk_scalar)
    
    # Set as active scalars if this is the first scalar
    if grid.GetPointData().GetScalars() is None:
        grid.GetPointData().SetScalars(vtk_scalar)


# =============================================================================
# Visualization
# =============================================================================

def create_streamline_actor(
    grid: vtk.vtkStructuredGrid,
    domain_size: float = DOMAIN_SIZE
) -> vtk.vtkActor:
    """
    Create streamlines visualization using VTK stream tracer.
    
    Streamlines are seeded from a line at the inlet (left boundary)
    and traced through the velocity field. They show the path that
    massless particles would follow in the flow.
    
    Args:
        grid: VTK grid containing velocity field.
        domain_size: Domain half-width for seed placement.
    
    Returns:
        VTK actor for streamline visualization.
    """
    # Create seed line at inlet (left side of domain)
    seeds = vtk.vtkLineSource()
    seeds.SetPoint1(-domain_size * 0.95, -domain_size * 0.8, 0.0)
    seeds.SetPoint2(-domain_size * 0.95, domain_size * 0.8, 0.0)
    seeds.SetResolution(25)  # Number of seed points
    
    # Stream tracer
    streamer = vtk.vtkStreamTracer()
    streamer.SetInputData(grid)
    streamer.SetSourceConnection(seeds.GetOutputPort())
    streamer.SetMaximumPropagation(domain_size * 4)
    streamer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
    streamer.SetInitialIntegrationStep(0.05)
    streamer.SetIntegrationDirectionToForward()
    streamer.SetComputeVorticity(True)
    
    # Create tubes from streamlines for better visibility
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(streamer.GetOutputPort())
    tube.SetRadius(0.03)
    tube.SetNumberOfSides(8)
    
    # Mapper - color by velocity magnitude
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("Velocity")
    mapper.SetScalarRange(0, U_INF * 2)
    
    # Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    return actor


def create_pressure_contour_actor(
    grid: vtk.vtkStructuredGrid,
    x: np.ndarray,
    y: np.ndarray,
    cp: np.ndarray
) -> vtk.vtkActor:
    """
    Create pressure coefficient contour visualization.
    
    Displays the pressure field as a colored surface. The color
    indicates the pressure coefficient:
    - Red/warm colors: high pressure (stagnation)
    - Blue/cool colors: low pressure (suction)
    
    Args:
        grid: VTK structured grid.
        x: x-coordinates.
        y: y-coordinates.
        cp: Pressure coefficient field.
    
    Returns:
        VTK actor for pressure visualization.
    """
    # Create a geometry filter to get renderable data
    geometry = vtk.vtkStructuredGridGeometryFilter()
    geometry.SetInputData(grid)
    
    # Create color lookup table
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)  # Blue to red
    lut.SetNumberOfColors(256)
    
    # Get Cp range (mask inside cylinder)
    cp_masked = mask_inside_cylinder(cp, x, y, CYLINDER_RADIUS, np.nan)
    cp_valid = cp_masked[~np.isnan(cp_masked)]
    cp_min, cp_max = np.min(cp_valid), np.max(cp_valid)
    lut.SetRange(cp_min, cp_max)
    lut.Build()
    
    # Mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometry.GetOutputPort())
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SelectColorArray("Cp")
    mapper.SetLookupTable(lut)
    mapper.SetScalarRange(cp_min, cp_max)
    
    # Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    return actor, lut


def create_cylinder_actor(radius: float = CYLINDER_RADIUS) -> vtk.vtkActor:
    """
    Create the cylinder (solid body) visualization.
    
    The cylinder is displayed as a filled gray circle to represent
    the solid boundary.
    
    Args:
        radius: Cylinder radius.
    
    Returns:
        VTK actor for cylinder.
    """
    # Create disk source (filled circle)
    disk = vtk.vtkDiskSource()
    disk.SetInnerRadius(0)
    disk.SetOuterRadius(radius)
    disk.SetCircumferentialResolution(64)
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(disk.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.5, 0.5, 0.5)  # Gray
    actor.SetPosition(0, 0, 0.01)  # Slightly in front
    
    return actor


def create_scalar_bar(
    lut: vtk.vtkLookupTable,
    title: str = "Pressure Coefficient (Cp)"
) -> vtk.vtkScalarBarActor:
    """
    Create a scalar bar (color legend) for the visualization.
    
    Args:
        lut: Lookup table used for coloring.
        title: Title for the scalar bar.
    
    Returns:
        VTK scalar bar actor.
    """
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle(title)
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetOrientationToVertical()
    scalar_bar.SetPosition(0.85, 0.1)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.8)
    
    return scalar_bar


def create_annotation_actor() -> vtk.vtkTextActor:
    """
    Create text annotation explaining the visualization.
    
    Returns:
        VTK text actor with simulation information.
    """
    text = vtk.vtkTextActor()
    
    # Calculate lift coefficient from circulation
    # L' = ρ * U_∞ * Γ (lift per unit span)
    # Cl = L' / (½ρU²*c) where c = 2R (diameter as chord)
    cl = CIRCULATION / (np.pi * CYLINDER_RADIUS * U_INF)
    
    info = (
        f"Potential Flow Around Cylinder\n"
        f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Cylinder Radius: {CYLINDER_RADIUS} m\n"
        f"Freestream Velocity: {U_INF} m/s\n"
        f"Circulation (Γ): {CIRCULATION:.2f} m²/s\n"
        f"Lift Coefficient: {cl:.2f}\n"
        f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        f"Streamlines show flow paths\n"
        f"Colors show pressure coefficient"
    )
    text.SetInput(info)
    text.GetTextProperty().SetFontSize(14)
    text.GetTextProperty().SetColor(1, 1, 1)
    text.GetTextProperty().SetFontFamilyToArial()
    text.SetPosition(10, 10)
    
    return text


def visualize_flow_field(
    grid: vtk.vtkStructuredGrid,
    x: np.ndarray,
    y: np.ndarray,
    cp: np.ndarray
) -> None:
    """
    Create and display the complete flow visualization.
    
    Combines:
    - Pressure coefficient contours (background)
    - Streamlines showing flow paths
    - Solid cylinder representation
    - Informative annotations
    
    Args:
        grid: VTK grid with velocity and scalar data.
        x: x-coordinates.
        y: y-coordinates.
        cp: Pressure coefficient field.
    """
    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.1, 0.1, 0.15)
    
    # Create and add actors
    cylinder_actor = create_cylinder_actor()
    renderer.AddActor(cylinder_actor)
    
    pressure_actor, lut = create_pressure_contour_actor(grid, x, y, cp)
    renderer.AddActor(pressure_actor)
    
    streamline_actor = create_streamline_actor(grid, DOMAIN_SIZE)
    renderer.AddActor(streamline_actor)
    
    scalar_bar = create_scalar_bar(lut)
    renderer.AddActor2D(scalar_bar)
    
    annotation = create_annotation_actor()
    renderer.AddActor2D(annotation)
    
    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 800)
    render_window.SetWindowName("Potential Flow Around Cylinder - CFD Visualization")
    
    # Create interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    
    # Set interaction style for 2D
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)
    
    # Position camera for 2D view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0, 0, 15)
    camera.SetFocalPoint(0, 0, 0)
    camera.SetViewUp(0, 1, 0)
    renderer.ResetCamera()
    
    # Start visualization
    render_window.Render()
    interactor.Start()
    
    # Cleanup
    render_window.Finalize()
    interactor.TerminateApp()


# =============================================================================
# Main Simulation
# =============================================================================

def run_flow_simulation() -> None:
    """
    Run the potential flow simulation and visualization.
    
    This function:
    1. Creates the computational grid
    2. Computes velocity field using potential flow theory
    3. Computes pressure coefficient from Bernoulli's equation
    4. Creates VTK data structures
    5. Launches interactive visualization
    
    Educational Outputs:
    - Demonstrates flow acceleration around cylinder sides
    - Shows stagnation points (front and rear for symmetric flow)
    - Illustrates pressure-velocity relationship
    - With circulation: shows asymmetric flow and lift generation
    """
    print("=" * 70)
    print("Potential Flow Around a Circular Cylinder")
    print("Educational CFD Visualization")
    print("=" * 70)
    print()
    print("Physical Parameters:")
    print(f"  Cylinder radius (R):     {CYLINDER_RADIUS} m")
    print(f"  Freestream velocity:     {U_INF} m/s")
    print(f"  Air density:             {RHO} kg/m³")
    print(f"  Circulation (Γ):         {CIRCULATION:.2f} m²/s")
    print()
    
    # Calculate and display derived quantities
    # Lift per unit span from Kutta-Joukowski theorem
    lift_per_span = RHO * U_INF * CIRCULATION
    print("Derived Quantities:")
    print(f"  Lift per unit span:      {lift_per_span:.2f} N/m")
    print(f"  (From Kutta-Joukowski: L' = ρ × U_∞ × Γ)")
    print()
    
    print(f"Grid: {GRID_RESOLUTION} × {GRID_RESOLUTION} points")
    print(f"Domain: [{-DOMAIN_SIZE}, {DOMAIN_SIZE}] × [{-DOMAIN_SIZE}, {DOMAIN_SIZE}]")
    print("=" * 70)
    print()
    
    # Step 1: Create computational grid
    print("Creating computational grid...")
    grid, x, y = create_flow_field_grid(GRID_RESOLUTION, DOMAIN_SIZE)
    
    # Step 2: Compute velocity field
    print("Computing velocity field from potential flow solution...")
    u, v = compute_velocity_field(x, y)
    
    # Mask values inside cylinder
    u = mask_inside_cylinder(u, x, y)
    v = mask_inside_cylinder(v, x, y)
    
    # Step 3: Compute pressure coefficient
    print("Computing pressure coefficient using Bernoulli's equation...")
    cp = compute_pressure_coefficient(u, v)
    cp = mask_inside_cylinder(cp, x, y)
    
    # Step 4: Compute stream function for additional analysis
    print("Computing stream function for streamline visualization...")
    psi = compute_stream_function(x, y)
    psi = mask_inside_cylinder(psi, x, y)
    
    # Step 5: Add data to VTK grid
    print("Building VTK data structures...")
    add_velocity_to_grid(grid, u, v)
    add_scalar_to_grid(grid, cp, "Cp")
    add_scalar_to_grid(grid, psi, "StreamFunction")
    
    # Compute velocity magnitude for visualization
    vel_mag = np.sqrt(np.nan_to_num(u, nan=0)**2 + np.nan_to_num(v, nan=0)**2)
    add_scalar_to_grid(grid, vel_mag, "VelocityMagnitude")
    
    print()
    print("Flow Field Statistics:")
    cp_valid = cp[~np.isnan(cp)]
    vel_valid = vel_mag[vel_mag > 0]
    print(f"  Pressure coefficient range: [{np.min(cp_valid):.3f}, {np.max(cp_valid):.3f}]")
    print(f"  Velocity magnitude range:   [{np.min(vel_valid):.3f}, {np.max(vel_valid):.3f}] m/s")
    print()
    
    # Find stagnation points (where velocity is minimum, near cylinder surface)
    r = np.sqrt(x**2 + y**2)
    near_surface = (r > CYLINDER_RADIUS * 0.95) & (r < CYLINDER_RADIUS * 1.1)
    vel_near_surface = vel_mag.copy()
    vel_near_surface[~near_surface] = np.inf
    vel_near_surface[np.isnan(vel_near_surface)] = np.inf
    if np.any(np.isfinite(vel_near_surface)):
        min_idx = np.unravel_index(np.argmin(vel_near_surface), vel_near_surface.shape)
        stag_x, stag_y = x[min_idx], y[min_idx]
        stag_theta = np.degrees(np.arctan2(stag_y, stag_x))
        print(f"  Stagnation point found near: ({stag_x:.2f}, {stag_y:.2f})")
        print(f"  Stagnation angle: {stag_theta:.1f}°")
    print()
    
    # Step 6: Launch visualization
    print("Launching VTK visualization...")
    print("  - Drag to rotate view")
    print("  - Scroll to zoom")
    print("  - Close window to exit")
    print("=" * 70)
    
    visualize_flow_field(grid, x, y, cp)
    
    print()
    print("Simulation complete.")


def main():
    """
    Main entry point for the potential flow simulation.
    
    This educational example demonstrates:
    1. Potential flow theory and its analytical solutions
    2. Pressure-velocity coupling through Bernoulli's equation
    3. Lift generation through circulation (Kutta-Joukowski theorem)
    4. VTK-based scientific visualization techniques
    
    The visualization shows:
    - Streamlines that flow around the cylinder
    - Pressure coefficient distribution (colored surface)
    - Stagnation points where flow velocity is zero
    - Asymmetric flow pattern when circulation is present
    """
    run_flow_simulation()


if __name__ == "__main__":
    main()
