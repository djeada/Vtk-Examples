"""
VTK Structured Mesh: CFD Educational Example with Heat Conduction Solver

This module provides an educational demonstration of structured meshes for
Computational Fluid Dynamics (CFD) applications. It combines mesh generation,
a numerical solver for steady-state heat conduction, and VTK visualization.

What is a Structured Mesh?
--------------------------
A structured mesh is a grid where:
- Points (nodes) are arranged in a regular i-j-k topology
- Each interior node has the same number of neighbors
- Connectivity is IMPLICIT (determined by indices, not stored explicitly)
- The mesh can be uniform or non-uniform in spacing

Mesh Components:
----------------
1. NODES (Points): Locations where solution variables are computed
   - In CFD: temperature, pressure, velocity are stored at nodes
   - Identified by indices (i, j) in 2D or (i, j, k) in 3D

2. EDGES: Connections between adjacent nodes
   - Form the "skeleton" of the mesh
   - Used for computing gradients and fluxes

3. CELLS (Elements): Smallest units bounded by edges
   - In 2D: quadrilaterals bounded by 4 nodes
   - In 3D: hexahedra bounded by 8 nodes
   - Control volumes for finite volume methods

Physical Problem: Steady-State Heat Conduction
----------------------------------------------
This example solves the 2D Laplace equation for heat conduction:

    ∂²T/∂x² + ∂²T/∂y² = 0

Boundary Conditions (Dirichlet):
- Left wall:   T = T_hot (100°C)
- Right wall:  T = T_cold (0°C)
- Top wall:    T = T_cold (0°C)
- Bottom wall: T = T_cold (0°C)

This represents a heated plate with one hot edge, simulating heat spreading
from a heat source through a conductive material.

Numerical Method: Finite Difference with Gauss-Seidel
-----------------------------------------------------
The Laplace equation is discretized using central differences:

    (T[i+1,j] - 2*T[i,j] + T[i-1,j])/dx² +
    (T[i,j+1] - 2*T[i,j] + T[i,j-1])/dy² = 0

For a uniform grid (dx = dy), this simplifies to:

    T[i,j] = (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1]) / 4

The Gauss-Seidel iterative method updates each node using the latest
available values from neighboring nodes, converging to the steady-state
solution.

CFD Applications:
-----------------
Structured meshes with similar solvers are used for:
- Heat exchangers and thermal management
- Electronic cooling simulations
- Building thermal analysis
- Geothermal heat flow modeling

Visualization:
--------------
The VTK visualization shows:
1. The structured mesh with nodes and edges
2. Temperature field as a color map (blue=cold, red=hot)
3. Interactive controls for mesh visibility
"""

import numpy as np
import vtk


def create_structured_mesh(
    nx: int, ny: int, lx: float, ly: float
) -> vtk.vtkStructuredGrid:
    """
    Create a 2D structured mesh for CFD simulation.

    The mesh is created with uniform spacing in both directions. Nodes are
    arranged in a regular i-j pattern where:
    - i increases in the x-direction (0 to nx-1)
    - j increases in the y-direction (0 to ny-1)

    Args:
        nx: Number of nodes in x-direction
        ny: Number of nodes in y-direction
        lx: Domain length in x-direction (meters)
        ly: Domain length in y-direction (meters)

    Returns:
        vtkStructuredGrid: The structured mesh with nodes positioned
    """
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(nx, ny, 1)  # 2D mesh (nz=1)

    # Calculate uniform spacing
    dx = lx / (nx - 1) if nx > 1 else lx
    dy = ly / (ny - 1) if ny > 1 else ly

    # Create points array
    points = vtk.vtkPoints()
    for j in range(ny):
        for i in range(nx):
            x = i * dx
            y = j * dy
            z = 0.0  # 2D mesh, z=0
            points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)
    return grid


def initialize_temperature_field(nx: int, ny: int, T_init: float = 0.0) -> np.ndarray:
    """
    Initialize the temperature field for the simulation.

    Creates a 2D numpy array representing temperature at each mesh node.
    The array uses row-major ordering where T[j, i] corresponds to the
    node at position (i, j) in the mesh.

    Args:
        nx: Number of nodes in x-direction
        ny: Number of nodes in y-direction
        T_init: Initial temperature value (default: 0.0)

    Returns:
        np.ndarray: 2D array of shape (ny, nx) with initial temperatures
    """
    return np.full((ny, nx), T_init, dtype=np.float64)


def apply_boundary_conditions(T: np.ndarray, T_hot: float, T_cold: float) -> None:
    """
    Apply Dirichlet boundary conditions to the temperature field.

    Boundary conditions define the temperature at domain boundaries:
    - Left wall (i=0): Hot temperature (heat source)
    - Right wall (i=nx-1): Cold temperature
    - Bottom wall (j=0): Cold temperature
    - Top wall (j=ny-1): Cold temperature

    In CFD, these represent:
    - Isothermal walls with fixed temperatures
    - Contact with heat source/sink

    Args:
        T: Temperature field array (modified in-place)
        T_hot: Temperature at left boundary (hot wall)
        T_cold: Temperature at other boundaries (cold walls)
    """
    # Left wall: hot boundary (heat source)
    T[:, 0] = T_hot

    # Right wall: cold boundary
    T[:, -1] = T_cold

    # Bottom wall: cold boundary
    T[0, :] = T_cold

    # Top wall: cold boundary
    T[-1, :] = T_cold

    # Corner nodes: average of adjacent boundaries
    T[0, 0] = (T_hot + T_cold) / 2
    T[-1, 0] = (T_hot + T_cold) / 2


def solve_laplace_gauss_seidel(
    T: np.ndarray,
    T_hot: float,
    T_cold: float,
    tolerance: float = 1e-6,
    max_iterations: int = 10000,
) -> tuple[np.ndarray, int, float]:
    """
    Solve the 2D Laplace equation using the Gauss-Seidel iterative method.

    The Laplace equation (∂²T/∂x² + ∂²T/∂y² = 0) is discretized using
    central finite differences on a uniform grid:

        T[i,j] = (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1]) / 4

    The Gauss-Seidel method updates each interior node using the most
    recent values from neighboring nodes. This in-place update provides
    faster convergence than the Jacobi method.

    Convergence criterion: Maximum absolute change in temperature between
    iterations falls below the specified tolerance.

    Args:
        T: Temperature field (modified in-place during iteration)
        T_hot: Hot boundary temperature (applied each iteration)
        T_cold: Cold boundary temperature (applied each iteration)
        tolerance: Convergence tolerance for temperature change
        max_iterations: Maximum number of iterations allowed

    Returns:
        tuple: (T, iterations, residual)
            - T: Final temperature field
            - iterations: Number of iterations performed
            - residual: Final maximum temperature change
    """
    ny, nx = T.shape
    residual = float("inf")

    for iteration in range(max_iterations):
        max_change = 0.0

        # Update interior nodes using Gauss-Seidel
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                T_old = T[j, i]
                # Five-point stencil average
                T[j, i] = 0.25 * (T[j, i + 1] + T[j, i - 1] + T[j + 1, i] + T[j - 1, i])
                change = abs(T[j, i] - T_old)
                if change > max_change:
                    max_change = change

        # Re-apply boundary conditions
        apply_boundary_conditions(T, T_hot, T_cold)

        residual = max_change

        # Check convergence
        if residual < tolerance:
            return T, iteration + 1, residual

    return T, max_iterations, residual


def add_temperature_to_mesh(grid: vtk.vtkStructuredGrid, T: np.ndarray) -> None:
    """
    Add the computed temperature field to the VTK structured mesh.

    Temperature values are stored as point data (vertex-centered data)
    on the structured grid. This is the standard approach for scalar
    fields in CFD visualization.

    Args:
        grid: VTK structured grid to add temperature data to
        T: 2D temperature array from the solver
    """
    ny, nx = T.shape
    temperature_array = vtk.vtkFloatArray()
    temperature_array.SetName("Temperature")
    temperature_array.SetNumberOfComponents(1)

    # VTK structured grids use i-fastest ordering (same as row-major flattening)
    for j in range(ny):
        for i in range(nx):
            temperature_array.InsertNextValue(T[j, i])

    grid.GetPointData().AddArray(temperature_array)
    grid.GetPointData().SetActiveScalars("Temperature")


def compute_mesh_quality(grid: vtk.vtkStructuredGrid) -> dict:
    """
    Compute quality metrics for the structured mesh.

    Grid quality is critical for CFD accuracy. This function calculates:
    - Aspect ratio: ratio of cell dimensions (ideally close to 1)
    - Cell size uniformity: variation in cell sizes

    For structured meshes, quality is typically excellent since cells
    are regular quadrilaterals or hexahedra.

    Args:
        grid: VTK structured grid to analyze

    Returns:
        dict: Dictionary containing quality metrics
    """
    dims = [0, 0, 0]
    grid.GetDimensions(dims)
    nx, ny, nz = dims

    bounds = grid.GetBounds()
    lx = bounds[1] - bounds[0]
    ly = bounds[3] - bounds[2]

    dx = lx / (nx - 1) if nx > 1 else lx
    dy = ly / (ny - 1) if ny > 1 else ly

    aspect_ratio = max(dx, dy) / min(dx, dy) if min(dx, dy) > 0 else 1.0

    return {
        "dimensions": (nx, ny, nz),
        "num_nodes": grid.GetNumberOfPoints(),
        "num_cells": grid.GetNumberOfCells(),
        "dx": dx,
        "dy": dy,
        "aspect_ratio": aspect_ratio,
        "domain_size": (lx, ly),
    }


def create_mesh_visualization_actors(
    grid: vtk.vtkStructuredGrid,
) -> tuple[vtk.vtkActor, vtk.vtkActor, vtk.vtkActor, vtk.vtkLookupTable]:
    """
    Create VTK actors for visualizing mesh structure.

    Creates visualization components:
    1. Surface actor: Colored by temperature field
    2. Edge actor: Shows mesh wireframe
    3. Node actor: Shows mesh nodes as points
    4. Lookup table: Color mapping for temperature values

    Args:
        grid: VTK structured grid with temperature data

    Returns:
        tuple: (surface_actor, edge_actor, node_actor, lut)
            - surface_actor: VTK actor showing temperature color map
            - edge_actor: VTK actor showing mesh edges
            - node_actor: VTK actor showing mesh nodes
            - lut: VTK lookup table for color mapping
    """
    # Create color lookup table (blue=cold, red=hot)
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)  # Blue to red
    lut.SetNumberOfTableValues(256)
    lut.Build()

    # Surface actor with temperature coloring
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(grid)
    surface_mapper.SetScalarModeToUsePointFieldData()
    surface_mapper.SelectColorArray("Temperature")
    temp_range = grid.GetPointData().GetArray("Temperature").GetRange()
    surface_mapper.SetScalarRange(temp_range)
    surface_mapper.SetLookupTable(lut)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetOpacity(0.9)

    # Edge actor (wireframe)
    edge_filter = vtk.vtkExtractEdges()
    edge_filter.SetInputData(grid)

    edge_mapper = vtk.vtkPolyDataMapper()
    edge_mapper.SetInputConnection(edge_filter.GetOutputPort())
    edge_mapper.ScalarVisibilityOff()

    edge_actor = vtk.vtkActor()
    edge_actor.SetMapper(edge_mapper)
    edge_actor.GetProperty().SetColor(0.2, 0.2, 0.2)
    edge_actor.GetProperty().SetLineWidth(1)

    # Node actor (points)
    point_filter = vtk.vtkVertexGlyphFilter()
    point_filter.SetInputData(grid)

    point_mapper = vtk.vtkPolyDataMapper()
    point_mapper.SetInputConnection(point_filter.GetOutputPort())
    point_mapper.ScalarVisibilityOff()

    node_actor = vtk.vtkActor()
    node_actor.SetMapper(point_mapper)
    node_actor.GetProperty().SetColor(0.0, 0.0, 0.0)
    node_actor.GetProperty().SetPointSize(4)

    return surface_actor, edge_actor, node_actor, lut


def create_scalar_bar(lut: vtk.vtkLookupTable) -> vtk.vtkScalarBarActor:
    """
    Create a color bar legend for the temperature field.

    Args:
        lut: VTK lookup table used for temperature coloring

    Returns:
        vtkScalarBarActor: Color bar actor for the visualization
    """
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(lut)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetPosition(0.85, 0.1)
    scalar_bar.SetWidth(0.1)
    scalar_bar.SetHeight(0.8)
    return scalar_bar


def create_text_annotation(text: str, position: tuple) -> vtk.vtkTextActor:
    """
    Create a text annotation for the visualization.

    Args:
        text: Text string to display
        position: (x, y) position in normalized viewport coordinates

    Returns:
        vtkTextActor: Text actor for the visualization
    """
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(text)
    text_actor.GetTextProperty().SetFontSize(14)
    text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    text_actor.SetPosition(position[0], position[1])
    return text_actor


def print_educational_summary(quality_metrics: dict, iterations: int, residual: float):
    """
    Print educational summary about the simulation.

    Args:
        quality_metrics: Mesh quality metrics dictionary
        iterations: Number of solver iterations
        residual: Final solver residual
    """
    print("\n" + "=" * 70)
    print("VTK Structured Mesh: CFD Heat Conduction Simulation")
    print("=" * 70)

    print("\n1. MESH STRUCTURE:")
    print(
        f"   - Dimensions: {quality_metrics['dimensions'][0]} x {quality_metrics['dimensions'][1]} nodes"
    )
    print(f"   - Total nodes: {quality_metrics['num_nodes']}")
    print(f"   - Total cells: {quality_metrics['num_cells']}")
    print(
        f"   - Cell spacing: dx={quality_metrics['dx']:.4f}m, dy={quality_metrics['dy']:.4f}m"
    )
    print(f"   - Aspect ratio: {quality_metrics['aspect_ratio']:.2f}")

    print("\n2. PHYSICAL PROBLEM:")
    print("   - Equation: 2D Laplace (steady-state heat conduction)")
    print("   - ∂²T/∂x² + ∂²T/∂y² = 0")
    print("   - Boundary conditions: Dirichlet (fixed temperatures)")
    print("   - Left wall: T = 100°C (heat source)")
    print("   - Other walls: T = 0°C (cold boundaries)")

    print("\n3. NUMERICAL METHOD:")
    print("   - Discretization: Finite Difference (5-point stencil)")
    print("   - Solver: Gauss-Seidel iterative method")
    print(f"   - Iterations to converge: {iterations}")
    print(f"   - Final residual: {residual:.2e}")

    print("\n4. KEY CFD CONCEPTS DEMONSTRATED:")
    print("   - Structured mesh topology (regular i-j indexing)")
    print("   - Implicit connectivity (no explicit cell-node mapping)")
    print("   - Vertex-centered data storage")
    print("   - Iterative solver convergence")
    print("   - Boundary condition enforcement")

    print("\n" + "=" * 70)


def main():
    """
    Main function demonstrating structured mesh for CFD heat conduction.

    This example:
    1. Creates a 2D structured mesh
    2. Solves the Laplace equation for steady-state heat conduction
    3. Visualizes the mesh and temperature field using VTK
    """
    # =========================================================================
    # SIMULATION PARAMETERS
    # =========================================================================
    # Domain dimensions (meters)
    lx, ly = 1.0, 1.0

    # Mesh resolution (number of nodes)
    nx, ny = 21, 21

    # Temperature boundary conditions (Celsius)
    T_hot = 100.0  # Left wall (heat source)
    T_cold = 0.0  # Other walls

    # Solver parameters
    tolerance = 1e-6
    max_iterations = 10000

    # =========================================================================
    # MESH GENERATION
    # =========================================================================
    print("Creating structured mesh...")
    grid = create_structured_mesh(nx, ny, lx, ly)
    quality_metrics = compute_mesh_quality(grid)

    # =========================================================================
    # SOLVER
    # =========================================================================
    print("Solving heat conduction equation...")
    T = initialize_temperature_field(nx, ny, T_init=T_cold)
    apply_boundary_conditions(T, T_hot, T_cold)
    T, iterations, residual = solve_laplace_gauss_seidel(
        T, T_hot, T_cold, tolerance, max_iterations
    )

    # Add temperature field to mesh
    add_temperature_to_mesh(grid, T)

    # Print educational summary
    print_educational_summary(quality_metrics, iterations, residual)

    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    print("\nStarting VTK visualization...")

    # Create visualization actors
    surface_actor, edge_actor, node_actor, lut = create_mesh_visualization_actors(grid)
    scalar_bar = create_scalar_bar(lut)

    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Add actors to renderer
    renderer.AddActor(surface_actor)
    renderer.AddActor(edge_actor)
    renderer.AddActor(node_actor)
    renderer.AddActor2D(scalar_bar)

    # Add title annotation
    title = create_text_annotation(
        "Structured Mesh: 2D Heat Conduction\nLeft wall: 100°C | Other walls: 0°C",
        (10, 550),
    )
    renderer.AddActor2D(title)

    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(900, 700)
    render_window.SetWindowName("VTK Structured Mesh: CFD Heat Conduction")

    # Create interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Set interaction style
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)

    # Add axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Setup camera for 2D view
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.SetParallelProjection(True)
    camera.SetPosition(0.5, 0.5, 5.0)
    camera.SetFocalPoint(0.5, 0.5, 0.0)
    camera.SetViewUp(0.0, 1.0, 0.0)
    renderer.ResetCamera()

    # Start visualization
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


if __name__ == "__main__":
    main()
