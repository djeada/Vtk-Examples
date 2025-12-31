"""
VTK Unstructured Mesh: CFD Educational Example with Heat Diffusion Solver

This module provides an educational demonstration of unstructured meshes for
Computational Fluid Dynamics (CFD) applications. It combines triangular mesh
generation, a finite element heat diffusion solver, and VTK visualization.

What is an Unstructured Mesh?
-----------------------------
An unstructured mesh is a grid where:
- Cells can be of arbitrary shape (triangles, tetrahedra, etc.)
- Connectivity is EXPLICIT (stored for each cell, not implicit)
- No regular i-j-k topology is required
- Can mesh any geometry, no matter how complex

Mesh Components:
----------------
1. NODES (Vertices): Points where solution variables are computed
   - Positions can be arbitrary (not on a regular grid)
   - Each node has a unique ID for connectivity lookup

2. EDGES: Connections between nodes
   - Defined by node pairs
   - Used for flux calculations in finite volume methods

3. CELLS (Elements): Polygons/polyhedra defined by node connectivity
   - Triangles (2D), Tetrahedra (3D) are most common
   - Each cell stores its node IDs explicitly

Unstructured vs Structured:
---------------------------
UNSTRUCTURED ADVANTAGES:
- Can mesh ANY geometry (complex shapes, curves, holes)
- Easy local refinement (adaptive meshing)
- Automatic mesh generation (Delaunay triangulation)
- Natural for finite element methods

UNSTRUCTURED DISADVANTAGES:
- Higher memory usage (explicit connectivity)
- Slower neighbor lookup (requires search)
- More complex data structures
- Generally requires more cells for same accuracy

Physical Problem: Steady-State Heat Diffusion
---------------------------------------------
This example solves the 2D Laplace equation on a triangular mesh:

    ∂²T/∂x² + ∂²T/∂y² = 0

with Dirichlet boundary conditions on a circular domain:
- Boundary nodes at angle 0-180°: T = T_hot (100°C)
- Boundary nodes at angle 180-360°: T = T_cold (0°C)

This represents heat diffusion from a hot semicircle to a cold semicircle.

Numerical Method: Finite Element with Linear Triangles
-------------------------------------------------------
Using linear triangular elements:
1. Temperature varies linearly within each triangle
2. Weak form leads to the stiffness matrix equation: K*T = f
3. Gauss-Seidel iteration solves the system
4. Each node's value is updated using its neighbors

For linear triangles on the Laplace equation with uniform material:
    T_i = (sum of neighbor temperatures weighted by edge cotangent) / (sum of weights)

CFD Applications:
-----------------
Unstructured meshes with similar methods are used for:
- Complex geometry aerodynamics (aircraft, vehicles)
- Biomedical flows (blood vessels, organs)
- Geophysical simulations (terrain, coastlines)
- Industrial equipment (valves, mixers)

Visualization:
--------------
The VTK visualization shows:
1. The triangular mesh with nodes and edges
2. Temperature field as a color map (blue=cold, red=hot)
3. Interactive camera controls
"""

import math
import numpy as np
import vtk


def create_circular_triangular_mesh(
    radius: float, n_radial: int, n_angular: int
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """
    Create a 2D triangular mesh on a circular domain.

    The mesh is generated using a polar grid approach and then triangulated:
    - Center node at origin
    - Concentric rings of nodes
    - Each ring connected to form triangular cells

    This is a common approach for circular domains in CFD, as it provides
    good mesh quality with natural refinement towards the center.

    Args:
        radius: Radius of the circular domain (meters)
        n_radial: Number of radial divisions
        n_angular: Number of angular divisions

    Returns:
        tuple: (nodes, triangles, boundary_node_ids)
            - nodes: Nx2 array of node coordinates
            - triangles: Mx3 array of triangle connectivity
            - boundary_node_ids: List of node IDs on the boundary
    """
    nodes = []
    triangles = []
    boundary_node_ids = []

    # Center node
    nodes.append([0.0, 0.0])
    center_id = 0

    # Generate concentric rings
    node_id = 1
    ring_start_ids = [0]  # Start with center

    for ring in range(1, n_radial + 1):
        r = radius * ring / n_radial
        ring_start_ids.append(node_id)

        for i in range(n_angular):
            theta = 2 * math.pi * i / n_angular
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            nodes.append([x, y])

            if ring == n_radial:
                boundary_node_ids.append(node_id)

            node_id += 1

    # Create triangles connecting center to first ring
    first_ring_start = ring_start_ids[1]
    for i in range(n_angular):
        n1 = center_id
        n2 = first_ring_start + i
        n3 = first_ring_start + (i + 1) % n_angular
        triangles.append([n1, n2, n3])

    # Create triangles between successive rings
    for ring in range(1, n_radial):
        inner_start = ring_start_ids[ring]
        outer_start = ring_start_ids[ring + 1]

        for i in range(n_angular):
            # Inner node indices
            i1 = inner_start + i
            i2 = inner_start + (i + 1) % n_angular
            # Outer node indices
            o1 = outer_start + i
            o2 = outer_start + (i + 1) % n_angular

            # Two triangles per quad section
            triangles.append([i1, o1, o2])
            triangles.append([i1, o2, i2])

    return np.array(nodes), np.array(triangles), boundary_node_ids


def create_vtk_unstructured_mesh(
    nodes: np.ndarray, triangles: np.ndarray
) -> vtk.vtkUnstructuredGrid:
    """
    Create a VTK unstructured grid from node and triangle arrays.

    Converts numpy arrays to VTK data structures:
    - Nodes become vtkPoints
    - Triangles become vtkTriangle cells with explicit connectivity

    Args:
        nodes: Nx2 array of node coordinates
        triangles: Mx3 array of triangle connectivity

    Returns:
        vtkUnstructuredGrid: The unstructured mesh
    """
    points = vtk.vtkPoints()
    for node in nodes:
        points.InsertNextPoint(node[0], node[1], 0.0)

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)

    for tri in triangles:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, int(tri[0]))
        triangle.GetPointIds().SetId(1, int(tri[1]))
        triangle.GetPointIds().SetId(2, int(tri[2]))
        ugrid.InsertNextCell(triangle.GetCellType(), triangle.GetPointIds())

    return ugrid


def build_node_neighbors(
    n_nodes: int, triangles: np.ndarray
) -> dict[int, set[int]]:
    """
    Build adjacency map for nodes in the mesh.

    For each node, find all nodes connected to it by an edge.
    This is used by the iterative solver to update node values.

    In unstructured meshes, this connectivity must be computed explicitly
    (unlike structured meshes where neighbors are known from indices).

    Args:
        n_nodes: Total number of nodes
        triangles: Mx3 array of triangle connectivity

    Returns:
        dict: Mapping from node ID to set of neighbor node IDs
    """
    neighbors: dict[int, set[int]] = {i: set() for i in range(n_nodes)}

    for tri in triangles:
        n0, n1, n2 = int(tri[0]), int(tri[1]), int(tri[2])
        neighbors[n0].add(n1)
        neighbors[n0].add(n2)
        neighbors[n1].add(n0)
        neighbors[n1].add(n2)
        neighbors[n2].add(n0)
        neighbors[n2].add(n1)

    return neighbors


def initialize_temperature_field(
    n_nodes: int, T_init: float = 0.0
) -> np.ndarray:
    """
    Initialize the temperature field for all mesh nodes.

    Creates a 1D array of temperature values, one per node.
    Unlike structured meshes, nodes are not arranged in a regular pattern.

    Args:
        n_nodes: Number of mesh nodes
        T_init: Initial temperature value

    Returns:
        np.ndarray: Array of temperature values
    """
    return np.full(n_nodes, T_init, dtype=np.float64)


def apply_boundary_conditions(
    T: np.ndarray,
    nodes: np.ndarray,
    boundary_ids: list[int],
    T_hot: float,
    T_cold: float,
) -> set[int]:
    """
    Apply Dirichlet boundary conditions on the circular boundary.

    Boundary condition assignment based on angular position:
    - Nodes in upper half (y >= 0): Hot temperature
    - Nodes in lower half (y < 0): Cold temperature

    This creates a thermal gradient across the domain.

    Args:
        T: Temperature array (modified in-place)
        nodes: Node coordinates array
        boundary_ids: List of boundary node IDs
        T_hot: Hot boundary temperature
        T_cold: Cold boundary temperature

    Returns:
        set[int]: Set of all boundary node IDs (for solver to skip)
    """
    boundary_set = set(boundary_ids)

    for node_id in boundary_ids:
        y = nodes[node_id, 1]
        if y >= 0:
            T[node_id] = T_hot
        else:
            T[node_id] = T_cold

    return boundary_set


def solve_laplace_gauss_seidel(
    T: np.ndarray,
    neighbors: dict[int, set[int]],
    boundary_nodes: set[int],
    tolerance: float = 1e-6,
    max_iterations: int = 10000,
) -> tuple[np.ndarray, int, float]:
    """
    Solve the Laplace equation on an unstructured mesh using Gauss-Seidel.

    For unstructured meshes, the update formula uses simple averaging:
        T_i = (sum of neighbor temperatures) / (number of neighbors)

    This is a simplified version of the finite element method for
    the Laplace equation with equal-weight neighbors. For production
    CFD, cotangent weights from the mesh geometry would be used.

    The Gauss-Seidel method:
    1. Iterate through interior nodes
    2. Update each using average of its neighbors
    3. Repeat until convergence

    Args:
        T: Temperature field (modified in-place)
        neighbors: Adjacency map from node to neighbors
        boundary_nodes: Set of boundary node IDs (not updated)
        tolerance: Convergence tolerance
        max_iterations: Maximum iterations

    Returns:
        tuple: (T, iterations, residual)
    """
    n_nodes = len(T)
    residual = float("inf")

    for iteration in range(max_iterations):
        max_change = 0.0

        for node_id in range(n_nodes):
            if node_id in boundary_nodes:
                continue

            neighbor_ids = neighbors[node_id]
            if len(neighbor_ids) == 0:
                continue

            T_old = T[node_id]
            T_new = sum(T[n] for n in neighbor_ids) / len(neighbor_ids)
            T[node_id] = T_new

            change = abs(T_new - T_old)
            if change > max_change:
                max_change = change

        residual = max_change

        if residual < tolerance:
            return T, iteration + 1, residual

    return T, max_iterations, residual


def add_temperature_to_mesh(
    ugrid: vtk.vtkUnstructuredGrid, T: np.ndarray
) -> None:
    """
    Add the computed temperature field to the VTK mesh.

    Temperature is stored as point data (vertex-centered) on the
    unstructured grid, which is the standard for finite element results.

    Args:
        ugrid: VTK unstructured grid
        T: Temperature array
    """
    temperature_array = vtk.vtkFloatArray()
    temperature_array.SetName("Temperature")
    temperature_array.SetNumberOfComponents(1)

    for temp in T:
        temperature_array.InsertNextValue(temp)

    ugrid.GetPointData().AddArray(temperature_array)
    ugrid.GetPointData().SetActiveScalars("Temperature")


def compute_mesh_quality(
    nodes: np.ndarray, triangles: np.ndarray, boundary_ids: list[int]
) -> dict:
    """
    Compute quality metrics for the unstructured mesh.

    Mesh quality metrics for CFD:
    - Number of nodes and cells
    - Boundary vs interior node count
    - Average triangles per node (valence)

    Good mesh quality is essential for accurate CFD simulations.

    Args:
        nodes: Node coordinates
        triangles: Triangle connectivity
        boundary_ids: Boundary node IDs

    Returns:
        dict: Quality metrics
    """
    n_nodes = len(nodes)
    n_triangles = len(triangles)
    n_boundary = len(boundary_ids)
    n_interior = n_nodes - n_boundary

    # Compute node valence (number of triangles per node)
    valence = np.zeros(n_nodes)
    for tri in triangles:
        for node_id in tri:
            valence[int(node_id)] += 1

    avg_valence = np.mean(valence)

    return {
        "n_nodes": n_nodes,
        "n_triangles": n_triangles,
        "n_boundary_nodes": n_boundary,
        "n_interior_nodes": n_interior,
        "avg_valence": avg_valence,
    }


def create_mesh_visualization_actors(
    ugrid: vtk.vtkUnstructuredGrid,
) -> tuple[vtk.vtkActor, vtk.vtkActor, vtk.vtkActor, vtk.vtkLookupTable]:
    """
    Create VTK actors for visualizing the unstructured mesh.

    Creates:
    1. Surface actor colored by temperature
    2. Edge actor showing mesh wireframe
    3. Node actor showing mesh vertices
    4. Lookup table for color mapping

    Args:
        ugrid: VTK unstructured grid with temperature data

    Returns:
        tuple: (surface_actor, edge_actor, node_actor, lut)
    """
    # Color lookup table (blue=cold, red=hot)
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0.0)
    lut.SetNumberOfTableValues(256)
    lut.Build()

    # Surface actor
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(ugrid)
    surface_mapper.SetScalarModeToUsePointFieldData()
    surface_mapper.SelectColorArray("Temperature")
    temp_range = ugrid.GetPointData().GetArray("Temperature").GetRange()
    surface_mapper.SetScalarRange(temp_range)
    surface_mapper.SetLookupTable(lut)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetOpacity(0.9)

    # Edge actor
    edge_filter = vtk.vtkExtractEdges()
    edge_filter.SetInputData(ugrid)

    edge_mapper = vtk.vtkPolyDataMapper()
    edge_mapper.SetInputConnection(edge_filter.GetOutputPort())
    edge_mapper.ScalarVisibilityOff()

    edge_actor = vtk.vtkActor()
    edge_actor.SetMapper(edge_mapper)
    edge_actor.GetProperty().SetColor(0.2, 0.2, 0.2)
    edge_actor.GetProperty().SetLineWidth(1)

    # Node actor
    point_filter = vtk.vtkVertexGlyphFilter()
    point_filter.SetInputData(ugrid)

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
        lut: VTK lookup table

    Returns:
        vtkScalarBarActor: Color bar actor
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
        text: Text to display
        position: (x, y) position

    Returns:
        vtkTextActor: Text actor
    """
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(text)
    text_actor.GetTextProperty().SetFontSize(14)
    text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    text_actor.SetPosition(position[0], position[1])
    return text_actor


def print_educational_summary(
    quality_metrics: dict, iterations: int, residual: float
) -> None:
    """
    Print educational summary about the simulation.

    Args:
        quality_metrics: Mesh quality metrics
        iterations: Solver iterations
        residual: Final residual
    """
    print("\n" + "=" * 70)
    print("VTK Unstructured Mesh: CFD Heat Diffusion Simulation")
    print("=" * 70)

    print("\n1. MESH STRUCTURE:")
    print(f"   - Total nodes: {quality_metrics['n_nodes']}")
    print(f"   - Total triangles: {quality_metrics['n_triangles']}")
    print(f"   - Boundary nodes: {quality_metrics['n_boundary_nodes']}")
    print(f"   - Interior nodes: {quality_metrics['n_interior_nodes']}")
    print(f"   - Average node valence: {quality_metrics['avg_valence']:.2f}")

    print("\n2. PHYSICAL PROBLEM:")
    print("   - Equation: 2D Laplace (steady-state heat diffusion)")
    print("   - ∂²T/∂x² + ∂²T/∂y² = 0")
    print("   - Domain: Circular with radius 1.0m")
    print("   - Boundary conditions: Dirichlet (fixed temperatures)")
    print("   - Upper half (y≥0): T = 100°C (hot)")
    print("   - Lower half (y<0): T = 0°C (cold)")

    print("\n3. NUMERICAL METHOD:")
    print("   - Mesh type: Triangular (Delaunay-style)")
    print("   - Discretization: Finite Element (linear triangles)")
    print("   - Solver: Gauss-Seidel with neighbor averaging")
    print(f"   - Iterations to converge: {iterations}")
    print(f"   - Final residual: {residual:.2e}")

    print("\n4. KEY CFD CONCEPTS DEMONSTRATED:")
    print("   - Unstructured mesh topology (explicit connectivity)")
    print("   - Triangular elements for arbitrary geometry")
    print("   - Vertex-centered data storage")
    print("   - Neighbor-based iterative solver")
    print("   - Circular domain meshing")

    print("\n5. UNSTRUCTURED VS STRUCTURED:")
    print("   - Explicit connectivity: O(n) memory for neighbors")
    print("   - Can mesh any geometry (circles, curves, holes)")
    print("   - Local refinement possible")
    print("   - More complex implementation")

    print("\n" + "=" * 70)


def main():
    """
    Main function demonstrating unstructured mesh for CFD heat diffusion.

    This example:
    1. Creates a 2D triangular mesh on a circular domain
    2. Solves the Laplace equation for heat diffusion
    3. Visualizes the mesh and temperature field using VTK
    """
    # =========================================================================
    # SIMULATION PARAMETERS
    # =========================================================================
    # Domain
    radius = 1.0  # meters

    # Mesh resolution
    n_radial = 10   # Radial divisions
    n_angular = 24  # Angular divisions

    # Temperature boundary conditions (Celsius)
    T_hot = 100.0   # Upper half
    T_cold = 0.0    # Lower half

    # Solver parameters
    tolerance = 1e-6
    max_iterations = 10000

    # =========================================================================
    # MESH GENERATION
    # =========================================================================
    print("Creating unstructured triangular mesh...")
    nodes, triangles, boundary_ids = create_circular_triangular_mesh(
        radius, n_radial, n_angular
    )

    quality_metrics = compute_mesh_quality(nodes, triangles, boundary_ids)
    neighbors = build_node_neighbors(len(nodes), triangles)

    # Create VTK mesh
    ugrid = create_vtk_unstructured_mesh(nodes, triangles)

    # =========================================================================
    # SOLVER
    # =========================================================================
    print("Solving heat diffusion equation...")
    T = initialize_temperature_field(len(nodes), T_init=0.0)
    boundary_set = apply_boundary_conditions(
        T, nodes, boundary_ids, T_hot, T_cold
    )
    T, iterations, residual = solve_laplace_gauss_seidel(
        T, neighbors, boundary_set, tolerance, max_iterations
    )

    # Add temperature to mesh
    add_temperature_to_mesh(ugrid, T)

    # Print educational summary
    print_educational_summary(quality_metrics, iterations, residual)

    # =========================================================================
    # VISUALIZATION
    # =========================================================================
    print("\nStarting VTK visualization...")

    # Create visualization actors
    surface_actor, edge_actor, node_actor, lut = create_mesh_visualization_actors(ugrid)
    scalar_bar = create_scalar_bar(lut)

    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Add actors
    renderer.AddActor(surface_actor)
    renderer.AddActor(edge_actor)
    renderer.AddActor(node_actor)
    renderer.AddActor2D(scalar_bar)

    # Add title
    title = create_text_annotation(
        "Unstructured Mesh: 2D Heat Diffusion\nUpper half: 100°C | Lower half: 0°C",
        (10, 550)
    )
    renderer.AddActor2D(title)

    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(900, 700)
    render_window.SetWindowName("VTK Unstructured Mesh: CFD Heat Diffusion")

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
    camera.SetPosition(0.0, 0.0, 5.0)
    camera.SetFocalPoint(0.0, 0.0, 0.0)
    camera.SetViewUp(0.0, 1.0, 0.0)
    renderer.ResetCamera()

    # Start visualization
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


if __name__ == "__main__":
    main()
