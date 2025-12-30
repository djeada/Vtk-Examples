"""
VTK Cells: The Building Blocks of Computational Meshes

This module provides an educational demonstration of VTK cells, which are fundamental
to Computational Fluid Dynamics (CFD) and Finite Element Analysis (FEA).

What are Cells?
---------------
In VTK (and numerical simulations), a cell represents the topology or connectivity
between points. Cells define how points are connected to form geometric primitives.
Crucially, cells do NOT store positional data - they only reference point indices
stored separately in a vtkPoints object.

Cell Types by Dimension:
------------------------
- 0D Cells: Vertex (single point) - used for particle tracking, point clouds
- 1D Cells: Line, PolyLine - used for streamlines, boundary edges
- 2D Cells: Triangle, Quad, Polygon - used for surface meshes, 2D CFD domains
- 3D Cells: Tetrahedron, Hexahedron, Wedge, Pyramid - used for volumetric CFD/FEA meshes

CFD/FEA Context:
----------------
In CFD simulations, the computational domain is discretized into cells:
- Each cell represents a control volume where conservation equations are solved
- Cell quality (aspect ratio, skewness) affects numerical accuracy
- Different cell types offer trade-offs between accuracy and computational cost:
  * Tetrahedra: Easy to mesh complex geometries, but require more cells
  * Hexahedra: More accurate for aligned flows, harder to mesh
  * Triangles/Quads: Used in 2D simulations or surface meshes

This example demonstrates creating, connecting, and visualizing different VTK cell
types using vtkUnstructuredGrid, which can hold heterogeneous cell types - essential
for complex CFD meshes with different element types in different regions.
"""

import math

import vtk


# Distinct colors for each cell type (shared between legend and color lookup table)
CELL_TYPE_COLORS = [
    (1.0, 0.0, 0.0),    # Red - Vertex
    (1.0, 0.5, 0.0),    # Orange - Line
    (1.0, 1.0, 0.0),    # Yellow - Triangle
    (0.0, 1.0, 0.0),    # Green - Quad
    (0.0, 1.0, 1.0),    # Cyan - Polygon
    (0.0, 0.0, 1.0),    # Blue - Tetrahedron
    (1.0, 0.0, 1.0),    # Magenta - Hexahedron
]


def create_vertex_cell(points, base_id, x_offset):
    """
    Create a 0D vertex cell (single point).

    A vertex is the simplest cell type - it represents a single point in space.
    In CFD, vertices are used for:
    - Particle positions in Lagrangian particle tracking
    - Point data visualization (e.g., probes, sensors)
    - Point clouds from LiDAR or experimental data

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    points.InsertNextPoint(x_offset, 0.5, 0)

    vertex = vtk.vtkVertex()
    vertex.GetPointIds().SetId(0, base_id)

    return vertex, vtk.VTK_VERTEX, 1


def create_line_cell(points, base_id, x_offset):
    """
    Create a 1D line cell (connects two points).

    Lines are fundamental 1D elements used for:
    - Boundary edges of 2D meshes
    - Streamlines and pathlines in flow visualization
    - Structural elements like beams or trusses in FEA
    - Surface mesh edges for boundary condition application

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    points.InsertNextPoint(x_offset, 0, 0)
    points.InsertNextPoint(x_offset + 1.5, 1, 0)

    line = vtk.vtkLine()
    line.GetPointIds().SetId(0, base_id)
    line.GetPointIds().SetId(1, base_id + 1)

    return line, vtk.VTK_LINE, 2


def create_triangle_cell(points, base_id, x_offset):
    """
    Create a 2D triangle cell (connects three points).

    Triangles are the most versatile 2D cell type:
    - Always planar (3 points define a plane)
    - Can mesh any 2D geometry, no matter how complex
    - Used extensively in surface meshing
    - Building blocks for tetrahedral volume meshes
    - Linear triangles have constant gradients within each cell

    In CFD, triangular meshes are common for:
    - Complex geometries where structured meshes are impractical
    - Adaptive mesh refinement (easy to subdivide)
    - Unstructured mesh solvers

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    points.InsertNextPoint(x_offset, 0, 0)
    points.InsertNextPoint(x_offset + 1, 0, 0)
    points.InsertNextPoint(x_offset + 0.5, 0.866, 0)

    triangle = vtk.vtkTriangle()
    for i in range(3):
        triangle.GetPointIds().SetId(i, base_id + i)

    return triangle, vtk.VTK_TRIANGLE, 3


def create_quad_cell(points, base_id, x_offset):
    """
    Create a 2D quadrilateral cell (connects four points).

    Quadrilaterals (quads) are 4-sided 2D cells:
    - More accurate than triangles for aligned flows
    - Require structured or semi-structured mesh topology
    - Points must be ordered consistently (counter-clockwise or clockwise)
    - Can be non-planar (but planarity improves accuracy)

    In CFD, quad meshes are preferred when:
    - Flow is aligned with mesh (boundary layers, channels)
    - Higher accuracy per cell is needed
    - Memory efficiency is important (fewer cells than triangles)

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    # Points ordered counter-clockwise
    points.InsertNextPoint(x_offset, 0, 0)
    points.InsertNextPoint(x_offset + 1, 0, 0)
    points.InsertNextPoint(x_offset + 1, 1, 0)
    points.InsertNextPoint(x_offset, 1, 0)

    quad = vtk.vtkQuad()
    for i in range(4):
        quad.GetPointIds().SetId(i, base_id + i)

    return quad, vtk.VTK_QUAD, 4


def create_polygon_cell(points, base_id, x_offset, n_sides=5):
    """
    Create a 2D polygon cell (connects n points).

    Polygons are general n-sided 2D cells:
    - Flexible but more computationally expensive
    - Used in polyhedral mesh formats
    - Can represent complex face geometries
    - OpenFOAM and other solvers support polyhedral cells

    In CFD, polygons are used in:
    - Polyhedral mesh solvers (better gradient approximation)
    - Adaptive mesh refinement (hanging nodes)
    - Mesh agglomeration for multigrid methods

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning
        n_sides: Number of polygon sides (default: 5 for pentagon)

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    center_y = 0.5
    radius = 0.5

    for i in range(n_sides):
        # Rotate by -pi/2 so first vertex is at top (visually appealing orientation)
        angle = 2 * math.pi * i / n_sides - math.pi / 2
        x = x_offset + 0.5 + radius * math.cos(angle)
        y = center_y + radius * math.sin(angle)
        points.InsertNextPoint(x, y, 0)

    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(n_sides)
    for i in range(n_sides):
        polygon.GetPointIds().SetId(i, base_id + i)

    return polygon, vtk.VTK_POLYGON, n_sides


def create_tetra_cell(points, base_id, x_offset):
    """
    Create a 3D tetrahedron cell (connects four points).

    Tetrahedra are the simplest 3D volumetric cells:
    - 4 vertices, 6 edges, 4 triangular faces
    - Can mesh any 3D geometry (like triangles in 2D)
    - Most common unstructured 3D mesh element
    - Generated by automatic mesh generators (e.g., Delaunay)

    In CFD, tetrahedral meshes are used when:
    - Complex geometries preclude structured meshes
    - Rapid mesh generation is needed
    - Adaptive refinement is required

    Trade-offs:
    - Pros: Geometric flexibility, automatic generation
    - Cons: More cells needed, numerical diffusion in some schemes

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    # Base triangle
    points.InsertNextPoint(x_offset, 0, 0)
    points.InsertNextPoint(x_offset + 1, 0, 0)
    points.InsertNextPoint(x_offset + 0.5, 0.866, 0)
    # Apex
    points.InsertNextPoint(x_offset + 0.5, 0.289, 0.816)

    tetra = vtk.vtkTetra()
    for i in range(4):
        tetra.GetPointIds().SetId(i, base_id + i)

    return tetra, vtk.VTK_TETRA, 4


def create_hexahedron_cell(points, base_id, x_offset):
    """
    Create a 3D hexahedron cell (connects eight points).

    Hexahedra are 6-faced 3D cells (like a cube):
    - 8 vertices, 12 edges, 6 quadrilateral faces
    - Highest accuracy per cell for aligned flows
    - Preferred for boundary layer meshing
    - More difficult to generate for complex geometries

    In CFD, hexahedral meshes are preferred when:
    - Flow direction is known (aligned meshes)
    - Accuracy is critical (boundary layers, heat transfer)
    - Computational efficiency is important

    Trade-offs:
    - Pros: Fewer cells needed, better accuracy for aligned flows
    - Cons: Difficult automatic generation, limited to simpler geometries

    Args:
        points: vtkPoints object to add coordinates to
        base_id: Starting point index for this cell's points
        x_offset: X-coordinate offset for positioning

    Returns:
        tuple: (cell object, cell type constant, number of points added)
    """
    # Bottom face (z=0), counter-clockwise from above
    points.InsertNextPoint(x_offset, 0, 0)
    points.InsertNextPoint(x_offset + 1, 0, 0)
    points.InsertNextPoint(x_offset + 1, 1, 0)
    points.InsertNextPoint(x_offset, 1, 0)
    # Top face (z=1), counter-clockwise from above
    points.InsertNextPoint(x_offset, 0, 1)
    points.InsertNextPoint(x_offset + 1, 0, 1)
    points.InsertNextPoint(x_offset + 1, 1, 1)
    points.InsertNextPoint(x_offset, 1, 1)

    hexa = vtk.vtkHexahedron()
    for i in range(8):
        hexa.GetPointIds().SetId(i, base_id + i)

    return hexa, vtk.VTK_HEXAHEDRON, 8


def build_cell_visualization():
    """
    Build an unstructured grid containing all cell types for visualization.

    This function demonstrates how to:
    1. Create different cell types with proper point connectivity
    2. Store heterogeneous cells in vtkUnstructuredGrid
    3. Attach scalar data to color-code cells by type
    4. Set up proper VTK visualization pipeline

    vtkUnstructuredGrid is essential for CFD because:
    - It can hold any combination of cell types
    - It supports both 2D and 3D cells in the same mesh
    - It's the basis for most unstructured CFD mesh formats

    Returns:
        tuple: (vtkUnstructuredGrid, cell_type_names list)
    """
    points = vtk.vtkPoints()
    ugrid = vtk.vtkUnstructuredGrid()

    # Scalar array to store cell type for coloring
    cell_type_scalars = vtk.vtkFloatArray()
    cell_type_scalars.SetName("CellType")
    cell_type_scalars.SetNumberOfComponents(1)

    # Cell type names for legend (ordered by scalar value)
    cell_type_names = [
        "Vertex (0D)",
        "Line (1D)",
        "Triangle (2D)",
        "Quad (2D)",
        "Polygon (2D)",
        "Tetrahedron (3D)",
        "Hexahedron (3D)",
    ]

    # Track current point index and x position for layout
    current_point_id = 0
    x_position = 0.0
    spacing = 2.0

    # Create all cell types
    cell_creators = [
        (create_vertex_cell, 0),
        (create_line_cell, 1),
        (create_triangle_cell, 2),
        (create_quad_cell, 3),
        (create_polygon_cell, 4),
        (create_tetra_cell, 5),
        (create_hexahedron_cell, 6),
    ]

    for creator_func, scalar_value in cell_creators:
        cell, cell_type, num_points = creator_func(points, current_point_id, x_position)
        ugrid.InsertNextCell(cell_type, cell.GetPointIds())
        cell_type_scalars.InsertNextValue(scalar_value)
        current_point_id += num_points
        x_position += spacing

    ugrid.SetPoints(points)
    ugrid.GetCellData().SetScalars(cell_type_scalars)

    return ugrid, cell_type_names


def create_color_lookup_table(num_colors):
    """
    Create a color lookup table for distinguishing cell types.

    Uses distinct colors for each cell type to make the visualization
    clear and educational.

    Args:
        num_colors: Number of distinct colors needed

    Returns:
        vtkLookupTable configured with distinct colors
    """
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(num_colors)
    lut.SetRange(0, num_colors - 1)

    for i, (r, g, b) in enumerate(CELL_TYPE_COLORS[:num_colors]):
        lut.SetTableValue(i, r, g, b, 1.0)

    lut.Build()
    return lut


def visualize_cells(ugrid, cell_type_names):
    """
    Visualize the unstructured grid with color-coded cell types.

    Sets up a complete VTK rendering pipeline with:
    - Color-mapped cells based on cell type
    - Visible edges to show cell structure
    - Point markers to show vertices
    - Legend explaining cell types

    Args:
        ugrid: vtkUnstructuredGrid containing the cells
        cell_type_names: List of cell type names for legend
    """
    # Create lookup table for cell type coloring
    lut = create_color_lookup_table(len(cell_type_names))

    # Main surface mapper with cell coloring
    surface_mapper = vtk.vtkDataSetMapper()
    surface_mapper.SetInputData(ugrid)
    surface_mapper.SetScalarModeToUseCellData()
    surface_mapper.SelectColorArray("CellType")
    surface_mapper.SetLookupTable(lut)
    surface_mapper.SetScalarRange(0, len(cell_type_names) - 1)

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(surface_mapper)
    surface_actor.GetProperty().SetEdgeVisibility(1)
    surface_actor.GetProperty().SetEdgeColor(0.2, 0.2, 0.2)
    surface_actor.GetProperty().SetLineWidth(2)

    # Point visualization using glyphs
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetRadius(0.08)
    sphere_source.SetThetaResolution(16)
    sphere_source.SetPhiResolution(16)

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(ugrid)
    glyph.SetSourceConnection(sphere_source.GetOutputPort())
    glyph.SetScaleModeToDataScalingOff()

    point_mapper = vtk.vtkPolyDataMapper()
    point_mapper.SetInputConnection(glyph.GetOutputPort())
    point_mapper.ScalarVisibilityOff()

    point_actor = vtk.vtkActor()
    point_actor.SetMapper(point_mapper)
    point_actor.GetProperty().SetColor(0.1, 0.1, 0.1)

    # Create legend
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(len(cell_type_names))
    legend.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPositionCoordinate().SetValue(0.02, 0.15)
    legend.GetPosition2Coordinate().SetCoordinateSystemToNormalizedViewport()
    legend.GetPosition2Coordinate().SetValue(0.25, 0.5)
    legend.UseBackgroundOn()
    legend.SetBackgroundColor(0.1, 0.1, 0.1)

    # Create a single sphere source for legend entries and update it
    legend_sphere = vtk.vtkSphereSource()
    legend_sphere.Update()
    legend_output = legend_sphere.GetOutput()

    for i, (name, color) in enumerate(zip(cell_type_names, CELL_TYPE_COLORS)):
        legend.SetEntry(i, legend_output, name, color)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(point_actor)
    renderer.AddActor(legend)
    renderer.SetBackground(0.15, 0.15, 0.2)

    # Render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1200, 600)
    render_window.SetWindowName("VTK Cells: Building Blocks of CFD Meshes")

    # Interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Camera positioning for good initial view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(6, -8, 8)
    camera.SetFocalPoint(6, 0.5, 0.3)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()

    # Add axes widget for orientation
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(interactor)
    axes_widget.SetViewport(0.0, 0.0, 0.15, 0.25)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Start interaction
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


def print_cell_info(ugrid, cell_type_names):
    """
    Print educational information about the cells in the grid.

    Demonstrates how to query cell properties in VTK, which is
    useful for mesh quality analysis in CFD preprocessing.

    Args:
        ugrid: vtkUnstructuredGrid containing the cells
        cell_type_names: List of cell type names
    """
    print("\n" + "=" * 70)
    print("VTK Cells: Educational Overview for CFD/FEA")
    print("=" * 70)

    print(f"\nTotal points in grid: {ugrid.GetNumberOfPoints()}")
    print(f"Total cells in grid:  {ugrid.GetNumberOfCells()}")

    print("\n" + "-" * 70)
    print("Cell Details:")
    print("-" * 70)

    for i in range(ugrid.GetNumberOfCells()):
        cell = ugrid.GetCell(i)
        cell_type = ugrid.GetCellType(i)
        num_points = cell.GetNumberOfPoints()
        num_edges = cell.GetNumberOfEdges()
        num_faces = cell.GetNumberOfFaces()

        print(f"\n{cell_type_names[i]}:")
        print(f"  VTK Cell Type ID: {cell_type}")
        print(f"  Number of vertices: {num_points}")
        print(f"  Number of edges: {num_edges}")
        print(f"  Number of faces: {num_faces}")

        # Print vertex coordinates
        point_ids = [cell.GetPointId(j) for j in range(num_points)]
        print(f"  Point IDs: {point_ids}")

    print("\n" + "=" * 70)
    print("Key Concepts for CFD:")
    print("=" * 70)
    print("""
1. CELL TOPOLOGY: Cells define connectivity between points, not positions.
   The same topology can represent different physical shapes.

2. CELL QUALITY: In CFD, cell quality affects numerical accuracy:
   - Aspect ratio: ratio of longest to shortest edge
   - Skewness: deviation from ideal element shape
   - Orthogonality: angle between face normal and cell center connections

3. CELL SELECTION: Choose cell types based on:
   - Geometry complexity (simple → hexahedra, complex → tetrahedra)
   - Flow alignment (aligned → quads/hexahedra, arbitrary → triangles/tetrahedra)
   - Accuracy requirements (higher → hexahedra, standard → tetrahedra)
   - Meshing time constraints (fast → tetrahedra, slower → hexahedra)

4. HYBRID MESHES: Production CFD meshes often combine cell types:
   - Hexahedra/prisms in boundary layers for accuracy
   - Tetrahedra in the bulk for geometric flexibility
   - Pyramids to transition between quad and triangle faces
""")


def main():
    """
    Main function demonstrating VTK cells for CFD/FEA applications.

    Creates and visualizes different VTK cell types to illustrate:
    - How cells define mesh topology
    - The relationship between points and cells
    - Different cell types used in CFD meshes
    - How to store heterogeneous cells in vtkUnstructuredGrid
    """
    # Build the unstructured grid with all cell types
    ugrid, cell_type_names = build_cell_visualization()

    # Print educational information to console
    print_cell_info(ugrid, cell_type_names)

    # Visualize the cells
    visualize_cells(ugrid, cell_type_names)


if __name__ == "__main__":
    main()
