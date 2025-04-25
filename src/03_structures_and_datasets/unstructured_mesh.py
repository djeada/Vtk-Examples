import random

import vtk


# Function to create an unstructured grid with a variety of cell types
def create_unstructured_grid():
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()

    # Define points for a larger structure
    point_ids = [
        points.InsertNextPoint(x, y, z)
        for x in range(4)
        for y in range(4)
        for z in range(4)
    ]

    for x in range(4):
        for y in range(4):
            for z in range(4):
                # Add some randomness to the node positions
                dx, dy, dz = (
                    random.uniform(-0.2, 0.2),
                    random.uniform(-0.2, 0.2),
                    random.uniform(-0.2, 0.2),
                )
                points.InsertNextPoint(x + dx, y + dy, z + dz)

    grid.SetPoints(points)

    # Add different types of cells to the unstructured grid
    for i in range(0, len(point_ids), 4):
        # Adding a hexahedron cell
        hexahedron = vtk.vtkHexahedron()
        for j in range(8):
            hexahedron.GetPointIds().SetId(j, (i + j) % len(point_ids))
        grid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())

        # Adding a tetrahedron cell
        tetra = vtk.vtkTetra()
        for j in range(4):
            tetra.GetPointIds().SetId(j, (i + j) % len(point_ids))
        grid.InsertNextCell(tetra.GetCellType(), tetra.GetPointIds())

    return grid


# Function to create a standard mapper and actor
def create_actor(source, color=None, point_size=None, line_width=None):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(source.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    if color:
        actor.GetProperty().SetColor(color)
    if point_size:
        actor.GetProperty().SetPointSize(point_size)
    if line_width:
        actor.GetProperty().SetLineWidth(line_width)

    return actor


# Create an unstructured grid
unstructuredGrid = create_unstructured_grid()

# Nodes visualization with smaller size
nodeGlyph = vtk.vtkGlyph3D()
nodeGlyph.SetInputData(unstructuredGrid)
sphereSource = vtk.vtkSphereSource()
sphereSource.SetRadius(0.05)  # Smaller radius for nodes
nodeGlyph.SetSourceConnection(sphereSource.GetOutputPort())
nodeActor = create_actor(nodeGlyph, point_size=1)

# Edges visualization with smaller width
edgeFilter = vtk.vtkExtractEdges()
edgeFilter.SetInputData(unstructuredGrid)
edgeActor = create_actor(edgeFilter, line_width=0.5)

# Highlight node
highlightNode = vtk.vtkSphereSource()
highlightNode.SetCenter(2, 2, 2)
highlightNode.SetRadius(0.1)
highlightNodeActor = create_actor(highlightNode, color=(1, 0, 0))

# Highlight edge
highlightEdge = vtk.vtkLineSource()
highlightEdge.SetPoint1(1, 1, 1)
highlightEdge.SetPoint2(2, 1, 1)
highlightEdgeActor = create_actor(highlightEdge, color=(0, 1, 0), line_width=2)

# Highlight cell
highlightCell = vtk.vtkCubeSource()
highlightCell.SetCenter(1.5, 1.5, 1.5)
highlightCell.SetXLength(1)
highlightCell.SetYLength(1)
highlightCell.SetZLength(1)
highlightCellActor = create_actor(highlightCell, color=(0, 0, 1))

# Legend creation
legend = vtk.vtkLegendBoxActor()
legend.SetNumberOfEntries(3)
legend.SetEntry(0, vtk.vtkPolyData(), "Node (Red)", [1, 0, 0])
legend.SetEntry(1, vtk.vtkPolyData(), "Edge (Green)", [0, 1, 0])
legend.SetEntry(2, vtk.vtkPolyData(), "Cell (Blue)", [0, 0, 1])

# Renderer and Render Window setup
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.SetSize(800, 600)
renderer.SetBackground(0.2, 0.3, 0.4)

# Render Window Interactor
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Adding actors to renderer
actors = [
    nodeActor,
    edgeActor,
    highlightNodeActor,
    highlightEdgeActor,
    highlightCellActor,
    legend,
]
for actor in actors:
    renderer.AddActor(actor)


# Visibility toggle functions
def toggle_visibility(actor):
    actor.SetVisibility(not actor.GetVisibility())
    renderWindow.Render()


# Button creation and setup
from src.common.button import Button2D

# Create a 2D button for toggling nodes
nodesButton = Button2D(renderWindowInteractor, "Nodes", (10, 10), (100, 40))
nodesButton.add_click_callback(lambda: toggle_visibility(nodeActor))

# Create a 2D button for toggling edges
edgesButton = Button2D(renderWindowInteractor, "Edges", (120, 10), (100, 40))
edgesButton.add_click_callback(lambda: toggle_visibility(edgeActor))

# Start the interaction
renderWindow.Render()
renderWindowInteractor.Start()
