import vtk


# Enhanced grid creation with function
def create_structured_grid(dimensions):
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(dimensions)

    points = vtk.vtkPoints()
    for z in range(dimensions[2]):
        for y in range(dimensions[1]):
            for x in range(dimensions[0]):
                points.InsertNextPoint(x, y, z)
    grid.SetPoints(points)
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


# Grid dimensions
dimensions = [10, 10, 10]
structuredGrid = create_structured_grid(dimensions)

# Outline creation
outline = vtk.vtkStructuredGridOutlineFilter()
outline.SetInputData(structuredGrid)
outlineActor = create_actor(outline)

# Points visualization
pointGlyph = vtk.vtkVertexGlyphFilter()
pointGlyph.SetInputData(structuredGrid)
pointActor = create_actor(pointGlyph, point_size=1)

# Highlight node
highlightNode = vtk.vtkSphereSource()
highlightNode.SetCenter(5, 5, 5)
highlightNode.SetRadius(0.2)
highlightNodeActor = create_actor(highlightNode, color=(1, 0, 0))

# Edges visualization
edgeFilter = vtk.vtkExtractEdges()
edgeFilter.SetInputData(structuredGrid)
edgeActor = create_actor(edgeFilter, color=(0.5, 0.5, 0.5), line_width=1)

# Highlight edge
highlightEdge = vtk.vtkLineSource()
highlightEdge.SetPoint1(5, 5, 5)
highlightEdge.SetPoint2(6, 5, 5)
highlightEdgeActor = create_actor(highlightEdge, color=(0, 1, 0), line_width=3)

# Highlight cell
highlightCell = vtk.vtkCubeSource()
highlightCell.SetCenter(7.5, 7.5, 7.5)
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

# Adding actors and legend to renderer
actors = [
    outlineActor,
    pointActor,
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

# Create a 2D button for toggling points
pointsButton = Button2D(renderWindowInteractor, "Points", (10, 10), (100, 40))
pointsButton.add_click_callback(lambda: toggle_visibility(pointActor))

# Create a 2D button for toggling edges
edgesButton = Button2D(renderWindowInteractor, "Edges", (120, 10), (100, 40))
edgesButton.add_click_callback(lambda: toggle_visibility(edgeActor))

# Start the interaction
renderWindow.Render()
renderWindowInteractor.Start()
