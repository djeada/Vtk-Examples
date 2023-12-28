import vtk

from src.common.button import Button2D

# Create a structured grid
dimensions = [10, 10, 10]
structuredGrid = vtk.vtkStructuredGrid()
structuredGrid.SetDimensions(dimensions)

# Create points for the grid
points = vtk.vtkPoints()
for z in range(dimensions[2]):
    for y in range(dimensions[1]):
        for x in range(dimensions[0]):
            points.InsertNextPoint(x, y, z)
structuredGrid.SetPoints(points)

# Create a mapper and actor for the grid outline
outline = vtk.vtkStructuredGridOutlineFilter()
outline.SetInputData(structuredGrid)

outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInputConnection(outline.GetOutputPort())

outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)

# Visualize all nodes with a smaller point size
pointGlyph = vtk.vtkVertexGlyphFilter()
pointGlyph.SetInputData(structuredGrid)

pointMapper = vtk.vtkPolyDataMapper()
pointMapper.SetInputConnection(pointGlyph.GetOutputPort())

pointActor = vtk.vtkActor()
pointActor.SetMapper(pointMapper)
pointActor.GetProperty().SetPointSize(1)  # Reduced point size

# Highlight a specific node with an increased radius and red color
highlightNode = vtk.vtkSphereSource()
highlightNode.SetCenter(5, 5, 5)
highlightNode.SetRadius(0.2)

highlightNodeMapper = vtk.vtkPolyDataMapper()
highlightNodeMapper.SetInputConnection(highlightNode.GetOutputPort())

highlightNodeActor = vtk.vtkActor()
highlightNodeActor.SetMapper(highlightNodeMapper)
highlightNodeActor.GetProperty().SetColor(1, 0, 0)  # Red color

# Visualize all edges with a subtle representation
edgeFilter = vtk.vtkExtractEdges()
edgeFilter.SetInputData(structuredGrid)

edgeMapper = vtk.vtkPolyDataMapper()
edgeMapper.SetInputConnection(edgeFilter.GetOutputPort())

edgeActor = vtk.vtkActor()
edgeActor.SetMapper(edgeMapper)
edgeActor.GetProperty().SetColor(0.5, 0.5, 0.5)  # Subtle gray color
edgeActor.GetProperty().SetLineWidth(1)  # Thin line width

# Highlight an edge with increased line width and green color
highlightEdge = vtk.vtkLineSource()
highlightEdge.SetPoint1(5, 5, 5)
highlightEdge.SetPoint2(6, 5, 5)

highlightEdgeMapper = vtk.vtkPolyDataMapper()
highlightEdgeMapper.SetInputConnection(highlightEdge.GetOutputPort())

highlightEdgeActor = vtk.vtkActor()
highlightEdgeActor.SetMapper(highlightEdgeMapper)
highlightEdgeActor.GetProperty().SetColor(0, 1, 0)  # Green color
highlightEdgeActor.GetProperty().SetLineWidth(3)  # Increased line width


# Highlight a cell with blue color
highlightCell = vtk.vtkCubeSource()
highlightCell.SetCenter(7.5, 7.5, 7.5)
highlightCell.SetXLength(1)
highlightCell.SetYLength(1)
highlightCell.SetZLength(1)

highlightCellMapper = vtk.vtkPolyDataMapper()
highlightCellMapper.SetInputConnection(highlightCell.GetOutputPort())

highlightCellActor = vtk.vtkActor()
highlightCellActor.SetMapper(highlightCellMapper)
highlightCellActor.GetProperty().SetColor(0, 0, 1)  # Blue color

# Create a legend with detailed descriptions
legend = vtk.vtkLegendBoxActor()
legend.SetNumberOfEntries(3)
legend.SetEntry(0, vtk.vtkPolyData(), "Node (Red)", [1, 0, 0])
legend.SetEntry(1, vtk.vtkPolyData(), "Edge (Green)", [0, 1, 0])
legend.SetEntry(2, vtk.vtkPolyData(), "Cell (Blue)", [0, 0, 1])

# Renderer and Render Window setup
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# Render Window Interactor
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Add actors and the legend to the renderer
renderer.AddActor(outlineActor)
renderer.AddActor(pointActor)
renderer.AddActor(edgeActor)  # Adding the actor for non-highlighted edges
renderer.AddActor(highlightNodeActor)
renderer.AddActor(highlightEdgeActor)
renderer.AddActor(highlightCellActor)
renderer.AddActor(legend)

# Set a more visually appealing background color
renderer.SetBackground(0.2, 0.3, 0.4)
renderWindow.SetSize(800, 600)


def toggle_points_visibility():
    isVisible = pointActor.GetVisibility()
    pointActor.SetVisibility(not isVisible)
    renderWindow.Render()


def toggle_edges_visibility():
    isVisible = edgeActor.GetVisibility()
    edgeActor.SetVisibility(not isVisible)
    highlightEdgeActor.SetVisibility(not isVisible)
    renderWindow.Render()


# Create a 2D button for toggling points
pointsButton = Button2D(renderWindowInteractor, "Points", (10, 10), (100, 40))
pointsButton.add_click_callback(toggle_points_visibility)

# Create a 2D button for toggling edges
edgesButton = Button2D(renderWindowInteractor, "Edges", (120, 10), (100, 40))
edgesButton.add_click_callback(toggle_edges_visibility)

# Start the interaction
renderWindow.Render()
renderWindowInteractor.Start()
