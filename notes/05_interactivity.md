## Interactivity and User Interface

### Overview
* VTK provides various tools for creating interactive visualizations and user interfaces
* Some popular techniques include:
  - Picking and Selection
  - Cutting Planes
  - 3D Widgets
  - Event Handling and Callbacks

### Picking and Selection
* Picking: interactively selecting an object or point within the visualization
* Selection: highlighting or extracting a subset of data based on specific criteria
* Example classes:
  - vtkCellPicker: selects cells or points within a 3D scene
  - vtkExtractSelection: extracts a subset of data based on the selection criteria

### Cutting Planes
* Cutting planes: interactively slicing through a 3D object to reveal its interior
* Useful for inspecting the internal structure of complex objects
* Example classes:
  - vtkPlane: defines a cutting plane
  - vtkCutter: slices through an object using a vtkPlane

### 3D Widgets
* 3D widgets: interactive tools for manipulating data or objects within a 3D scene
* Can be used for tasks such as positioning lights or cameras, or adjusting clipping planes
* Example classes:
  - vtkBoxWidget: interactively resizes a bounding box
  - vtkPlaneWidget: manipulates a plane within a 3D scene
  - vtkSliderWidget: adjusts scalar values using a slider

### Event Handling and Callbacks
* Event handling: responding to user input or other events within a VTK application
* Callbacks: user-defined functions that are executed in response to specific events
* Example classes:
  - vtkRenderWindowInteractor: manages user input and event handling
  - vtkCallbackCommand: associates a user-defined callback function with a specific event

## Example: Cutting Plane Interaction
```python
import vtk

# Create a sphere source
sphere = vtk.vtkSphereSource()
sphere.SetRadius(1)

# Create a mapper and actor for the sphere
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphere.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Create a cutting plane and cutter
plane = vtk.vtkPlane()
cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(sphere.GetOutputPort())

# Create a mapper and actor for the cut geometry
cut_mapper = vtk.vtkPolyDataMapper()
cut_mapper.SetInputConnection(cutter.GetOutputPort())
cut_actor = vtk.vtkActor()
cut_actor.SetMapper(cut_mapper)

# Set up the renderer and interactor
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.AddActor(cut_actor)
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# Set up the plane widget for interactive cutting
plane_widget = vtk.vtkPlaneWidget()
plane_widget.SetInteractor(interactor)
plane_widget.SetInputConnection(sphere.GetOutputPort())
plane_widget.PlaceWidget()
plane_widget.AddObserver("InteractionEvent", lambda obj, event: plane_widget.GetPlane(plane))

# Start the interactor
interactor.Initialize()
interactor.Start()
```
