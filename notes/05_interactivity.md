## Interactivity and User Interface

VTK comes equipped with a range of tools designed to help developers create interactive visualizations and user interfaces. Some of the popular techniques employed for this purpose include:

1. Picking and Selection
2. Cutting Planes
3. 3D Widgets
4. Event Handling and Callbacks

## Picking and Selection

Picking refers to the process of interactively selecting a particular object or point within the visualization. Selection, on the other hand, involves highlighting or extracting a subset of data based on certain defined criteria.

Here are two of the many classes VTK provides for this purpose:

- `vtkCellPicker`: This class allows you to select cells or points within a 3D scene.
- `vtkExtractSelection`: This class helps to extract a subset of data based on the selection criteria.

For instance, to pick a cell within a 3D scene, you would use the `vtkCellPicker` class as follows:

```python
picker = vtk.vtkCellPicker()
picker.Pick(10, 10, 0, renderer)
pickedCell = picker.GetCellId()
```

![picker](https://github.com/djeada/Vtk-Examples/assets/37275728/51f7b6ed-9086-4fae-b976-8ac0dc139be2)

## Cutting Planes

Cutting planes are used to interactively slice through a 3D object, revealing its internal structure. This technique is particularly useful when you need to inspect the inner structure of complex objects.

The two important classes involved in this process are:

- `vtkPlane`: This class is used to define a cutting plane.
- `vtkCutter`: This class uses a vtkPlane to slice through an object.

To illustrate, if we wanted to create a cutting plane through a 3D object, we might use:

```python
plane = vtk.vtkPlane()
plane.SetNormal(1, 0, 0)

cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputData(object)
```

![plane_cutter](https://github.com/djeada/Vtk-Examples/assets/37275728/8d6e1588-f198-4ccf-b920-8a6762ea92ce)

## 3D Widgets

3D widgets are interactive tools used for manipulating data or objects within a 3D scene. They can be used for tasks such as positioning lights or cameras, or adjusting clipping planes.

Some common 3D widget classes include:

- `vtkBoxWidget`: This widget allows for interactive resizing of a bounding box.
- `vtkPlaneWidget`: This widget is used to manipulate a plane within a 3D scene.
- `vtkSliderWidget`: This widget is used to adjust scalar values via a slider.

For instance, a vtkSliderWidget could be used to control the opacity of an actor like this:

```python
sliderWidget = vtk.vtkSliderWidget()
sliderWidget.SetInteractor(interactor)
sliderWidget.SetRepresentation(sliderRep)
sliderWidget.AddObserver("InteractionEvent", lambda obj, event: actor.GetProperty().SetOpacity(obj.GetRepresentation().GetValue()))
```

![slider](https://github.com/djeada/Vtk-Examples/assets/37275728/bf595429-3a41-4e47-b996-36d216138cde)

## Event Handling and Callbacks

Event handling involves responding to user input or other events within a VTK application, whereas callbacks are user-defined functions that get executed in response to specific events.

Key classes for this purpose include:

- `vtkRenderWindowInteractor`: This class manages user input and event handling.
- `vtkCallbackCommand`: This class associates a user-defined callback function with a specific event.

### Mouse Events

Mouse events can be captured and used to control aspects of the VTK visualization. The `vtkRenderWindowInteractor` class provides the functionality to handle mouse events.

Here is an example of handling a mouse click event:

```python
def onMouseClick(obj, event):
    x, y = obj.GetEventPosition()
    print(f"Mouse clicked at coordinates {x}, {y}")

interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
interactorStyle.AddObserver("LeftButtonPressEvent", onMouseClick)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetInteractorStyle(interactorStyle)
```

### Keyboard Input Events

Similarly, keyboard input events can also be captured using the vtkRenderWindowInteractor class. These events can be used to interact with the visualization or control the program flow.

Here is an example of handling a keyboard press event:

```python
def onKeyPress(obj, event):
    key = obj.GetKeySym()
    print(f"Key {key} was pressed")

interactorStyle = vtk.vtkInteractorStyleTrackballCamera()
interactorStyle.AddObserver("KeyPressEvent", onKeyPress)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetInteractorStyle(interactorStyle)
```
