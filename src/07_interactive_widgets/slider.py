import vtk


def main():
    # Create a sphere
    sphereSource = vtk.vtkSphereSource()
    sphereSource.Update()

    # Create a mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphereSource.GetOutputPort())

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Create a renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)

    # Create a render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Create a render window interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Initialize the interactor and start the rendering loop
    renderWindow.Render()
    renderWindowInteractor.Initialize()

    # Setup picker
    picker = vtk.vtkCellPicker()

    # Create a text actor to display the cell ID
    textActor = vtk.vtkTextActor()
    textActor.SetInput("Cell ID: ")
    textActor.GetTextProperty().SetColor(1.0, 1.0, 1.0)  # White color
    textActor.SetPosition(10, 10)  # Position in screen coordinates
    renderer.AddActor(textActor)

    def onClick(event, obj):
        x, y = renderWindowInteractor.GetEventPosition()
        picker.Pick(x, y, 0, renderer)
        pickedCell = picker.GetCellId()
        if pickedCell != -1:
            textActor.SetInput(f"Cell ID: {pickedCell}")
            renderWindow.Render()

    # Add click event listener
    renderWindowInteractor.AddObserver("LeftButtonPressEvent", onClick)

    # Start interaction
    renderWindowInteractor.Start()


if __name__ == "__main__":
    main()
