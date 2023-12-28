import vtk

from src.common.button import Button2D


def on_button_clicked():
    print("Button clicked!")


def main():
    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1, 1, 1)  # White background
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(300, 300)
    render_window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Create a 2D button
    button = Button2D(interactor, "Click Me", (50, 50), (100, 50))
    button.add_click_callback(on_button_clicked)

    # Start the interaction
    render_window.Render()
    interactor.Initialize()
    interactor.Start()


if __name__ == "__main__":
    main()
