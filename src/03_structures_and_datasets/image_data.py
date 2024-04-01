from PIL import Image
import numpy as np
import vtk

def create_gradient_image(width, height, filename='gradient_image.png'):
    """
    Create a gradient image and save it.
    """
    # Create a gradient from blue to green
    gradient = np.zeros((height, width, 3), dtype=np.uint8)
    for i in range(height):
        color = int(255 * i / height)  # Gradient factor
        gradient[i, :, 0] = 255 - color
        gradient[i, :, 1] = color

    image = Image.fromarray(gradient, 'RGB')
    image.save(filename)
    return filename
def load_image_as_vtk_image_data(filename):
    """
    Load an image and convert it into vtkImageData.
    """
    reader = vtk.vtkPNGReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()

def visualize_vtk_image_data(vtk_image_data):
    """
    Visualize vtkImageData.
    """
    # Map the image data through a lookup table
    color = vtk.vtkImageMapToColors()
    color.SetOutputFormatToRGB()
    color.SetInputData(vtk_image_data)

    # Create an actor
    actor = vtk.vtkImageActor()
    actor.GetMapper().SetInputConnection(color.GetOutputPort())

    # Create a renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0, 0, 0)  # Set black background

    # Create a render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Create a render window interactor
    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Start the interaction
    render_window.Render()
    render_window_interactor.Start()

def main():
    filename = create_gradient_image(256, 256)
    vtk_image_data = load_image_as_vtk_image_data(filename)
    visualize_vtk_image_data(vtk_image_data)

if __name__ == "__main__":
    main()
