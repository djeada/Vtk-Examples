import vtk
from vtkmodules.vtkCommonCore import vtkLookupTable
from vtkmodules.vtkCommonDataModel import vtkPlane
from vtkmodules.vtkImagingColor import vtkImageMapToColors
from vtkmodules.vtkImagingCore import vtkImageReslice
from vtkmodules.vtkRenderingCore import (
    vtkImageActor,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)


class VolumeSlicer:
    def __init__(self, volume: vtk.vtkVolume):
        """
        Initialize the VolumeSlicer with a volume.

        Args:
            volume (vtk.vtkVolume): The volume to be sliced.
        """
        self.volume = volume

    def create_slicing_plane(self, normal: tuple, color: tuple) -> vtkImageActor:
        """
        Create and setup a slicing plane.

        Args:
            normal (tuple): Normal vector of the slicing plane.
            color (tuple): Color of the slicing plane.

        Returns:
            vtkImageActor: Actor representing the slice.
        """
        # Create the slicing plane
        plane = vtkPlane()
        plane.SetOrigin(self.volume.GetCenter())
        plane.SetNormal(*normal)

        # Create the reslice mapper
        reslice = vtkImageReslice()
        reslice.SetInputConnection(self.volume.GetOutputPort())
        reslice.SetOutputDimensionality(2)
        reslice.SetReslicePlane(plane)
        reslice.SetInterpolationModeToLinear()

        # Create a lookup table to map the data to colors
        lookup_table = vtkLookupTable()
        lookup_table.SetRange(0, 1000)  # Adjust according to your data
        lookup_table.SetValueRange(0.0, 1.0)
        lookup_table.SetSaturationRange(0.0, 0.0)
        lookup_table.SetRampToLinear()
        lookup_table.Build()

        # Map the image through the lookup table
        color_mapper = vtkImageMapToColors()
        color_mapper.SetLookupTable(lookup_table)
        color_mapper.SetInputConnection(reslice.GetOutputPort())

        # Create an actor to display the image
        actor = vtkImageActor()
        actor.GetMapper().SetInputConnection(color_mapper.GetOutputPort())
        actor.GetProperty().SetColor(color)

        return actor


def setup_render_window():
    """
    Set up the render window and renderer.

    Returns:
        Tuple[vtkRenderWindow, vtkRenderer]: A tuple containing the render window and renderer.
    """
    renderer = vtkRenderer()
    render_window = vtkRenderWindow()
    render_window.AddRenderer(renderer)
    return render_window, renderer


def main():
    # Load a volume dataset (replace with actual volume loading code)
    volume = vtk.vtkVolume()

    # Initialize the VolumeSlicer
    slicer = VolumeSlicer(volume)

    # Setup the slicing planes
    actor_x = slicer.create_slicing_plane((1, 0, 0), (1, 0, 0))  # Red slice in X
    actor_y = slicer.create_slicing_plane((0, 1, 0), (0, 1, 0))  # Green slice in Y
    actor_z = slicer.create_slicing_plane((0, 0, 1), (0, 0, 1))  # Blue slice in Z

    # Setup the render window and renderer
    render_window, renderer = setup_render_window()

    # Add actors to the renderer
    renderer.AddActor(actor_x)
    renderer.AddActor(actor_y)
    renderer.AddActor(actor_z)

    # Create a render window interactor and start the interaction loop
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    render_window.Render()
    interactor.Start()


if __name__ == "__main__":
    main()
