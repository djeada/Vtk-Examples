"""
This example demonstrates the basics of volume rendering in VTK.
Volume rendering is a method used for displaying a 2D projection
of a 3D discretely sampled data set, typically a 3D scalar field.

This script creates a synthetic 3D scalar field where the scalar value
is a function of the voxel's position. The scalar field is then visualized
using volume rendering, where the color and opacity of each voxel is
determined by the scalar value and the transfer functions.
"""

import vtk


def main():
    # Create a synthetic volume
    volume_data = vtk.vtkImageData()
    volume_data.SetDimensions(50, 50, 50)
    volume_data.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)

    # Assign a scalar value to each voxel in the volume
    for z in range(50):
        for y in range(50):
            for x in range(50):
                volume_data.SetScalarComponentFromDouble(x, y, z, 0, x + y + z)

    # Define opacity transfer function
    # Opacity transfer function is a mapping of voxel scalar value to voxel opacity
    opacity_transfer_function = vtk.vtkPiecewiseFunction()
    opacity_transfer_function.AddPoint(
        0, 0.0
    )  # voxels with scalar value 0 have 0 opacity
    opacity_transfer_function.AddPoint(
        255, 1.0
    )  # voxels with scalar value 255 have 1 opacity

    # Define color transfer function
    color_transfer_function = vtk.vtkColorTransferFunction()
    color_transfer_function.AddRGBPoint(
        0.0, 0.0, 0.0, 1.0
    )  # voxels with scalar value 0 are blue
    color_transfer_function.AddRGBPoint(
        255.0, 1.0, 0.0, 0.0
    )  # voxels with scalar value 255 are red

    # Define volume properties using the transfer functions
    volume_properties = vtk.vtkVolumeProperty()
    volume_properties.SetColor(color_transfer_function)
    volume_properties.SetScalarOpacity(opacity_transfer_function)
    volume_properties.ShadeOn()  # enable shading

    # Define volume mapper
    # Volume mapper takes a 3D volume and maps it to 2D screen space
    volume_mapper = vtk.vtkSmartVolumeMapper()
    volume_mapper.SetInputData(volume_data)

    # Define volume actor
    # Volume actor is used to set the position and orientation of the volume
    volume_actor = vtk.vtkVolume()
    volume_actor.SetMapper(volume_mapper)
    volume_actor.SetProperty(volume_properties)

    # Define renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window_interactor = vtk.vtkRenderWindowInteractor()

    render_window.AddRenderer(renderer)
    render_window_interactor.SetRenderWindow(render_window)

    # Add volume actor to the renderer
    renderer.AddVolume(volume_actor)

    # Start the rendering loop
    render_window_interactor.Initialize()
    render_window.Render()
    render_window.SetWindowName("Volume Rendering - Cube")
    render_window_interactor.Start()


if __name__ == "__main__":
    main()
