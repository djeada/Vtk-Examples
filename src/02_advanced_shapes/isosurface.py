"""
This module demonstrates the process of creating and visualizing a 3D scalar dataset using the Visualization Toolkit (VTK) and NumPy. It is structured into several key functions, each performing a distinct step in the visualization workflow.

Workflow:
1. Creation of Sample Data (`create_sample_data`):
   - Generates a 3D scalar dataset based on a mathematical function.
   - Utilizes NumPy for efficient computation of scalar values within the dataset.
   - The dataset is represented as VTK's vtkImageData, which stores the 3D scalar values.

2. Isosurface Creation (`create_isosurface`):
   - Generates an isosurface from the vtkImageData.
   - An isosurface represents a 3D surface that connects points of equal scalar value (isovalue).
   - This function takes a range of isovalues and creates corresponding contours in the dataset.

3. Color Mapping (`apply_color_map`):
   - Applies a color map to the isosurface actor.
   - The color mapping is based on the scalar values, enhancing the visualization of the isosurface.
   - This step involves creating a vtkColorTransferFunction that maps scalar values to specific colors.

4. Visualization (`visualize_actor`):
   - Sets up a VTK renderer, render window, and interactor to display the isosurface.
   - The isosurface actor, created and colored in previous steps, is added to the VTK renderer for visualization.
   - The background color and other rendering settings are configured for optimal display.
"""

from typing import List, Tuple

import numpy as np
import vtk

# Constants
BACKGROUND_COLOR = (0.1, 0.2, 0.3)  # Background color for the visualization


def create_sample_data(dimensions: Tuple[int, int, int]) -> vtk.vtkImageData:
    """
    Create a sample 3D scalar dataset using a mathematical function.

    Args:
        dimensions (Tuple[int, int, int]): Dimensions of the dataset (x, y, z).

    Returns:
        vtk.vtkImageData: The created 3D scalar dataset.
    """
    if not all(d > 0 for d in dimensions):
        raise ValueError("Dimensions must be positive integers")

    # Create an empty image data object
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(dimensions)
    imageData.SetSpacing(1.0, 1.0, 1.0)
    imageData.SetOrigin(0, 0, 0)

    # Efficiently fill the image data with scalar values using NumPy
    x, y, z = np.indices(dimensions)
    values = np.sin(np.sqrt(x**2 + y**2 + z**2)).ravel()

    # Convert numpy array to VTK array
    scalars = vtk.vtkFloatArray()
    scalars.SetNumberOfValues(len(values))
    scalars.SetVoidArray(values, len(values), 1)

    imageData.GetPointData().SetScalars(scalars)
    return imageData


def create_isosurface(
    data: vtk.vtkImageData, isovalue_range: List[float]
) -> vtk.vtkActor:
    """
    Create an isosurface with a range of isovalues from the given data.

    Args:
        data (vtk.vtkImageData): The input dataset.
        isovalue_range (List[float]): Range of isovalues.

    Returns:
        vtk.vtkActor: The actor representing the isosurface.
    """
    # Generate an isosurface
    contour = vtk.vtkContourFilter()
    contour.SetInputData(data)
    contour.GenerateValues(len(isovalue_range), isovalue_range[0], isovalue_range[-1])

    # Map the contours to graphical primitives
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(contour.GetOutputPort())
    mapper.SetScalarRange(isovalue_range[0], isovalue_range[-1])

    # Create an actor for the contours
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


def apply_color_map(actor: vtk.vtkActor) -> None:
    """
    Apply a color map to the actor.

    Args:
        actor (vtk.vtkActor): The actor to apply the color map to.
    """
    # Create a color transfer function
    colorTransferFunction = vtk.vtkColorTransferFunction()
    colorTransferFunction.AddRGBPoint(-1.0, 0, 0, 1)  # Blue for low values
    colorTransferFunction.AddRGBPoint(0, 0, 1, 0)  # Green for mid values
    colorTransferFunction.AddRGBPoint(1.0, 1, 0, 0)  # Red for high values

    # Apply the color transfer function to the actor
    actor.GetMapper().SetLookupTable(colorTransferFunction)


def visualize_actor(actor: vtk.vtkActor) -> None:
    """
    Visualize the given actor in a render window.

    Args:
        actor (vtk.vtkActor): The actor to visualize.
    """
    # Create a renderer and render window
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Create a render window interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Add the actor to the scene
    renderer.AddActor(actor)
    renderer.SetBackground(*BACKGROUND_COLOR)

    # Start the visualization
    renderWindow.Render()
    renderWindowInteractor.Start()


# Main script
if __name__ == "__main__":
    try:
        data = create_sample_data((50, 50, 50))
        iso_actor = create_isosurface(
            data, np.linspace(-1, 1, 10).tolist()
        )  # Isovalue range
        apply_color_map(iso_actor)
        visualize_actor(iso_actor)
    except Exception as e:
        print(f"Error: {e}")
