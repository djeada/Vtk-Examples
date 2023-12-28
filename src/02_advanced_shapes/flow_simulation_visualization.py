from typing import Tuple

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


def create_vector_field(grid_size: int = 20) -> np.ndarray:
    """
    Create a synthetic vector field representing airflow.

    Args:
        grid_size (int): The size of the grid for the vector field.

    Returns:
        np.ndarray: A 3D numpy array representing the vector field.
    """
    spacing = 0.1  # Spacing between points in the grid
    x, y, z = (
        np.linspace(-1, 1, grid_size),
        np.linspace(-1, 1, grid_size),
        np.linspace(-1, 1, grid_size),
    )
    x, y, z = np.meshgrid(x, y, z, indexing="ij")

    # Create a rotational vector field
    u = -np.sin(np.pi * y) * np.cos(np.pi * z)
    v = np.sin(np.pi * x) * np.cos(np.pi * z)
    w = np.zeros_like(u)

    return np.stack((u, v, w), axis=-1)


def create_structured_grid(vectors: np.ndarray) -> vtk.vtkStructuredGrid:
    """
    Convert a numpy vector field to a VTK structured grid.

    Args:
        vectors (np.ndarray): A 3D numpy array representing the vector field.

    Returns:
        vtk.vtkStructuredGrid: A VTK structured grid.
    """
    grid_size = vectors.shape[0]
    vtk_vectors = numpy_to_vtk(vectors.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_vectors.SetNumberOfComponents(3)
    vtk_vectors.SetName("Vectors")

    structured_grid = vtk.vtkStructuredGrid()
    structured_grid.SetDimensions(grid_size, grid_size, grid_size)

    points = vtk.vtkPoints()
    for i in range(grid_size):
        for j in range(grid_size):
            for k in range(grid_size):
                points.InsertNextPoint(i, j, k)

    structured_grid.SetPoints(points)
    structured_grid.GetPointData().SetVectors(vtk_vectors)
    return structured_grid


def visualize_vector_field(structured_grid: vtk.vtkStructuredGrid):
    """
    Visualize the given vector field using VTK.

    Args:
        structured_grid (vtk.vtkStructuredGrid): A VTK structured grid representing the vector field.
    """
    glyph_source = vtk.vtkArrowSource()
    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(glyph_source.GetOutputPort())
    glyph.SetInputData(structured_grid)
    glyph.SetScaleFactor(0.05)

    glyph_mapper = vtk.vtkPolyDataMapper()
    glyph_mapper.SetInputConnection(glyph.GetOutputPort())

    glyph_actor = vtk.vtkActor()
    glyph_actor.SetMapper(glyph_mapper)

    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    renderer.AddActor(glyph_actor)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    render_window.Render()
    interactor.Start()


def main():
    vectors = create_vector_field()
    structured_grid = create_structured_grid(vectors)
    visualize_vector_field(structured_grid)


if __name__ == "__main__":
    main()
