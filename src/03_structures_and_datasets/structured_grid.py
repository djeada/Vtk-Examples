"""
Structured grids in VTK represent a mesh in which points are arranged in a regular lattice in 3D space. It is best suited for problems where the domain can be discretized as a regular grid. For example, in computational fluid dynamics, it is common to use a structured grid to represent the fluid domain.

Workflow Overview:

1. Structured Grid Creation (create_structured_grid):
   - Initializes a vtkStructuredGrid with specified dimensions (nx, ny, nz), representing the grid's resolution in each direction.
   - Populates the grid with points in a regular, lattice-like arrangement using vtkPoints, illustrating how to define spatial coordinates in a structured grid.
   - Returns a fully defined vtkStructuredGrid, ready for visualization.

2. Visualization (visualize_structured_grid):
   - Takes the created structured grid as input and sets up a vtkDataSetMapper for it. The mapper translates the grid data into a format suitable for rendering.
   - Utilizes a predefined 'VisualisationPipeline' to handle the rendering process. This pipeline is configured to make the edges of the grid visible, enhancing the perception of the grid structure in 3D space.

3. Main Function (main):
   - Defines the dimensions for the structured grid.
   - Calls the grid creation and visualization functions, orchestrating the process from grid setup to rendering.
"""

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_structured_grid(nx, ny, nz):
    """
    Create a structured grid with specified dimensions.

    Args:
    nx, ny, nz (int): Dimensions of the grid in the x, y, and z directions.

    Returns:
    vtkStructuredGrid: A structured grid with the specified dimensions.
    """
    structuredGrid = vtk.vtkStructuredGrid()
    structuredGrid.SetDimensions(nx, ny, nz)

    points = vtk.vtkPoints()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                points.InsertNextPoint(i, j, k)

    structuredGrid.SetPoints(points)
    return structuredGrid


def visualize_structured_grid(structuredGrid):
    """
    Visualize the given structured grid using a predefined visualization pipeline.

    Args:
    structuredGrid (vtkStructuredGrid): The structured grid to be visualized.
    """
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(structuredGrid)

    # Configure and run the visualization pipeline
    pipeline = VisualisationPipeline(mappers=[mapper], edges_visible=True)
    pipeline.run()


def main():
    # Define dimensions for the structured grid
    nx, ny, nz = 6, 6, 2

    # Create and visualize the structured grid
    structuredGrid = create_structured_grid(nx, ny, nz)
    visualize_structured_grid(structuredGrid)


if __name__ == "__main__":
    main()
