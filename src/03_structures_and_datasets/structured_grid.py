"""
Structured grids in VTK represent a mesh in which points are arranged in a regular lattice in 3D space. It is best suited for problems where the domain can be discretized as a regular grid. For example, in computational fluid dynamics, it is common to use a structured grid to represent the fluid domain.
"""

import vtk
import numpy as np

from src.simple_pipeline import VisualisationPipeline

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
