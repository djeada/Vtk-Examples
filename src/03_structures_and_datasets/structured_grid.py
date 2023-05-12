"""
Structured grids in VTK represent a mesh in which points are arranged in a regular lattice in 3D space. It is best suited for problems where the domain can be discretized as a regular grid. For example, in computational fluid dynamics, it is common to use a structured grid to represent the fluid domain.
"""

import vtk
import numpy as np

from src.simple_pipeline import VisualisationPipeline

# Define dimensions
nx, ny, nz = 6, 6, 2

# Create a structured grid
structuredGrid = vtk.vtkStructuredGrid()
structuredGrid.SetDimensions(nx, ny, nz)

# Create points
points = vtk.vtkPoints()
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(i, j, k)

# Set the points to the grid
structuredGrid.SetPoints(points)

# Create a mapper and actor for visualization
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(structuredGrid)

pipeline = VisualisationPipeline(mappers=[mapper], edges_visible=True)
pipeline.run()
