"""
Fields in VTK can be attached to both points and cells. In this example, I will demonstrate how to attach a scalar field to points and a vector field to cells. We can create a simple situation where we assign random scalar values to points and random vector values to cells.
in VTK, you can store physical properties such as pressure, temperature, or any scalar or vector quantity associated with the points or cells in your dataset. These are often referred to as point data or cell data, respectively.

    Point Data: These are data attributes that are associated with the points of a dataset. For example, in a computational fluid dynamics (CFD) simulation, you might store the velocity vector at each point in the mesh as point data.

    Cell Data: These are data attributes that are associated with the cells of a dataset. For example, you might store the pressure or temperature within each cell of the mesh as cell data.
"""
import vtk
import numpy as np

from src.simple_pipeline import VisualisationPipeline

# Create the points
points = vtk.vtkPoints()
points.InsertNextPoint(0, 0, 0)  # Point 0
points.InsertNextPoint(1, 0, 0)  # Point 1
points.InsertNextPoint(0, 1, 0)  # Point 2
points.InsertNextPoint(1, 1, 0)  # Point 3

# Create a quad cell
quad = vtk.vtkQuad()
quad.GetPointIds().SetId(0, 0)  # the first point of the quad is point 0
quad.GetPointIds().SetId(1, 1)  # the second point is point 1
quad.GetPointIds().SetId(2, 3)  # the third point is point 3
quad.GetPointIds().SetId(3, 2)  # the fourth point is point 2

# Create a cell array and add the cells to it
cells = vtk.vtkCellArray()
cells.InsertNextCell(quad)

# Create polydata to hold the points and cells
polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.SetPolys(cells)

# Attach a scalar field to the points
scalars = vtk.vtkFloatArray()
scalars.SetName("Scalars")  # This name will be used later for coloring
for i in range(points.GetNumberOfPoints()):
    scalars.InsertNextValue(np.random.rand())
polydata.GetPointData().SetScalars(scalars)

# Attach a vector field to the cells
vectors = vtk.vtkFloatArray()
vectors.SetNumberOfComponents(3)  # We have 3D vectors
vectors.SetName("Vectors")  # This name will be used later for coloring
for i in range(cells.GetNumberOfCells()):
    vectors.InsertNextTuple3(np.random.rand(), np.random.rand(), np.random.rand())
polydata.GetCellData().SetVectors(vectors)

# Visualizing the fields
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(polydata)
mapper.SetScalarModeToUsePointData()  # Color by the scalar data
mapper.SetColorModeToMapScalars()  # Color by the scalar data
mapper.SelectColorArray("Scalars")  # Color by the scalar data

pipeline = VisualisationPipeline(mappers=[mapper])
pipeline.run()
