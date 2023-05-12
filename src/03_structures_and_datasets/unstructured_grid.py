import vtk
import numpy as np

from src.simple_pipeline import VisualisationPipeline

# Create a hexahedron
hexahedron = vtk.vtkHexahedron()
hexahedron.GetPointIds().SetId(0, 0)
hexahedron.GetPointIds().SetId(1, 1)
hexahedron.GetPointIds().SetId(2, 2)
hexahedron.GetPointIds().SetId(3, 3)
hexahedron.GetPointIds().SetId(4, 4)
hexahedron.GetPointIds().SetId(5, 5)
hexahedron.GetPointIds().SetId(6, 6)
hexahedron.GetPointIds().SetId(7, 7)

# Create a tetrahedron
tetrahedron = vtk.vtkTetra()
tetrahedron.GetPointIds().SetId(0, 8)
tetrahedron.GetPointIds().SetId(1, 9)
tetrahedron.GetPointIds().SetId(2, 10)
tetrahedron.GetPointIds().SetId(3, 11)

# Create the points for the two cells
points = vtk.vtkPoints()
points.InsertNextPoint(0, 0, 0)
points.InsertNextPoint(1, 0, 0)
points.InsertNextPoint(1, 1, 0)
points.InsertNextPoint(0, 1, 0)
points.InsertNextPoint(0, 0, 1)
points.InsertNextPoint(1, 0, 1)
points.InsertNextPoint(1, 1, 1)
points.InsertNextPoint(0, 1, 1)
points.InsertNextPoint(2, 2, 2)
points.InsertNextPoint(2, 3, 2)
points.InsertNextPoint(3, 2, 2)
points.InsertNextPoint(2, 2, 3)

# Create an unstructured grid
ugrid = vtk.vtkUnstructuredGrid()
ugrid.SetPoints(points)
ugrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())
ugrid.InsertNextCell(tetrahedron.GetCellType(), tetrahedron.GetPointIds())

# Attach a scalar field to the points
scalars = vtk.vtkFloatArray()
scalars.SetName("Scalars")
for i in range(points.GetNumberOfPoints()):
    scalars.InsertNextValue(np.random.rand())
ugrid.GetPointData().SetScalars(scalars)

# Attach a vector field to the cells
vectors = vtk.vtkFloatArray()
vectors.SetNumberOfComponents(3)
vectors.SetName("Vectors")
for i in range(ugrid.GetNumberOfCells()):
    vectors.InsertNextTuple3(np.random.rand(), np.random.rand(), np.random.rand())
ugrid.GetCellData().SetVectors(vectors)

# Visualizing the fields
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(ugrid)
mapper.SetScalarModeToUsePointData()
mapper.SetColorModeToMapScalars()
mapper.SelectColorArray("Scalars")


pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
pipeline.run()
