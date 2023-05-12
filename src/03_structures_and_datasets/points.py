"""
n VTK, a point is a location in space. Points are used to define the geometry of graphical primitives. You can think of points as the vertices of graphical primitives such as lines, polygons, and volumes.

To store and manipulate points, VTK provides the vtkPoints class. This class stores an array of 3D points, and provides methods for inserting points, retrieving points, and other point-related operations.
"""
import vtk
import vtk.util.numpy_support as vtk_np
import numpy as np

from src.simple_pipeline import VisualisationPipeline

# Create a vtkPoints object.
points = vtk.vtkPoints()

# Insert some points.
points.InsertNextPoint(0.0, 0.0, 0.0)
points.InsertNextPoint(1.0, 0.0, 0.0)
points.InsertNextPoint(0.0, 1.0, 0.0)

# Print the number of points.
print("Number of points:", points.GetNumberOfPoints())

# Retrieve the points.
for i in range(points.GetNumberOfPoints()):
    point = points.GetPoint(i)
    print(f"Point {i}: {point}")

# Convert the vtkPoints object to a numpy array.
numpy_points = vtk_np.vtk_to_numpy(points.GetData())
print("Numpy array:\n", numpy_points)

# Convert a numpy array to vtkPoints.
new_numpy_points = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
new_vtk_points = vtk.vtkPoints()
new_vtk_points.SetData(vtk_np.numpy_to_vtk(new_numpy_points))
print("New vtkPoints:")
for i in range(new_vtk_points.GetNumberOfPoints()):
    point = new_vtk_points.GetPoint(i)
    print(f"Point {i}: {point}")

# Create polydata to store everything
polydata = vtk.vtkPolyData()

# Add the points to the dataset
polydata.SetPoints(points)

# Glyph the points to make them visible
glyphFilter = vtk.vtkVertexGlyphFilter()
glyphFilter.SetInputData(polydata)
glyphFilter.Update()

# Create a mapper and actor
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(glyphFilter.GetOutputPort())

pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
pipeline.run()
