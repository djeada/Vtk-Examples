"""
This module demonstrates the concept of Isosurfaces in VTK.

Isosurfaces are three-dimensional analogs of contour lines. They represent
points of a constant value (the isovalue) within a volume of space. In
other words, an isosurface is a surface that represents points in the 3D
data volume that are all the same value.
"""

import vtk

# Create a grid
grid = vtk.vtkImageData()
grid.SetDimensions(25, 25, 25)
grid.SetOrigin(-12.0, -12.0, -12.0)
grid.SetSpacing(1.0, 1.0, 1.0)

# Create some scalars
scalars = vtk.vtkFloatArray()
scalars.SetNumberOfValues(25 * 25 * 25)

for i in range(25 * 25 * 25):
    scalars.SetValue(i, i % 255)

grid.GetPointData().SetScalars(scalars)

# Generate an isosurface at a particular value
contours = vtk.vtkContourFilter()
contours.SetInputData(grid)
contours.GenerateValues(1, 125, 125)

# Map the contours to graphical primitives
contourMapper = vtk.vtkPolyDataMapper()
contourMapper.SetInputConnection(contours.GetOutputPort())
contourMapper.SetScalarRange(0, 255)

# Create an actor for the contours
contourActor = vtk.vtkActor()
contourActor.SetMapper(contourMapper)

# Visualize
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(contourActor)
renderer.SetBackground(0, 0, 0)  # Background color black

renderWindowInteractor.Initialize()
renderWindow.Render()
renderWindowInteractor.Start()
