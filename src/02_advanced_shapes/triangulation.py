"""
This module demonstrates triangulation in VTK. Triangulation is a process
by which a surface or an object is divided into triangles. This is done
because a triangle is a simple shape that is easy to render, and a complex
shape can be represented as a collection of triangles.

This example generates a set of random points and uses the vtkDelaunay2D
filter to create a triangulation from those points. In addition, it also
shows the edges of the triangulation.
"""

import random

import vtk

# create a point cloud
points = vtk.vtkPoints()

# generate 25 random points
for _ in range(25):
    points.InsertNextPoint(random.uniform(-1, 1), random.uniform(-1, 1), 0)

# create a polydata object
input_polydata = vtk.vtkPolyData()
input_polydata.SetPoints(points)

# create a delaunay2D filter and set the input polydata
delaunay = vtk.vtkDelaunay2D()
delaunay.SetInputData(input_polydata)
delaunay.Update()

# create a triangle filter
triangle_filter = vtk.vtkTriangleFilter()
triangle_filter.SetInputConnection(delaunay.GetOutputPort())
triangle_filter.PassLinesOn()
triangle_filter.PassVertsOn()
triangle_filter.Update()

# create a mapper and actor for the triangulation
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(triangle_filter.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetRepresentationToWireframe()

# create a renderer and add the actor to it
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(0, 0, 0)

# create a render window and add the renderer to it
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

# create an interactor and connect it to the render window
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# initialize the interactor and start the rendering loop
interactor.Initialize()
interactor.Start()
