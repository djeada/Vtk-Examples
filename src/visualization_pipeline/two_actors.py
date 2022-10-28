import random

import vtk


def setup_data():
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(0, 0, 0)
    sphere.SetRadius(5)
    sphere.SetThetaResolution(100)
    return sphere


def setup_mapper(sphere):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    return mapper


def setup_actor(mapper, pos):
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(pos)
    actor.GetProperty().SetColor(random.random(), random.random(), random.random())
    return actor


def setup_renderer(actors):
    renderer = vtk.vtkRenderer()
    for actor in actors:
        renderer.AddActor(actor)
    renderer.SetBackground(0, 0, 0)
    return renderer


def setup_window(renderer):
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(800, 600)
    return window


def setup_interactor(window):
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.Initialize()
    return interactor


if __name__ == "__main__":
    sphere1 = setup_data()
    mapper1 = setup_mapper(sphere1)
    actor1 = setup_actor(mapper1, (0, 0, 0))

    sphere2 = setup_data()
    mapper2 = setup_mapper(sphere2)
    actor2 = setup_actor(mapper2, (10, 0, 0))

    renderer = setup_renderer([actor1, actor2])
    window = setup_window(renderer)
    interactor = setup_interactor(window)

    interactor.Start()
