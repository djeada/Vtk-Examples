import random

import vtk


def setup_cone(height, radius, resolution):
    cone = vtk.vtkConeSource()
    cone.SetHeight(height)
    cone.SetRadius(radius)
    cone.SetResolution(resolution)
    return cone


if __name__ == "__main__":

    # define sources
    cone1 = setup_cone(5, 2, 10)
    cone2 = setup_cone(3, 2, 5)
    cone3 = setup_cone(2, 2, 2)

    sources = (cone1, cone2, cone3)

    mappers = []
    for source in sources:
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        mappers.append(mapper)

    actors = []
    for mapper in mappers:
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actors.append(actor)

    renderers = []
    for i, actor in enumerate(actors):
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(random.random(), random.random(), random.random())
        renderer.SetViewport(
            0 + i * 1 / len(actors), 0.0, (i + 1) * 1 / len(actors), 1.0
        )  # x, y, width, height
        renderers.append(renderer)

    coneWin = vtk.vtkRenderWindow()
    coneWin.SetSize(800, 400)

    for renderer in renderers:
        coneWin.AddRenderer(renderer)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(coneWin)
    iren.Initialize()
    iren.Start()
