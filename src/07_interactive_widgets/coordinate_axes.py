import random

import vtk

if __name__ == "__main__":

    plane1 = vtk.vtkPlaneSource()
    plane1.SetOrigin(0, 0, 0)
    plane1.SetPoint1(0, 1, 0)
    plane1.SetPoint1(1, 0, 0)
    plane1.SetNormal(1.0, 0.0, 0.0)
    plane1.Update()

    plane2 = vtk.vtkPlaneSource()
    plane2.SetCenter(0, 1.0, 0)
    plane2.SetNormal(0.0, 1.0, 0.0)
    plane2.Update()

    sources = (plane1, plane2)

    mappers = []
    for source in sources:
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(source.GetOutput())
        mappers.append(mapper)

    actors = []
    for mapper in mappers:
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(random.random(), random.random(), random.random())
        actors.append(actor)

    renderer = vtk.vtkRenderer()
    for actor in actors:
        renderer.AddActor(actor)

    renderer.ResetCamera()

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(800, 600)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    interactor.SetRenderWindow(window)
    interactor.Initialize()

    # add axes widget
    axes = vtk.vtkAxesActor()
    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(interactor)
    widget.SetViewport(0.0, 0.0, 0.4, 0.4)
    widget.SetEnabled(1)
    widget.InteractiveOn()

    interactor.Start()
