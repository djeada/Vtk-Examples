"""Animated plane slicing a translucent sphere with vtkCutter."""
import math

import vtk


def build_scene():
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1)
    sphere.SetThetaResolution(72)
    sphere.SetPhiResolution(56)
    plane = vtk.vtkPlane()
    plane.SetNormal(1, 0, 0)
    plane.SetOrigin(0, 0, 0)
    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(sphere.GetOutputPort())
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(cutter.GetOutputPort())
    tube.SetRadius(0.027)
    tube.SetNumberOfSides(10)
    shell_mapper = vtk.vtkPolyDataMapper()
    shell_mapper.SetInputConnection(sphere.GetOutputPort())
    shell = vtk.vtkActor()
    shell.SetMapper(shell_mapper)
    shell.GetProperty().SetColor(0.2, 0.65, 0.95)
    shell.GetProperty().SetOpacity(0.22)
    cut_mapper = vtk.vtkPolyDataMapper()
    cut_mapper.SetInputConnection(tube.GetOutputPort())
    cut = vtk.vtkActor()
    cut.SetMapper(cut_mapper)
    cut.GetProperty().SetColor(1, 0.24, 0.07)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(shell)
    renderer.AddActor(cut)
    renderer.SetBackground(0.02, 0.025, 0.055)
    frame = [0]

    def step():
        frame[0] += 1
        plane.SetOrigin(0.95 * math.sin(frame[0] * 0.025), 0, 0)
        plane.Modified()

    return renderer, step


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Slice a Sphere")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(35)
    renderer.GetActiveCamera().Elevation(15)
    renderer.ResetCameraClippingRange()
    renderer.GetActiveCamera().Zoom(1.35)

    def tick(obj, event):
        step()
        window.Render()

    interactor.AddObserver("TimerEvent", tick)
    window.Render()
    interactor.Initialize()
    interactor.CreateRepeatingTimer(30)
    interactor.Start()


if __name__ == "__main__":
    main()
