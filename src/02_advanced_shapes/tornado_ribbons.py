"""A stylized twisting tornado made from helices, not CFD."""
import math

import vtk


def make_tornado(strands=18, samples=170):
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    colors = vtk.vtkFloatArray()
    colors.SetName("height")
    for strand in range(strands):
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(samples)
        phase = 2 * math.pi * strand / strands
        for i in range(samples):
            t = i / (samples - 1)
            theta = phase + 10 * math.pi * t + 0.6 * math.sin(8 * t + phase)
            radius = 0.12 + 1.45 * t**1.35
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            z = 5 * t
            idx = points.InsertNextPoint(x, y, z)
            line.GetPointIds().SetId(i, idx)
            colors.InsertNextValue(t)
        lines.InsertNextCell(line)
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetLines(lines)
    mesh.GetPointData().SetScalars(colors)
    return mesh


def build_scene():
    tubes = vtk.vtkTubeFilter()
    tubes.SetInputData(make_tornado())
    tubes.SetRadius(0.022)
    tubes.SetNumberOfSides(7)
    colors = vtk.vtkColorTransferFunction()
    colors.AddRGBPoint(0.0, 0.11, 0.42, 0.95)
    colors.AddRGBPoint(0.48, 0.1, 0.92, 0.93)
    colors.AddRGBPoint(1.0, 1.0, 0.38, 0.14)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tubes.GetOutputPort())
    mapper.SetLookupTable(colors)
    mapper.SetScalarRange(0, 1)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.012, 0.018, 0.045)
    renderer.AddActor(actor)
    return renderer, lambda: actor.RotateZ(1.0)


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Tornado Ribbons")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.25)

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
