"""Procedural lava planet: displace a sphere and color its elevation."""
import math

import vtk


def make_planet():
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1.0)
    sphere.SetThetaResolution(110)
    sphere.SetPhiResolution(70)
    sphere.Update()

    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(sphere.GetOutput())
    points = mesh.GetPoints()
    elevation = vtk.vtkFloatArray()
    elevation.SetName("elevation")
    for i in range(points.GetNumberOfPoints()):
        x, y, z = points.GetPoint(i)
        noise = (
            0.52 * math.sin(9 * x + 5 * y) * math.cos(11 * z - 3 * x)
            + 0.30 * math.sin(22 * x - 13 * z) * math.cos(16 * y)
            + 0.18 * math.cos(37 * y + 29 * z)
        )
        height = 0.10 * noise
        points.SetPoint(i, x * (1 + height), y * (1 + height), z * (1 + height))
        elevation.InsertNextValue(height)
    mesh.GetPointData().SetScalars(elevation)
    mesh.Modified()
    return mesh


def build_scene():
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(make_planet())
    normals.SplittingOff()
    color = vtk.vtkColorTransferFunction()
    color.AddRGBPoint(-0.10, 0.06, 0.025, 0.025)
    color.AddRGBPoint(-0.025, 0.32, 0.045, 0.012)
    color.AddRGBPoint(0.02, 0.94, 0.17, 0.025)
    color.AddRGBPoint(0.10, 1.0, 0.8, 0.18)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.SetLookupTable(color)
    mapper.SetScalarRange(-0.10, 0.10)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.25)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.015, 0.01, 0.035)
    renderer.AddActor(actor)
    return renderer, lambda: actor.RotateY(0.7)


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Lava Planet")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.4)

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
