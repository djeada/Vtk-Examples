"""Generate mountain terrain from layered trigonometric height fields."""
import math

import vtk


def make_terrain():
    plane = vtk.vtkPlaneSource()
    plane.SetOrigin(-5, -5, 0)
    plane.SetPoint1(5, -5, 0)
    plane.SetPoint2(-5, 5, 0)
    plane.SetXResolution(100)
    plane.SetYResolution(100)
    plane.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(plane.GetOutput())
    heights = vtk.vtkFloatArray()
    heights.SetName("elevation")
    for x, y, _ in (mesh.GetPoint(i) for i in range(mesh.GetNumberOfPoints())):
        ridge = math.sin(x * 0.9 + math.sin(y * 0.8)) * math.cos(y * 0.65)
        detail = 0.22 * math.sin(3.5 * x + 1.7 * y) * math.cos(2.8 * y)
        peak = 1.7 * math.exp(-((x - 0.7) ** 2 + (y + 0.5) ** 2) / 4)
        heights.InsertNextValue(1.25 * ridge + detail + peak)
    mesh.GetPointData().SetScalars(heights)
    warp = vtk.vtkWarpScalar()
    warp.SetInputData(mesh)
    warp.SetScaleFactor(1.0)
    warp.UseNormalOn()
    warp.SetNormal(0, 0, 1)
    warp.Update()
    return warp.GetOutput()


def build_scene():
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(make_terrain())
    normals.SplittingOff()
    colors = vtk.vtkColorTransferFunction()
    colors.AddRGBPoint(-1.5, 0.055, 0.16, 0.25)
    colors.AddRGBPoint(0.0, 0.07, 0.34, 0.24)
    colors.AddRGBPoint(1.0, 0.57, 0.48, 0.30)
    colors.AddRGBPoint(2.5, 0.95, 0.96, 0.98)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.SetLookupTable(colors)
    mapper.SetScalarRange(-1.5, 2.5)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.03, 0.045, 0.08)
    renderer.AddActor(actor)
    renderer.GetActiveCamera().SetPosition(10, -14, 10)
    renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
    renderer.GetActiveCamera().SetViewUp(0, 0, 1)
    return renderer, lambda: actor.RotateZ(0.2)


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Procedural Terrain")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

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
