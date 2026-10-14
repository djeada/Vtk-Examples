"""Illustrative gravity-well grid, not a relativistic simulation."""
import math

import vtk


def well_height(x, y):
    return -1.45 / math.sqrt(x * x + y * y + 0.35)


def make_well():
    plane = vtk.vtkPlaneSource()
    plane.SetOrigin(-4, -4, 0)
    plane.SetPoint1(4, -4, 0)
    plane.SetPoint2(-4, 4, 0)
    plane.SetXResolution(58)
    plane.SetYResolution(58)
    plane.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(plane.GetOutput())
    points = mesh.GetPoints()
    for i in range(points.GetNumberOfPoints()):
        x, y, _ = points.GetPoint(i)
        points.SetPoint(i, x, y, well_height(x, y))
    points.Modified()
    return mesh


def build_scene():
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(make_well())
    grid = vtk.vtkActor()
    grid.SetMapper(mapper)
    grid.GetProperty().SetRepresentationToWireframe()
    grid.GetProperty().SetColor(0.19, 0.77, 0.92)
    grid.GetProperty().SetOpacity(0.85)
    ball = vtk.vtkSphereSource()
    ball.SetRadius(0.23)
    ball.SetThetaResolution(28)
    ball.SetPhiResolution(20)
    ball_mapper = vtk.vtkPolyDataMapper()
    ball_mapper.SetInputConnection(ball.GetOutputPort())
    satellite = vtk.vtkActor()
    satellite.SetMapper(ball_mapper)
    satellite.GetProperty().SetColor(1, 0.63, 0.15)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(0.018, 0.025, 0.05)
    renderer.AddActor(grid)
    renderer.AddActor(satellite)
    camera = renderer.GetActiveCamera()
    camera.SetPosition(9, -12, 9)
    camera.SetFocalPoint(0, 0, -1)
    camera.SetViewUp(0, 0, 1)
    frame = [0]

    def step():
        frame[0] += 1
        theta = frame[0] * 0.025
        x, y = 2.2 * math.cos(theta), 2.2 * math.sin(theta)
        satellite.SetPosition(x, y, well_height(x, y) + 0.35)

    step()
    return renderer, step


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Gravity Well (Illustration)")
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
