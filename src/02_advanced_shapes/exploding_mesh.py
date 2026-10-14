"""Break a triangulated sphere into independent, animated faces."""
import math

import vtk


def make_fragments():
    sphere = vtk.vtkSphereSource()
    sphere.SetThetaResolution(28)
    sphere.SetPhiResolution(20)
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputConnection(sphere.GetOutputPort())
    triangles.Update()
    source = triangles.GetOutput()
    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    origins = []
    directions = []
    for cell_index in range(source.GetNumberOfCells()):
        cell = source.GetCell(cell_index)
        if cell.GetNumberOfPoints() != 3:
            continue
        vertices = [source.GetPoint(cell.GetPointId(j)) for j in range(3)]
        center = [sum(p[j] for p in vertices) / 3 for j in range(3)]
        for vertex in vertices:
            origins.append(vertex)
            directions.append(tuple(v * 1.6 for v in center))
        polygon = vtk.vtkTriangle()
        for j in range(3):
            polygon.GetPointIds().SetId(j, points.InsertNextPoint(*vertices[j]))
        polys.InsertNextCell(polygon)
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(polys)

    def update(amount):
        for i, (origin, direction) in enumerate(zip(origins, directions)):
            points.SetPoint(i, *(origin[j] + amount * direction[j] for j in range(3)))
        points.Modified()
        mesh.Modified()

    return mesh, update


def build_scene():
    mesh, update = make_fragments()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(mesh)
    normals.SplittingOn()
    normals.SetFeatureAngle(1)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.27, 0.83, 0.98)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.03, 0.10, 0.24)
    actor.GetProperty().SetSpecular(0.3)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.01, 0.02, 0.05)
    frame = [0]

    def step():
        frame[0] += 1
        amount = 0.5 - 0.5 * math.cos(frame[0] * 0.025)
        update(amount)
        actor.RotateY(0.6)

    return renderer, step


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Exploding Mesh")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(0.8)

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
