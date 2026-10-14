"""A recursive arrangement of cubes combined into one VTK mesh."""
import vtk


def make_sculpture(depth=3):
    cube = vtk.vtkCubeSource()
    append = vtk.vtkAppendPolyData()

    def add_branch(center, size, level):
        transform = vtk.vtkTransform()
        transform.PostMultiply()
        transform.Scale(size, size, size)
        transform.Translate(*center)
        instance = vtk.vtkTransformPolyDataFilter()
        instance.SetTransform(transform)
        instance.SetInputConnection(cube.GetOutputPort())
        append.AddInputConnection(instance.GetOutputPort())
        if level >= depth:
            return
        distance = size * 0.83
        for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                           (0, -1, 0), (0, 0, 1), (0, 0, -1)):
            child = (center[0] + dx * distance,
                     center[1] + dy * distance,
                     center[2] + dz * distance)
            add_branch(child, size * 0.43, level + 1)

    add_branch((0, 0, 0), 1.8, 0)
    append.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(append.GetOutput())
    return mesh


def build_scene():
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(make_sculpture())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.28, 0.72, 1.0)
    actor.GetProperty().SetSpecular(0.55)
    actor.GetProperty().SetSpecularPower(28)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.015, 0.018, 0.04)
    return renderer, lambda: actor.RotateY(0.5)


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Fractal Sculpture")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(25)
    renderer.GetActiveCamera().Elevation(25)
    renderer.ResetCameraClippingRange()
    renderer.GetActiveCamera().Zoom(1.1)

    def tick(obj, event):
        step()
        window.Render()

    interactor.AddObserver("TimerEvent", tick)
    window.Render()
    interactor.Initialize()
    interactor.CreateRepeatingTimer(35)
    interactor.Start()


if __name__ == "__main__":
    main()
