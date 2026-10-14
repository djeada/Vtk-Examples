"""Cut a cylindrical tunnel through a sphere using implicit CSG."""
import vtk


def make_difference(resolution=72):
    sphere = vtk.vtkSphere()
    sphere.SetRadius(1.0)
    cylinder = vtk.vtkCylinder()
    cylinder.SetRadius(0.36)
    difference = vtk.vtkImplicitBoolean()
    difference.SetOperationTypeToDifference()
    difference.AddFunction(sphere)
    difference.AddFunction(cylinder)
    sampler = vtk.vtkSampleFunction()
    sampler.SetImplicitFunction(difference)
    sampler.SetModelBounds(-1.15, 1.15, -1.15, 1.15, -1.15, 1.15)
    sampler.SetSampleDimensions(resolution, resolution, resolution)
    sampler.ComputeNormalsOff()
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputConnection(sampler.GetOutputPort())
    contour.SetValue(0, 0)
    contour.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(contour.GetOutput())
    return mesh


def build_scene():
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(make_difference())
    normals.SetFeatureAngle(55)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.84, 0.34, 0.14)
    actor.GetProperty().SetSpecular(0.5)
    actor.GetProperty().SetSpecularPower(28)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.018, 0.025, 0.05)
    return renderer, lambda: actor.RotateY(0.65)


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Boolean Geometry")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30)
    renderer.GetActiveCamera().Elevation(20)
    renderer.ResetCameraClippingRange()
    renderer.GetActiveCamera().Zoom(1.35)

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
