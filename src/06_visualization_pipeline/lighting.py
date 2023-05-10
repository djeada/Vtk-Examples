import vtk


def setup_data():
    cone = vtk.vtkConeSource()
    cone.SetResolution(100)
    return cone


def setup_mapper(cone):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cone.GetOutputPort())
    return mapper


def setup_actor(mapper):
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(0, 0, 0)
    actor.GetProperty().SetColor(0.8, 0.6, 0.1)
    actor.GetProperty().SetAmbient(0.1)
    actor.GetProperty().SetDiffuse(0.6)
    actor.GetProperty().SetSpecular(0.8)
    actor.GetProperty().SetSpecularPower(100)
    return actor


def setup_camera():
    camera = vtk.vtkCamera()
    camera.SetFocalPoint(0.0, 0.0, 0.0)
    camera.SetPosition(0.0, 0.0, 8.0)
    camera.SetViewUp(0.0, 1.0, 0.0)
    return camera


def setup_light():
    light = vtk.vtkLight()
    light.SetFocalPoint(0.0, 0.0, 0.0)
    light.SetPosition(0.0, 0.0, 1.0)
    light.SetColor(1.0, 1.0, 1.0)
    light.SetIntensity(1.0)
    return light


def setup_spotlight():
    spot = vtk.vtkLight()
    spot.SetFocalPoint(0.0, 0.0, 0.0)
    spot.SetPosition(0.0, 0.0, 1.0)
    spot.SetColor(1.0, 1.0, 1.0)
    spot.SetIntensity(1.0)
    spot.SetConeAngle(30)
    spot.SetPositional(True)
    spot.SetAttenuationValues(1, 0, 0)
    return spot


def setup_renderer(actor, camera, light, spot):
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetActiveCamera(camera)
    renderer.AddLight(light)
    renderer.AddLight(spot)
    renderer.SetBackground(0.0, 0.0, 0.0)
    renderer.LightFollowCameraOn()
    return renderer


def setup_render_window(renderer):
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)
    return render_window


def setup_interactor(render_window):
    istyle = vtk.vtkInteractorStyleSwitch()
    istyle.SetCurrentStyleToTrackballCamera()

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    interactor.SetInteractorStyle(istyle)

    interactor.Initialize()
    interactor.Start()


if __name__ == "__main__":
    cone = setup_data()
    mapper = setup_mapper(cone)
    actor = setup_actor(mapper)
    camera = setup_camera()
    light = setup_light()
    spot = setup_spotlight()
    renderer = setup_renderer(actor, camera, light, spot)
    render_window = setup_render_window(renderer)
    setup_interactor(render_window)
