import vtk
import time


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
    return actor


def setup_camera():
    camera = vtk.vtkCamera()
    camera.SetFocalPoint(0.0, 0.0, 0.0)
    camera.SetPosition(0.0, 0.0, 8.0)
    camera.SetViewUp(0.0, 1.0, 0.0)
    return camera


def setup_renderer(actor, camera):
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetActiveCamera(camera)
    return renderer


def setup_render_window(renderer):
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    return render_window


if __name__ == "__main__":
    cone = setup_data()
    mapper = setup_mapper(cone)
    actor = setup_actor(mapper)
    camera = setup_camera()
    renderer = setup_renderer(actor, camera)
    render_window = setup_render_window(renderer)

    azimuth = 0
    while 1:
        render_window.Render()
        if azimuth >= 360:
            azimuth = 0
        azimuth += 0.1
        camera.Azimuth(azimuth)

        time.sleep(0.1)
