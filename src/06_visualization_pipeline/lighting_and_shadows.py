import vtk


def create_scene():
    # Create a renderer and render window
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)

    # Create an interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Create a light source
    light = vtk.vtkLight()
    light.SetFocalPoint(0, 0, 0)
    light.SetPosition(0, 0, 10)
    light.SetColor(1.0, 1.0, 1.0)
    renderer.AddLight(light)

    # Create a sphere
    sphere_source = vtk.vtkSphereSource()
    sphere_source.SetRadius(1.0)
    sphere_mapper = vtk.vtkPolyDataMapper()
    sphere_mapper.SetInputConnection(sphere_source.GetOutputPort())
    sphere_actor = vtk.vtkActor()
    sphere_actor.SetMapper(sphere_mapper)
    sphere_actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # Set sphere color to red
    sphere_actor.SetPosition(-1.5, 0, 0)  # Set sphere position
    renderer.AddActor(sphere_actor)

    # Create a cube
    cube_source = vtk.vtkCubeSource()
    cube_source.SetXLength(1.0)
    cube_source.SetYLength(1.0)
    cube_source.SetZLength(1.0)
    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube_source.GetOutputPort())
    cube_actor = vtk.vtkActor()
    cube_actor.SetMapper(cube_mapper)
    cube_actor.GetProperty().SetColor(0.0, 1.0, 0.0)  # Set cube color to green
    cube_actor.SetPosition(0, 0, 0)  # Set cube position
    renderer.AddActor(cube_actor)

    # Create a cone
    cone_source = vtk.vtkConeSource()
    cone_source.SetHeight(2.0)
    cone_source.SetRadius(1.0)
    cone_mapper = vtk.vtkPolyDataMapper()
    cone_mapper.SetInputConnection(cone_source.GetOutputPort())
    cone_actor = vtk.vtkActor()
    cone_actor.SetMapper(cone_mapper)
    cone_actor.GetProperty().SetColor(0.0, 0.0, 1.0)  # Set cone color to blue
    cone_actor.SetPosition(1.5, 0, 0)  # Set cone position
    renderer.AddActor(cone_actor)

    # Set camera position and orientation
    camera = vtk.vtkCamera()
    camera.SetPosition(0, -5, 5)
    camera.SetFocalPoint(0, 0, 0)
    renderer.SetActiveCamera(camera)

    # Set background color to black
    renderer.SetBackground(0.0, 0.0, 0.0)

    # Initialize the interactor and start the rendering loop
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


if __name__ == "__main__":
    create_scene()
