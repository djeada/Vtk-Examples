import math

import vtk

# Create a VTK renderer and render window
renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

# Create a VTK render window interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(render_window)

# Create a VTK cone source and mapper
cone_source = vtk.vtkConeSource()
cone_mapper = vtk.vtkPolyDataMapper()
cone_mapper.SetInputConnection(cone_source.GetOutputPort())

# Create a VTK actor and set the mapper
cone_actor = vtk.vtkActor()
cone_actor.SetMapper(cone_mapper)

# Add the actor to the renderer
renderer.AddActor(cone_actor)

# Set up the camera
camera = vtk.vtkCamera()
renderer.SetActiveCamera(camera)

# Set initial camera parameters
camera.SetPosition(0, 0, 5)
camera.SetFocalPoint(0, 0, 0)
camera.SetViewUp(0, 1, 0)

# Define camera animation parameters
num_frames = 100
delta_angle = 360.0 / num_frames


def animate_camera():
    for frame in range(num_frames):
        angle = frame * delta_angle
        rad_angle = math.radians(angle)  # Convert angle to radians

        # Update camera position
        camera.SetPosition(5 * math.cos(rad_angle), 0, 5 * math.sin(rad_angle))

        # Update the render window
        render_window.Render()


animate_camera()

# Start the interactor
interactor.Start()
