## Animations and Time-Varying Data

### Overview
* VTK provides tools for creating animations and visualizing time-varying data
* Some popular techniques include:
  - Keyframe Animation
  - Temporal Data Visualization
  - Animation Export

### Keyframe Animation
* Keyframe animation: creating smooth transitions between different states of a visualization
* Useful for visual storytelling or creating explanatory visualizations
* Example classes:
  - vtkAnimationCue: represents an animation sequence with start and end times
  - vtkAnimationScene: manages multiple animation cues and orchestrates the animation playback

### Temporal Data Visualization
* Temporal data visualization: displaying data that varies over time
* Useful for visualizing simulations, time series data, or dynamic systems
* Example classes:
  - vtkTemporalDataSet: represents a collection of datasets associated with different time steps
  - vtkTemporalInterpolator: interpolates between datasets at different time steps to create smooth transitions

### Animation Export
* Animation export: saving an animation as a video or image sequence
* Allows for sharing or offline viewing of animations
* Example classes:
  - vtkWindowToImageFilter: captures the contents of a render window as an image
  - vtkAVIWriter, vtkOggTheoraWriter, vtkFFMPEGWriter: save image sequences as video files

## Example: Keyframe Animation
```python
import vtk

# Create a sphere source
sphere = vtk.vtkSphereSource()
sphere.SetRadius(1)

# Create a mapper and actor for the sphere
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(sphere.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# Set up the renderer and render window
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

# Create an animation scene
animation_scene = vtk.vtkAnimationScene()
animation_scene.SetModeToSequence()
animation_scene.SetLoop(0)
animation_scene.SetFrameRate(24)

# Create an animation cue for the sphere radius
radius_cue = vtk.vtkAnimationCue()
radius_cue.SetStartTime(0)
radius_cue.SetEndTime(2)
animation_scene.AddCue(radius_cue)

# Define the callback function for updating the sphere radius
def update_sphere_radius(cue, event):
    t = cue.GetAnimationTime()
    sphere.SetRadius(1 + 0.5 * (1 + vtk.vtkMath.Sin(2 * vtk.vtkMath.Pi() * t)))

# Set up the callback and start the animation
radius_cue.AddObserver("AnimationCueTickEvent", update_sphere_radius)
animation_scene.Play()

# Clean up the animation scene
animation_scene.RemoveAllCues()
```
