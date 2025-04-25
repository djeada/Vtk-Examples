import vtk
import numpy as np
from vtkmodules.vtkInteractionWidgets import vtkSliderWidget, vtkSliderRepresentation2D
from vtkmodules.vtkFiltersSources import vtkCylinderSource
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkPolyDataMapper,
)
from vtkmodules.vtkCommonColor import vtkNamedColors

# Create the cylinder
cylinder = vtkCylinderSource()
cylinder.SetHeight(10)
cylinder.SetRadius(2)
cylinder.SetResolution(100)
cylinder.Update()

# Get the cylinder's points
points = cylinder.GetOutput().GetPoints()

# Store the original point coordinates for reference
original_points = np.array(
    [points.GetPoint(i) for i in range(points.GetNumberOfPoints())]
)

# Mapper
cylinderMapper = vtkPolyDataMapper()
cylinderMapper.SetInputConnection(cylinder.GetOutputPort())

# Actor
cylinderActor = vtkActor()
cylinderActor.SetMapper(cylinderMapper)

# Renderer
renderer = vtkRenderer()
renderer.AddActor(cylinderActor)
renderer.SetBackground(0, 0, 0)  # Background color black

# Render Window
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.SetSize(1200, 800)  # Increased size for more space

# Interactor
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Global variables to store slider values
current_bending = 0
current_compression = 0
current_torsion = 0
current_shear = 0


def apply_deformation():
    for i in range(points.GetNumberOfPoints()):
        x, y, z = original_points[i]

        # Apply compression along the Z-axis
        z = z * (1 - current_compression)

        # Apply bending (displacement along Y-axis based on Z, similar to image)
        if current_bending != 0:
            curvature = current_bending / 100.0  # Adjust curvature sensitivity
            y = y + curvature * z**2

        # Apply torsion (twisting, rotation around Z-axis based on Z, like a corkscrew)
        if current_torsion != 0:
            theta = current_torsion * (z / 10.0)
            x_new = x * np.cos(np.radians(theta)) - y * np.sin(np.radians(theta))
            y_new = x * np.sin(np.radians(theta)) + y * np.cos(np.radians(theta))
            x, y = x_new, y_new

        # Apply shear (lateral displacement along XZ-plane)
        x = x + current_shear * z

        points.SetPoint(i, x, y, z)

    # Update the cylinder with new point data
    points.Modified()
    renderWindow.Render()


# Slider callback functions
def update_bending(obj, event):
    global current_bending
    current_bending = obj.GetRepresentation().GetValue()
    apply_deformation()


def update_compression(obj, event):
    global current_compression
    current_compression = obj.GetRepresentation().GetValue()
    apply_deformation()


def update_torsion(obj, event):
    global current_torsion
    current_torsion = obj.GetRepresentation().GetValue()
    apply_deformation()


def update_shear(obj, event):
    global current_shear
    current_shear = obj.GetRepresentation().GetValue()
    apply_deformation()


# Creating Sliders
def create_slider(min_val, max_val, title, y_position):
    sliderRep = vtkSliderRepresentation2D()
    sliderRep.SetMinimumValue(min_val)
    sliderRep.SetMaximumValue(max_val)
    sliderRep.SetValue(0.0)
    sliderRep.SetTitleText(title)
    sliderRep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    sliderRep.GetPoint1Coordinate().SetValue(0.01, y_position)  # Farther to the left
    sliderRep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    sliderRep.GetPoint2Coordinate().SetValue(
        0.3, y_position
    )  # Enough distance from cylinder

    sliderWidget = vtkSliderWidget()
    sliderWidget.SetInteractor(renderWindowInteractor)
    sliderWidget.SetRepresentation(sliderRep)
    sliderWidget.EnabledOn()

    return sliderWidget


# Adjusted positions for sliders to avoid overlapping
bending_slider = create_slider(-30, 30, "Bending", 0.8)
bending_slider.AddObserver("EndInteractionEvent", update_bending)

compression_slider = create_slider(0, 0.5, "Compression", 0.6)
compression_slider.AddObserver("EndInteractionEvent", update_compression)

torsion_slider = create_slider(-180, 180, "Torsion", 0.4)
torsion_slider.AddObserver("EndInteractionEvent", update_torsion)

shear_slider = create_slider(-1, 1, "Shear", 0.2)
shear_slider.AddObserver("EndInteractionEvent", update_shear)

# Start interaction
renderWindow.Render()
renderWindowInteractor.Start()
