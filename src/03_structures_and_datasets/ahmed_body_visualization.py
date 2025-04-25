import vtk
import numpy as np


# Load the STL file into VTK
def load_stl_geometry(file_path):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(file_path)
    reader.Update()
    return reader.GetOutput()


# Apply scalar values to the geometry
def apply_scalar_values(geometry):
    scalar_values = vtk.vtkDoubleArray()
    scalar_values.SetName("StaticPressure")
    scalar_values.SetNumberOfComponents(1)

    num_points = geometry.GetNumberOfPoints()
    for i in range(num_points):
        scalar_value = np.sin(i / num_points * np.pi * 2) * 300  # Example values
        scalar_values.InsertNextValue(scalar_value)

    geometry.GetPointData().SetScalars(scalar_values)


# Create a color legend bar
def create_color_legend(mapper):
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(mapper.GetLookupTable())
    scalar_bar.SetTitle("Static Pressure")
    scalar_bar.GetLabelTextProperty().SetColor(0, 0, 0)
    scalar_bar.GetTitleTextProperty().SetColor(0, 0, 0)
    scalar_bar.SetNumberOfLabels(5)
    return scalar_bar


# Visualize the geometry with scalar values
def visualize_geometry(geometry):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(geometry)
    mapper.SetScalarRange(-300, 300)
    mapper.SetScalarModeToUsePointData()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1)

    # Create and add the color legend bar
    scalar_bar = create_color_legend(mapper)
    renderer.AddActor2D(scalar_bar)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)

    render_window_interactor = vtk.vtkRenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)

    # Set up the title and labels
    title = vtk.vtkTextActor()
    title.SetInput("Ahmed Body CFD Simulation")
    title.GetTextProperty().SetFontSize(24)
    title.GetTextProperty().SetColor(0, 0, 0)
    title.SetPosition2(10, 570)
    renderer.AddActor2D(title)

    x_label = vtk.vtkTextActor()
    x_label.SetInput("X-Axis")
    x_label.GetTextProperty().SetFontSize(14)
    x_label.GetTextProperty().SetColor(0, 0, 0)
    x_label.SetPosition2(700, 30)
    renderer.AddActor2D(x_label)

    y_label = vtk.vtkTextActor()
    y_label.SetInput("Y-Axis")
    y_label.GetTextProperty().SetFontSize(14)
    y_label.GetTextProperty().SetColor(0, 0, 0)
    y_label.SetPosition2(30, 550)
    renderer.AddActor2D(y_label)

    render_window.Render()
    render_window_interactor.Start()


# Main function to download, load, and visualize the Ahmed body
def main():
    save_path = "../../data/stls/ahmed_body.stl"

    # Load the STL file
    car_geometry = load_stl_geometry(save_path)

    # Apply scalar values to the geometry
    apply_scalar_values(car_geometry)

    # Visualize the geometry with scalar values
    visualize_geometry(car_geometry)


if __name__ == "__main__":
    main()
