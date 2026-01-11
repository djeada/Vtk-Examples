import vtk
import math
from datetime import datetime, UTC


class DistanceToPointFilter(vtk.vtkPolyDataAlgorithm):
    """
    A VTK filter that computes the Euclidean distance from each point in the input
    vtkPolyData to a specified target point.
    """

    def __init__(self):
        super().__init__()
        self.__TargetPoint = [0.0, 0.0, 0.0]
        self.GetExecutive().SetOutputData(0, vtk.vtkPolyData())

    def SetTargetPoint(self, x, y, z):
        try:
            self.__TargetPoint = [float(x), float(y), float(z)]
            self.Modified()
            return True
        except ValueError as e:
            print(f"Error setting target point: {e}")
            return False

    def GetTargetPoint(self):
        return tuple(self.__TargetPoint)

    def ProcessDataObject(self, input_data):
        if not input_data:
            print("No input data provided")
            return None

        output_data = vtk.vtkPolyData()
        output_data.DeepCopy(input_data)

        num_points = input_data.GetNumberOfPoints()
        print(f"Processing {num_points} points...")

        distances = vtk.vtkDoubleArray()
        distances.SetName("DistanceToTarget")
        distances.SetNumberOfComponents(1)
        distances.SetNumberOfTuples(num_points)

        tx, ty, tz = self.__TargetPoint
        min_dist = float("inf")
        max_dist = float("-inf")

        for i in range(num_points):
            point = input_data.GetPoint(i)
            dx = point[0] - tx
            dy = point[1] - ty
            dz = point[2] - tz
            distance = math.sqrt(dx * dx + dy * dy + dz * dz)
            distances.SetValue(i, distance)
            min_dist = min(min_dist, distance)
            max_dist = max(max_dist, distance)

        print(f"Distance range: [{min_dist:.3f}, {max_dist:.3f}]")
        output_data.GetPointData().SetScalars(distances)
        return output_data


def create_actor(poly_data, scalar_range=None, color=None):
    """Create a VTK actor from polydata with optional scalar coloring or solid color."""
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(poly_data)

    if scalar_range and poly_data.GetPointData().GetScalars():
        mapper.SetScalarRange(scalar_range)

        # Create lookup table for distance visualization
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)  # blue to red
        lut.SetSaturationRange(1.0, 1.0)
        lut.SetValueRange(1.0, 1.0)
        lut.SetTableRange(scalar_range)
        lut.Build()

        mapper.SetLookupTable(lut)
        mapper.SetScalarVisibility(1)
    else:
        mapper.SetScalarVisibility(0)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    if color:
        actor.GetProperty().SetColor(color)

    return actor


def create_visualization():
    try:
        # Create source
        print("Creating sphere source...")
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetThetaResolution(30)
        sphere.SetPhiResolution(30)
        sphere.Update()

        input_data = sphere.GetOutput()
        print(f"Sphere created with {input_data.GetNumberOfPoints()} points")

        # Process with filter
        print("Creating and executing distance filter...")
        filter = DistanceToPointFilter()
        output_data = filter.ProcessDataObject(input_data)

        if not output_data:
            raise RuntimeError("Filter failed to process data")

        print(f"Filter output has {output_data.GetNumberOfPoints()} points")

        # Get scalar range for coloring
        scalars = output_data.GetPointData().GetScalars()
        if not scalars:
            raise RuntimeError("No scalar data in filter output")

        scalar_range = scalars.GetRange()
        print(f"Scalar range: [{scalar_range[0]:.3f}, {scalar_range[1]:.3f}]")

        # Create actors for both spheres
        original_actor = create_actor(input_data, color=(0.8, 0.8, 1.0))  # Light blue
        filtered_actor = create_actor(output_data, scalar_range=scalar_range)

        # Create renderers for split view
        left_renderer = vtk.vtkRenderer()
        right_renderer = vtk.vtkRenderer()

        # Set viewport (x1, y1, x2, y2) - split screen in half
        left_renderer.SetViewport(0.0, 0.0, 0.5, 1.0)
        right_renderer.SetViewport(0.5, 0.0, 1.0, 1.0)

        # Add actors to respective renderers
        left_renderer.AddActor(original_actor)
        right_renderer.AddActor(filtered_actor)

        # Set background colors
        left_renderer.SetBackground(0.1, 0.1, 0.1)
        right_renderer.SetBackground(0.2, 0.2, 0.2)

        # Create render window
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(left_renderer)
        render_window.AddRenderer(right_renderer)
        render_window.SetSize(1200, 600)  # Wider window for side-by-side view

        # Create interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)

        # Add orientation widgets to both views
        for renderer, viewport_coords in [
            (left_renderer, (0.0, 0.0, 0.2, 0.2)),
            (right_renderer, (0.5, 0.0, 0.7, 0.2)),
        ]:
            axes = vtk.vtkAxesActor()
            widget = vtk.vtkOrientationMarkerWidget()
            widget.SetOrientationMarker(axes)
            widget.SetInteractor(interactor)
            widget.SetViewport(*viewport_coords)
            widget.SetEnabled(1)
            widget.InteractiveOn()

        # Add titles
        left_text = vtk.vtkTextActor()
        left_text.SetInput("Original Sphere")
        left_text.GetTextProperty().SetFontSize(24)
        left_text.SetPosition(10, 10)
        left_renderer.AddActor2D(left_text)

        right_text = vtk.vtkTextActor()
        right_text.SetInput("Distance-Colored Sphere")
        right_text.GetTextProperty().SetFontSize(24)
        right_text.SetPosition(10, 10)
        right_renderer.AddActor2D(right_text)

        # Initialize view
        left_renderer.ResetCamera()
        right_renderer.ResetCamera()

        # Link cameras
        right_renderer.SetActiveCamera(left_renderer.GetActiveCamera())

        render_window.Render()
        interactor.Initialize()

        print("Visualization setup complete")
        return interactor

    except Exception as e:
        print(f"Error in visualization setup: {str(e)}")
        import traceback

        print(f"Traceback:\n{traceback.format_exc()}")
        return None


if __name__ == "__main__":
    print(
        "Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted):",
        datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S"),
    )
    print("Current User's Login: djeada")

    interactor = create_visualization()
    if interactor:
        print("Starting interaction...")
        interactor.Start()
    else:
        print("Failed to setup visualization.")
