import numpy as np
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def create_airfoil_geometry():
    """
    Create a basic airfoil geometry using an elliptical cross-section.

    Returns:
        vtk.vtkPolyData: The airfoil geometry as a VTK PolyData object.
    """
    # Parameters for the ellipse (airfoil cross-section)
    radius_x, radius_y = 1.0, 0.2  # Major and minor radii of the ellipse
    num_points = 100  # Number of points to define the ellipse

    # Create points for an ellipse
    ellipse_points = vtk.vtkPoints()
    for i in range(num_points):
        angle = 2.0 * np.pi * i / num_points
        x = radius_x * np.cos(angle)
        y = radius_y * np.sin(angle)
        ellipse_points.InsertNextPoint(x, y, 0)

    # Create a polygon that represents the ellipse
    polygon = vtk.vtkPolygon()
    polygon.GetPointIds().SetNumberOfIds(num_points)
    for i in range(num_points):
        polygon.GetPointIds().SetId(i, i)

    # Create a PolyData object
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(ellipse_points)

    # Create a cell array to store the polygon
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polygon)
    polydata.SetPolys(cells)

    # Extrude the polygon to create a 3D airfoil shape
    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.SetInputData(polydata)
    extrude.SetExtrusionTypeToNormalExtrusion()
    extrude.SetVector(0.0, 0.0, 1.0)
    extrude.SetScaleFactor(5.0)  # Length of the airfoil

    # Convert the extruded shape to PolyData
    extrude.Update()
    airfoil_geometry = extrude.GetOutput()

    return airfoil_geometry


def generate_flow_field(airfoil_geometry):
    """
    Generate a vector field representing airflow around the airfoil.

    Args:
        airfoil_geometry (vtk.vtkPolyData): The airfoil geometry.

    Returns:
        vtk.vtkStructuredGrid: A structured grid representing the vector field.
    """
    # Define the bounds and resolution of the flow field
    bounds = airfoil_geometry.GetBounds()
    grid_resolution = 20  # Adjust for finer or coarser grid

    # Create a mesh grid for vector field
    x = np.linspace(bounds[0] - 1, bounds[1] + 1, grid_resolution)
    y = np.linspace(bounds[2] - 1, bounds[3] + 1, grid_resolution)
    z = np.linspace(bounds[4] - 1, bounds[5] + 1, grid_resolution)
    x, y, z = np.meshgrid(x, y, z, indexing="ij")

    # Initialize the vector field (flow primarily along x-axis)
    flow_vectors = np.zeros((grid_resolution, grid_resolution, grid_resolution, 3))
    flow_vectors[..., 0] = (
        1  # Setting x-component of vectors to simulate horizontal flow
    )

    # This is a simple approximation and can be improved for more realistic simulations
    for i in range(grid_resolution):
        for j in range(grid_resolution):
            for k in range(grid_resolution):
                # Check if the point is close to the airfoil and modify the flow vector accordingly
                if (
                    np.sqrt((x[i, j, k] - 0) ** 2 + (y[i, j, k] - 0) ** 2) < 1.0
                ):  # Example condition
                    flow_vectors[i, j, k, 1] += 0.5  # Modify y-component as an example

    # Convert numpy array to VTK structured grid
    vectors_flat = flow_vectors.ravel()
    vtk_vectors = numpy_to_vtk(vectors_flat, deep=True, array_type=vtk.VTK_FLOAT)
    vtk_vectors.SetNumberOfComponents(3)

    structured_grid = vtk.vtkStructuredGrid()
    structured_grid.SetDimensions(grid_resolution, grid_resolution, grid_resolution)

    points = vtk.vtkPoints()
    for i in range(grid_resolution):
        for j in range(grid_resolution):
            for k in range(grid_resolution):
                points.InsertNextPoint(x[i, j, k], y[i, j, k], z[i, j, k])

    structured_grid.SetPoints(points)
    structured_grid.GetPointData().SetVectors(vtk_vectors)

    return structured_grid


def calculate_pressure_distribution(airfoil_geometry, flow_field):
    """
    Approximate the pressure distribution on the airfoil based on the flow field.

    Args:
        airfoil_geometry (vtk.vtkPolyData): The airfoil geometry.
        flow_field (vtk.vtkStructuredGrid): The flow field around the airfoil.

    Returns:
        np.ndarray: An array representing the pressure distribution on the airfoil's surface.
    """
    # Assume air density (rho) in kg/m^3 (1.225 kg/m^3 at sea level, 15 degrees Celsius)
    rho = 1.225
    # Assume free stream speed (V_inf) in m/s
    V_inf = 1.0
    # Assume free stream pressure (P_inf) in Pascals
    P_inf = 101325  # Standard atmospheric pressure at sea level

    # Extract the flow vectors from the flow field
    vtk_vectors = flow_field.GetPointData().GetVectors()
    vectors = vtk_to_numpy(vtk_vectors)

    # Calculate the speed at each point in the flow field
    speeds = np.linalg.norm(vectors, axis=1)

    # Calculate pressure at each point using Bernoulli's equation
    # P = P_inf - 0.5 * rho * (V^2 - V_inf^2)
    pressures = P_inf - 0.5 * rho * (speeds**2 - V_inf**2)

    # Ensure no negative pressures (optional, depending on your application)
    pressures = np.maximum(pressures, 0)

    # Pressure distribution along the airfoil surface
    # For a more accurate pressure distribution, you would project the pressures onto the airfoil surface
    # This may involve finding the nearest points on the airfoil to the flow field points and interpolating the pressures

    # Here we provide a simplified mapping assuming the pressures array aligns with the airfoil geometry points
    # In a real CFD application, you would perform an interpolation based on the airfoil's surface projection

    return pressures


def visualize_simulation(airfoil_geometry, flow_field, pressure_distribution):
    """
    Use VTK to visualize the airfoil, the flow field, and the pressure distribution.

    Args:
        airfoil_geometry (vtk.vtkPolyData): The airfoil geometry.
        flow_field (vtk.vtkStructuredGrid): The flow field around the airfoil.
        pressure_distribution (np.ndarray): The pressure distribution on the airfoil's surface.
    """
    # Create a mapper and actor for the airfoil
    airfoil_mapper = vtk.vtkPolyDataMapper()
    airfoil_mapper.SetInputData(airfoil_geometry)

    airfoil_actor = vtk.vtkActor()
    airfoil_actor.SetMapper(airfoil_mapper)

    # Apply a color map to the airfoil based on the pressure distribution
    pressure_colors = vtk.vtkFloatArray()
    pressure_colors.SetName("Pressure")
    for pressure in pressure_distribution:
        pressure_colors.InsertNextValue(pressure)
    airfoil_geometry.GetPointData().SetScalars(pressure_colors)
    airfoil_mapper.SetScalarRange(
        np.min(pressure_distribution), np.max(pressure_distribution)
    )

    # Increase the number of seed points and adjust their placement
    seeds = vtk.vtkPointSource()
    seeds.SetNumberOfPoints(100)  # Increase the number of seed points
    seeds.SetRadius(5)  # Adjust the radius as needed
    seeds.SetCenter(0, 0, 0)  # Place the seed points in an area with significant flow

    # Streamline tracer setup with adjusted parameters
    streamer = vtk.vtkStreamTracer()
    streamer.SetInputData(flow_field)
    streamer.SetSourceConnection(seeds.GetOutputPort())
    streamer.SetMaximumPropagation(200)  # Increase the propagation length
    streamer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
    streamer.SetInitialIntegrationStep(0.1)
    streamer.SetIntegrationDirectionToBoth()
    streamer.SetComputeVorticity(False)

    # Streamline visualization
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(streamer.GetOutputPort())

    streamline_actor = vtk.vtkActor()
    streamline_actor.SetMapper(mapper)
    streamline_actor.GetProperty().SetColor(1, 0, 0)  # Set a visible color

    # Renderer and Render Window setup
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    # Add actors to the renderer
    renderer.AddActor(airfoil_actor)
    renderer.AddActor(streamline_actor)  # Add streamline actor

    # Set background color and reset camera
    renderer.SetBackground(0.1, 0.2, 0.3)
    renderer.ResetCamera()

    # Start visualization
    render_window.Render()
    interactor.Start()


def main():
    airfoil_geometry = create_airfoil_geometry()
    flow_field = generate_flow_field(airfoil_geometry)
    pressure_distribution = calculate_pressure_distribution(
        airfoil_geometry, flow_field
    )
    visualize_simulation(airfoil_geometry, flow_field, pressure_distribution)


if __name__ == "__main__":
    main()
