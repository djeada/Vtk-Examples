import os
import sys
import numpy as np
import pandas as pd
import vtk
import math
from datetime import datetime
import logging
from typing import Tuple, List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def generate_sample_geospatial_data(n_points: int = 10000) -> pd.DataFrame:
    """
    Generate synthetic geospatial data with the following attributes:
    - latitude: -90 to 90
    - longitude: -180 to 180
    - population_density: simulated based on latitude (higher near equator)
    - elevation: simulated with some mountain ranges
    - pollution_index: simulated based on population density
    """
    np.random.seed(42)

    # Generate random coordinates with more points near populated areas
    lat = np.random.normal(loc=0, scale=30, size=n_points)  # More points near equator
    lat = np.clip(lat, -90, 90)
    lon = np.random.uniform(-180, 180, n_points)

    # Simulate population density (higher near equator and certain longitudes)
    population_density = 1000 * np.exp(-np.abs(lat) / 30) * (1 + 0.5 * np.sin(np.radians(lon)))
    population_density *= np.random.lognormal(0, 0.5, n_points)

    # Simulate elevation (create some mountain ranges)
    elevation = np.zeros(n_points)
    for mountain_lat, mountain_lon in [(35, 140), (45, -120), (30, 80)]:  # Example mountain ranges
        distance = np.sqrt(
            (lat - mountain_lat) ** 2 +
            (lon - mountain_lon) ** 2
        )
        elevation += 8000 * np.exp(-distance ** 2 / 1000)
    elevation += np.random.normal(0, 200, n_points)  # Add noise
    elevation = np.maximum(0, elevation)  # No negative elevation

    # Simulate pollution (correlated with population density but with random variation)
    pollution_index = 0.7 * population_density / population_density.max()
    pollution_index += 0.3 * np.random.beta(2, 5, n_points)
    pollution_index = np.clip(pollution_index, 0, 1)

    return pd.DataFrame({
        'latitude': lat,
        'longitude': lon,
        'population_density': population_density,
        'elevation': elevation,
        'pollution_index': pollution_index
    })


def lat_lon_to_xyz(lat: float, lon: float, elevation: float = 0, radius: float = 6371000) -> Tuple[float, float, float]:
    """
    Convert latitude and longitude to 3D coordinates.
    Args:
        lat: Latitude in degrees
        lon: Longitude in degrees
        elevation: Elevation above sea level in meters
        radius: Earth's radius in meters (default: 6,371,000 meters)
    Returns:
        Tuple of (x, y, z) coordinates
    """
    # Convert to radians
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Calculate radius with elevation
    r = radius + elevation

    # Convert to Cartesian coordinates
    x = r * math.cos(lat_rad) * math.cos(lon_rad)
    y = r * math.cos(lat_rad) * math.sin(lon_rad)
    z = r * math.sin(lat_rad)

    return (x, y, z)


def create_earth_surface() -> vtk.vtkPolyData:
    """Create a basic Earth surface sphere."""
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(6371000)  # Earth's radius in meters
    sphere.SetThetaResolution(100)
    sphere.SetPhiResolution(100)
    sphere.Update()
    return sphere.GetOutput()


def create_vtk_points(df: pd.DataFrame) -> Tuple[vtk.vtkPoints, List[vtk.vtkDoubleArray]]:
    """
    Convert DataFrame to VTK points and arrays.
    Returns points and arrays for population density, elevation, and pollution.
    """
    points = vtk.vtkPoints()

    # Create arrays for each attribute
    population_array = vtk.vtkDoubleArray()
    population_array.SetName("Population Density")

    elevation_array = vtk.vtkDoubleArray()
    elevation_array.SetName("Elevation")

    pollution_array = vtk.vtkDoubleArray()
    pollution_array.SetName("Pollution Index")

    # Convert each point and add data
    for idx, row in df.iterrows():
        x, y, z = lat_lon_to_xyz(row['latitude'], row['longitude'], row['elevation'])
        points.InsertNextPoint(x, y, z)

        population_array.InsertNextValue(row['population_density'])
        elevation_array.InsertNextValue(row['elevation'])
        pollution_array.InsertNextValue(row['pollution_index'])

    return points, [population_array, elevation_array, pollution_array]


def create_visualization(earth_surface: vtk.vtkPolyData,
                         points: vtk.vtkPoints,
                         data_arrays: List[vtk.vtkDoubleArray]) -> None:
    """Create the visualization with Earth surface and data points."""
    # Create polydata for points
    point_data = vtk.vtkPolyData()
    point_data.SetPoints(points)

    # Add arrays to point data
    for array in data_arrays:
        point_data.GetPointData().AddArray(array)

    # Create vertex cells for each point
    vertices = vtk.vtkCellArray()
    for i in range(points.GetNumberOfPoints()):
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, i)
        vertices.InsertNextCell(vertex)
    point_data.SetVerts(vertices)

    # Create mappers and actors
    # Earth surface
    earth_mapper = vtk.vtkPolyDataMapper()
    earth_mapper.SetInputData(earth_surface)

    earth_actor = vtk.vtkActor()
    earth_actor.SetMapper(earth_mapper)
    earth_actor.GetProperty().SetColor(0.3, 0.3, 0.3)  # Dark gray
    earth_actor.GetProperty().SetOpacity(0.5)

    # Data points
    point_mapper = vtk.vtkPolyDataMapper()
    point_mapper.SetInputData(point_data)
    point_mapper.ScalarVisibilityOn()
    point_mapper.SetScalarModeToUsePointData()
    point_mapper.SelectColorArray("Population Density")
    point_mapper.SetScalarRange(
        point_data.GetPointData().GetArray("Population Density").GetRange()
    )

    point_actor = vtk.vtkActor()
    point_actor.SetMapper(point_mapper)
    point_actor.GetProperty().SetPointSize(3)

    # Create renderer and window
    renderer = vtk.vtkRenderer()
    renderer.AddActor(earth_actor)
    renderer.AddActor(point_actor)
    renderer.SetBackground(0.1, 0.1, 0.1)

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(1024, 768)

    # Create interactor and add key bindings
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

    # Add custom interactor style
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)

    # Add key bindings for switching between data arrays
    def create_key_press_callback(arrays: List[vtk.vtkDoubleArray]):
        def key_press_callback(obj, event):
            key = obj.GetKeySym()
            if key == "1":
                point_mapper.SelectColorArray("Population Density")
                point_mapper.SetScalarRange(arrays[0].GetRange())
            elif key == "2":
                point_mapper.SelectColorArray("Elevation")
                point_mapper.SetScalarRange(arrays[1].GetRange())
            elif key == "3":
                point_mapper.SelectColorArray("Pollution Index")
                point_mapper.SetScalarRange(arrays[2].GetRange())
            window.Render()

        return key_press_callback

    interactor.AddObserver("KeyPressEvent", create_key_press_callback(data_arrays))

    # Initialize and start
    interactor.Initialize()
    window.Render()

    # Add on-screen instructions
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(
        "Controls:\n"
        "Left mouse: Rotate\n"
        "Middle mouse: Pan\n"
        "Right mouse: Zoom\n"
        "1: Population Density\n"
        "2: Elevation\n"
        "3: Pollution Index"
    )
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)
    text_actor.SetPosition(10, 10)
    renderer.AddActor2D(text_actor)

    interactor.Start()


def main():
    try:
        # Generate sample data
        logger.info("Generating sample geospatial data...")
        df = generate_sample_geospatial_data()

        # Save to parquet (simulating Spark output)
        parquet_path = "geospatial_data.parquet"
        logger.info(f"Saving data to {parquet_path}")
        df.to_parquet(parquet_path)

        # Read from parquet (simulating loading Spark results)
        logger.info("Reading data from parquet...")
        df = pd.read_parquet(parquet_path)

        # Create VTK objects
        logger.info("Creating VTK visualization...")
        earth_surface = create_earth_surface()
        points, data_arrays = create_vtk_points(df)

        # Visualize
        logger.info("Starting visualization...")
        create_visualization(earth_surface, points, data_arrays)

        # Cleanup
        if os.path.exists(parquet_path):
            os.remove(parquet_path)
            logger.info("Cleaned up temporary files")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()