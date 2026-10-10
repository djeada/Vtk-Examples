"""Demonstrate a pandas fallback for a Spark-to-VTK data workflow.

The sample creates and filters point data with pandas, then converts it to VTK.
It does not require a Spark or Java installation.
"""

import os
import sys
import numpy as np
import pandas as pd
import vtk
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def create_sample_data_without_spark(num_points=1000):
    """Generate sample 3D points with a scalar value using only Pandas."""
    np.random.seed(42)
    data = {
        "x": np.random.normal(0, 1, num_points),
        "y": np.random.normal(0, 1, num_points),
        "z": np.random.normal(0, 1, num_points),
    }
    df = pd.DataFrame(data)
    df["distance"] = np.sqrt(df["x"] ** 2 + df["y"] ** 2 + df["z"] ** 2)
    return df


def process_data_without_spark(df):
    """Process data using Pandas instead of Spark."""
    return df[df["distance"] < 2.0].copy()


def convert_to_vtk(df):
    """Convert DataFrame to VTK PolyData."""
    try:
        # Create points array
        points = vtk.vtkPoints()
        vertices = vtk.vtkCellArray()

        # Create array for scalar values
        scalars = vtk.vtkDoubleArray()
        scalars.SetName("Distance")

        # Add points and scalar values
        for idx, row in df.iterrows():
            point_id = points.InsertNextPoint(row["x"], row["y"], row["z"])
            scalars.InsertNextValue(row["distance"])

            # Create a vertex cell for each point
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, point_id)
            vertices.InsertNextCell(vertex)

        # Create PolyData
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts(vertices)
        polydata.GetPointData().SetScalars(scalars)

        return polydata
    except Exception as e:
        logger.error(f"Error in converting to VTK: {str(e)}")
        raise


def visualize_vtk(polydata):
    """Create a basic VTK visualization."""
    try:
        # Create mapper and actor
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polydata)
        mapper.ScalarVisibilityOn()
        mapper.SetScalarRange(polydata.GetPointData().GetScalars().GetRange())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # Create renderer and window
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(0.1, 0.1, 0.1)

        window = vtk.vtkRenderWindow()
        window.AddRenderer(renderer)
        window.SetSize(800, 600)

        # Add interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(window)

        # Initialize and start
        interactor.Initialize()
        window.Render()
        interactor.Start()
    except Exception as e:
        logger.error(f"Error in visualization: {str(e)}")
        raise


def main():
    try:
        logger.info("Starting pandas fallback for the Spark-to-VTK workflow")

        # Create and process data using Pandas
        logger.info("Generating sample point data with pandas")
        df = create_sample_data_without_spark()

        logger.info("Processing data")
        processed_df = process_data_without_spark(df)

        # Save to parquet (optional)
        parquet_path = "processed_data.parquet"
        logger.info(f"Saving data to {parquet_path}")
        processed_df.to_parquet(parquet_path)

        # Convert to VTK
        logger.info("Converting data to VTK format")
        vtk_data = convert_to_vtk(processed_df)

        # Visualize
        logger.info("Starting visualization")
        visualize_vtk(vtk_data)

        # Clean up
        if os.path.exists(parquet_path):
            os.remove(parquet_path)
            logger.info("Cleaned up temporary files")

    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
