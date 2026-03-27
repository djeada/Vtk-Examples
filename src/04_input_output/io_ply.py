import os
import sys
from pathlib import Path

import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.join(SCRIPT_DIR, "../../data/plys/sphere.ply")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../data/plys/sphere_transformed.ply")


def write_ply(data: vtk.vtkPolyData, filename: str):
    """Write a vtkPolyData to a PLY file."""
    writer = vtk.vtkPLYWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_ply(filename: str) -> vtk.vtkPolyData:
    """Read a PLY file and return a vtkPolyData."""
    reader = vtk.vtkPLYReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def create_sample_ply(filename: str):
    """Create a sample PLY file from a sphere source for demonstration."""
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1.0)
    sphere.SetPhiResolution(20)
    sphere.SetThetaResolution(20)
    sphere.Update()

    os.makedirs(os.path.dirname(filename), exist_ok=True)
    write_ply(sphere.GetOutput(), filename)
    print(f"Created sample PLY file: {filename}")


def transform_polydata(polydata: vtk.vtkPolyData, transform: vtk.vtkTransform):
    """Apply a transform to a vtkPolyData."""
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(polydata)
    transform_filter.SetTransform(transform)
    transform_filter.Update()
    return transform_filter.GetOutput()


if __name__ == "__main__":

    # Create sample data if it doesn't exist
    if not os.path.exists(FILE_NAME):
        create_sample_ply(FILE_NAME)

    # Read PLY file
    ply_data = read_ply(FILE_NAME)

    print(f"Number of points: {ply_data.GetNumberOfPoints()}")
    print(f"Number of cells: {ply_data.GetNumberOfCells()}")

    # Display for verification
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(ply_data)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    ply_data = transform_polydata(ply_data, transform)

    # Write PLY file
    write_ply(ply_data, OUTPUT_FILE)
