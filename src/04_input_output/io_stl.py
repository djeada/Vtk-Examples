import os

import vtk

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.join(SCRIPT_DIR, "../../data/stls/cube.stl")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../data/stls/cube_transformed.stl")


def write_stl(data: vtk.vtkPolyData, filename: str):
    """Write a vtkPolyData to an STL file."""
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_stl(filename: str) -> vtk.vtkPolyData:
    """Read an STL file and return a vtkPolyData."""
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def transform_polydata(polydata: vtk.vtkPolyData, transform: vtk.vtkTransform):
    """Apply a transform to a vtkPolyData."""
    transform_filter = vtk.vtkTransformPolyDataFilter()
    transform_filter.SetInputData(polydata)
    transform_filter.SetTransform(transform)
    transform_filter.Update()
    return transform_filter.GetOutput()


if __name__ == "__main__":

    # Read STL file
    stl_data = read_stl(FILE_NAME)

    # Display for verification
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(stl_data)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    stl_data = transform_polydata(stl_data, transform)

    # Write STL file
    write_stl(stl_data, OUTPUT_FILE)
