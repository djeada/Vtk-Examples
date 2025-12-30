import os

import vtk

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.join(SCRIPT_DIR, "../../data/objs/cube.obj")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../data/objs/cube_transformed.obj")


def write_obj(data: vtk.vtkPolyData, filename: str):
    """Write a vtkPolyData to an OBJ file."""
    writer = vtk.vtkOBJWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_obj(filename: str) -> vtk.vtkPolyData:
    """Read an OBJ file and return a vtkPolyData."""
    reader = vtk.vtkOBJReader()
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

    # Read OBJ file
    obj_data = read_obj(FILE_NAME)

    # Display for verification
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(obj_data)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    obj_data = transform_polydata(obj_data, transform)

    # Write OBJ file
    write_obj(obj_data, OUTPUT_FILE)
