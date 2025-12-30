import os

import vtk

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FOAM_CASE_PATH = os.path.join(SCRIPT_DIR, "../../data/open_foam/temp.case")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../data/open_foam/output.vtm")


def write_vtk(data: vtk.vtkMultiBlockDataSet, filename: str):
    """Write a vtkMultiBlockDataSet to a VTM file."""
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_foam(case_path: str) -> vtk.vtkMultiBlockDataSet:
    """
    Read an OpenFOAM case file.

    The case_path should point to a .foam file or a case directory.
    For a .foam file, it should be located in the case directory.
    """
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(case_path)
    reader.Update()
    return reader.GetOutput()


def transform_block(block: vtk.vtkMultiBlockDataSet, transform: vtk.vtkTransform):
    """Apply a transform to a vtkMultiBlockDataSet."""
    transform_filter = vtk.vtkTransformFilter()
    transform_filter.SetInputData(block)
    transform_filter.SetTransform(transform)
    transform_filter.Update()
    return transform_filter.GetOutput()


if __name__ == "__main__":

    # Read OpenFOAM case
    foam_data = read_foam(FOAM_CASE_PATH)

    # Display for verification
    mapper = vtk.vtkCompositePolyDataMapper2()
    mapper.SetInputData(foam_data)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    foam_data = transform_block(foam_data, transform)

    # Write VTK file
    write_vtk(foam_data, OUTPUT_FILE)
