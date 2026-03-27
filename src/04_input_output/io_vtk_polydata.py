import os
import sys
from pathlib import Path

import vtk

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.join(SCRIPT_DIR, "../../data/vtks/cube_polydata.vtk")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "../../data/vtks/cube_polydata_transformed.vtk")


def write_vtk_polydata(data: vtk.vtkPolyData, filename: str):
    """Write a vtkPolyData to a legacy VTK file."""
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_vtk_polydata(filename: str) -> vtk.vtkPolyData:
    """Read a legacy VTK PolyData file and return a vtkPolyData."""
    reader = vtk.vtkPolyDataReader()
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


def print_dataset_info(dataset: vtk.vtkPolyData):
    """Print information about the dataset's arrays."""
    point_data = dataset.GetPointData()
    cell_data = dataset.GetCellData()

    print(f"Number of points: {dataset.GetNumberOfPoints()}")
    print(f"Number of cells: {dataset.GetNumberOfCells()}")
    print(f"Number of polys: {dataset.GetNumberOfPolys()}")
    print(f"Number of lines: {dataset.GetNumberOfLines()}")

    print(f"Number of point data arrays: {point_data.GetNumberOfArrays()}")
    for j in range(point_data.GetNumberOfArrays()):
        array = point_data.GetArray(j)
        print(f"  Point data array {j}: {point_data.GetArrayName(j)}")
        print(f"    Number of components: {array.GetNumberOfComponents()}")
        print(f"    Number of tuples: {array.GetNumberOfTuples()}")
        print(f"    Range: {array.GetRange()}")

    print(f"Number of cell data arrays: {cell_data.GetNumberOfArrays()}")
    for j in range(cell_data.GetNumberOfArrays()):
        array = cell_data.GetArray(j)
        print(f"  Cell data array {j}: {cell_data.GetArrayName(j)}")
        print(f"    Number of components: {array.GetNumberOfComponents()}")
        print(f"    Number of tuples: {array.GetNumberOfTuples()}")
        print(f"    Range: {array.GetRange()}")


if __name__ == "__main__":

    # Read legacy VTK PolyData file
    polydata = read_vtk_polydata(FILE_NAME)

    # Print dataset info
    print_dataset_info(polydata)

    # Display for verification
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    polydata = transform_polydata(polydata, transform)

    # Write legacy VTK PolyData file
    write_vtk_polydata(polydata, OUTPUT_FILE)
