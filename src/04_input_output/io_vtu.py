import os

import vtk

from src.common.simple_pipeline import VisualisationPipeline

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_NAME = os.path.join(SCRIPT_DIR, "../../data/vtus/grid_of_triangles.vtu")
OUTPUT_FILE = os.path.join(
    SCRIPT_DIR, "../../data/vtus/grid_of_triangles_transformed.vtu"
)


def write_vtu(data: vtk.vtkUnstructuredGrid, filename: str):
    """Write a vtkUnstructuredGrid to a VTU file."""
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_vtu(filename: str) -> vtk.vtkUnstructuredGrid:
    """Read a VTU file and return a vtkUnstructuredGrid."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def transform_grid(grid: vtk.vtkUnstructuredGrid, transform: vtk.vtkTransform):
    """Apply a transform to a vtkUnstructuredGrid."""
    transform_filter = vtk.vtkTransformFilter()
    transform_filter.SetInputData(grid)
    transform_filter.SetTransform(transform)
    transform_filter.Update()
    return transform_filter.GetOutput()


def print_dataset_info(dataset: vtk.vtkUnstructuredGrid):
    """Print information about the dataset's arrays."""
    point_data = dataset.GetPointData()
    cell_data = dataset.GetCellData()

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

    # Read VTU file
    vtu_data = read_vtu(FILE_NAME)

    # Print dataset info
    print_dataset_info(vtu_data)

    # Display for verification
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(vtu_data)

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Apply a transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    vtu_data = transform_grid(vtu_data, transform)

    # Write VTU file
    write_vtu(vtu_data, OUTPUT_FILE)
