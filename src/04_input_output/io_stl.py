import vtk

from src.simple_pipeline import VisualisationPipeline

FILE_NAME = "../../data/stls/cube.stl"


def write_stl(data: vtk.vtkUnstructuredGrid, filename: str):
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_stl(filename: str) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def transform_grid(grid: vtk.vtkUnstructuredGrid, transform: vtk.vtkTransform):
    transform_filter = vtk.vtkTransformFilter()
    transform_filter.SetInputData(grid)
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

    # Apply a translation transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    stl_data = transform_grid(stl_data, transform)

    # Write STL file
    write_stl(stl_data, "test.stl")
