import vtk

from src.simple_pipeline import VisualisationPipeline

FILE_NAME = "../../data/objs/cube.obj"


def write_obj(data: vtk.vtkPolyData, filename: str):
    writer = vtk.vtkOBJWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_obj(filename: str) -> vtk.vtkPolyData:
    reader = vtk.vtkOBJReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def transform_polydata(polydata: vtk.vtkPolyData, transform: vtk.vtkTransform):
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

    # Apply a translation transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    obj_data = transform_polydata(obj_data, transform)

    # Write OBJ file
    write_obj(obj_data, "../../data/objs/cube.obj")
