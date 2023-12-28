import vtk

from src.common.simple_pipeline import VisualisationPipeline

FOAM_CASE_PATH = "../../data/open_foam/temp.case"


def write_vtk(data: vtk.vtkMultiBlockDataSet, filename: str):
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_foam(case_path: str) -> vtk.vtkMultiBlockDataSet:
    reader = vtk.vtkOpenFOAMReader()
    reader.SetFileName(case_path + "/system/controlDict")
    reader.Update()
    return reader.GetOutput()


def transform_block(block: vtk.vtkMultiBlockDataSet, transform: vtk.vtkTransform):
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

    # Apply a translation transform
    transform = vtk.vtkTransform()
    transform.Scale(2.0, 1.0, 1.0)  # Double the size along the x-axis
    transform.RotateZ(45.0)  # Rotate 45 degrees about the z-axis
    foam_data = transform_block(foam_data, transform)

    # Write VTK file
    write_vtk(foam_data, "test.vtk")
