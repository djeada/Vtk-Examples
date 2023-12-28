import vtk
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet

from src.common.simple_pipeline import VisualisationPipeline

FILE_NAME = "../../data/exodus/tetrahedron.exo"


def read_exodus(filename: str) -> vtkMultiBlockDataSet:
    """Read an Exodus II file."""
    reader = vtk.vtkExodusIIReader()
    reader.SetFileName(filename)
    reader.UpdateInformation()
    reader.SetAllArrayStatus(
        vtk.vtkExodusIIReader.NODAL, 1
    )  # Enables the reading of all nodal fields.
    reader.Update()
    return reader.GetOutput()


def write_exodus(data: vtkMultiBlockDataSet, filename: str) -> None:
    """Write an Exodus II file."""
    writer = vtk.vtkExodusIIWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def translate_mesh(input_block, dx, dy, dz):
    if isinstance(input_block, vtk.vtkMultiBlockDataSet):
        for i in range(input_block.GetNumberOfBlocks()):
            sub_block = input_block.GetBlock(i)
            if sub_block is not None:
                translate_mesh(sub_block, dx, dy, dz)

    elif isinstance(input_block, vtk.vtkPointSet):
        points = input_block.GetPoints()
        for i in range(points.GetNumberOfPoints()):
            x, y, z = points.GetPoint(i)
            points.SetPoint(i, x + dx, y + dy, z + dz)


if __name__ == "__main__":
    # Read Exodus II file
    exodus_data = read_exodus(FILE_NAME)

    # Use vtkCompositeDataGeometryFilter to convert the multiblock dataset to polydata
    geometryFilter = vtk.vtkCompositeDataGeometryFilter()
    geometryFilter.SetInputDataObject(exodus_data)
    geometryFilter.Update()

    # Create a mapper and actor for the polydata
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometryFilter.GetOutputPort())
    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Translate the mesh along the x-axis by 10 units
    translate_mesh(exodus_data, 10.0, 0.0, 0.0)

    # Write Exodus II file
    write_exodus(exodus_data, "output.exo")
