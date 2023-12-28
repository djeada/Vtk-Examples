import vtk

from src.common.simple_pipeline import VisualisationPipeline

FILE_NAME = "../../data/vtms/grid_of_triangles.vtm"


def write_vtm(data: vtk.vtkMultiBlockDataSet, filename: str):
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(data)
    writer.Write()


def read_vtm(filename: str) -> vtk.vtkMultiBlockDataSet:
    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


def print_dataset_info(dataset: vtk.vtkMultiBlockDataSet):
    for i in range(dataset.GetNumberOfBlocks()):
        block = dataset.GetBlock(i)
        if block:
            point_data = block.GetPointData()
            cell_data = block.GetCellData()

            print(f"Block {i}:")
            print(f"  Number of point data arrays: {point_data.GetNumberOfArrays()}")
            for j in range(point_data.GetNumberOfArrays()):
                array = point_data.GetArray(j)
                print(f"  Point data array {j}: {point_data.GetArrayName(j)}")
                print(f"    Number of components: {array.GetNumberOfComponents()}")
                print(f"    Number of tuples: {array.GetNumberOfTuples()}")
                print(f"    Range: {array.GetRange()}")

            print(f"  Number of cell data arrays: {cell_data.GetNumberOfArrays()}")
            for j in range(cell_data.GetNumberOfArrays()):
                array = cell_data.GetArray(j)
                print(f"  Cell data array {j}: {cell_data.GetArrayName(j)}")
                print(f"    Number of components: {array.GetNumberOfComponents()}")
                print(f"    Number of tuples: {array.GetNumberOfTuples()}")
                print(f"    Range: {array.GetRange()}")


if __name__ == "__main__":

    # Read VTM file
    vtm_data = read_vtm(FILE_NAME)

    # Print dataset info
    print_dataset_info(vtm_data)

    # Use vtkCompositeDataGeometryFilter to convert the multiblock dataset to polydata
    geometryFilter = vtk.vtkCompositeDataGeometryFilter()
    geometryFilter.SetInputDataObject(vtm_data)
    geometryFilter.Update()

    # Create a mapper and actor for the polydata
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometryFilter.GetOutputPort())

    # Display
    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

    # Write VTM file
    write_vtm(vtm_data, "test.vtm")
