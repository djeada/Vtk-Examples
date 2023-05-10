import numpy
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support
from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
from vtkmodules.vtkFiltersCore import vtkAppendFilter
from vtkmodules.vtkIOExodus import vtkExodusIIReader
from vtkmodules.vtkIOLegacy import (
    vtkUnstructuredGridReader,
    vtkStructuredPointsWriter,
    vtkStructuredGridWriter,
    vtkGenericDataObjectWriter,
)
from vtkmodules.vtkRenderingCore import vtkPolyDataMapper, vtkDataSetMapper

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "../../data/disk_out_ref.ex2"


def read_exodus_ii(file_name):
    reader = vtkExodusIIReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    return output


def display_block_info(block):
    print("Number of points: ", block.GetNumberOfPoints())
    print("Number of cells: ", block.GetNumberOfCells())
    print("Number of point data arrays: ", block.GetPointData().GetNumberOfArrays())
    print("Number of cell data arrays: ", block.GetCellData().GetNumberOfArrays())
    print("Number of field data arrays: ", block.GetFieldData().GetNumberOfArrays())
    print("")

    # get point data
    point_data = block.GetPointData()
    # get number of point data arrays
    n_arrays = point_data.GetNumberOfArrays()
    # iterate over point data arrays
    for j in range(n_arrays):
        # get point data array
        array = point_data.GetArray(j)
        # get point data array name
        array_name = array.GetName()
        # get point data array type
        array_type = array.GetDataTypeAsString()
        # get point data array size
        array_size = array.GetNumberOfComponents()
        # get point data array range
        array_range = array.GetRange()
        # get point data array values
        array_values = numpy_support.vtk_to_numpy(array)
        # print point data array info
        print("Point data array name: ", array_name)
        print("Point data array type: ", array_type)
        print("Point data array size: ", array_size)
        print("Point data array range: ", array_range)
        print("Point data array values: ", array_values)
        print("")

    # get cell data
    cell_data = block.GetCellData()
    # get number of cell data arrays
    n_arrays = cell_data.GetNumberOfArrays()
    # iterate over cell data arrays
    for j in range(n_arrays):
        # get cell data array
        array = cell_data.GetArray(j)
        # get cell data array name
        array_name = array.GetName()
        # get cell data array type
        array_type = array.GetDataTypeAsString()
        # get cell data array size
        array_size = array.GetNumberOfComponents()
        # get cell data array range
        array_range = array.GetRange()
        # get cell data array values
        array_values = numpy_support.vtk_to_numpy(array)
        # print cell data array info
        print("Cell data array name: ", array_name)
        print("Cell data array type: ", array_type)
        print("Cell data array size: ", array_size)
        print("Cell data array range: ", array_range)
        print("Cell data array values: ", array_values)
        print("")


def investigate_multi_block_dataset(multi_block_dataset):
    # get number of blocks
    n_blocks = multi_block_dataset.GetNumberOfBlocks()
    # iterate over blocks
    for i in range(n_blocks):
        # get block
        block = multi_block_dataset.GetBlock(i)
        if block is None:
            return
        # get block name
        block_name = multi_block_dataset.GetMetaData(i).Get(vtkCompositeDataSet.NAME())
        # print block info
        print("Block name: ", block_name)
        # if block is multi-block dataset call this function recursively
        if block.IsA("vtkMultiBlockDataSet"):
            investigate_multi_block_dataset(block)
        else:
            display_block_info(block)


def setup_mappers(multi_block_dataset):
    # get number of blocks
    n_blocks = multi_block_dataset.GetNumberOfBlocks()
    # iterate over blocks
    for i in range(n_blocks):
        # get block
        block = multi_block_dataset.GetBlock(i)
        if block is None:
            return
        # get block name
        block_name = multi_block_dataset.GetMetaData(i).Get(vtkCompositeDataSet.NAME())
        # if block is multi-block dataset call this function recursively
        if block.IsA("vtkMultiBlockDataSet"):
            setup_mappers(block)
        else:
            mapper = vtkDataSetMapper()
            mapper.SetInputData(block)
            mapper.SetScalarRange(block.GetScalarRange())
            mappers.append(mapper)


if __name__ == "__main__":

    output = read_exodus_ii(file_name)
    investigate_multi_block_dataset(output)

    # now we have to setup mappers recursively for multi-block dataset
    mappers = []
    setup_mappers(output)

    pipeline = VisualisationPipeline(mappers)
    pipeline.run()
