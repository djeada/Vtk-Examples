from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import vtkPolyData, vtkMultiBlockDataSet
from vtkmodules.vtkFiltersCore import vtkAppendFilter
from vtkmodules.vtkIOLegacy import vtkGenericDataObjectWriter

output_path = "output.vtk"


def generate_multi_block_dataset():

    # create a multi-block dataset
    multi_block_dataset = vtkMultiBlockDataSet()

    # create a block consisting of a collection of points in cube shape
    block = vtkPolyData()
    # create a points array
    points = vtkPoints()
    for i in range(100):
        points.InsertNextPoint(i, i, i)
    # set points array to block
    block.SetPoints(points)
    # set block to multi-block dataset
    multi_block_dataset.SetBlock(0, block)
    return multi_block_dataset


def append_blocks_recursive(multi_block_dataset):
    # get number of blocks
    n_blocks = multi_block_dataset.GetNumberOfBlocks()
    # iterate over blocks
    for i in range(n_blocks):
        # get block
        block = multi_block_dataset.GetBlock(i)
        if block is None:
            return
        # if block is multi-block dataset call this function recursively
        if block.IsA("vtkMultiBlockDataSet"):
            append_blocks_recursive(block)
        else:
            append_filter.AddInputData(block)


def write_vtk(file_name, data, writer_type=vtkGenericDataObjectWriter):
    writer = writer_type()
    writer.SetFileName(file_name)
    writer.SetInputData(data)
    writer.Write()


if __name__ == "__main__":
    multi_block_dataset = generate_multi_block_dataset()

    append_filter = vtkAppendFilter()
    append_blocks_recursive(multi_block_dataset)
    append_filter.Update()
    data = append_filter.GetOutput()

    write_vtk(output_path, data)
