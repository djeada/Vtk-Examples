import vtk
from converter_interface import Converter


class VTKtoVTMConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(input_filename)
        reader.Update()

        # Create a multiblock dataset and add the unstructured grid to it
        mbds = vtk.vtkMultiBlockDataSet()
        mbds.SetNumberOfBlocks(1)
        mbds.SetBlock(0, reader.GetOutput())

        writer = vtk.vtkXMLMultiBlockDataWriter()
        writer.SetFileName(output_filename)
        writer.SetInputData(mbds)
        writer.Write()


class VTMtoVTKConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkXMLMultiBlockDataReader()
        reader.SetFileName(input_filename)
        reader.Update()

        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(output_filename)
        writer.SetInputConnection(reader.GetOutputPort())
        writer.Write()
