import vtk

from converter_interface import Converter


class STLtoOBJConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkSTLReader()
        reader.SetFileName(input_filename)
        reader.Update()

        writer = vtk.vtkOBJWriter()
        writer.SetFileName(output_filename)
        writer.SetInputConnection(reader.GetOutputPort())
        writer.Write()
