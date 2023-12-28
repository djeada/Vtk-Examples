import vtk
from converter_interface import Converter


class VTKtoSTLConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(input_filename)
        reader.Update()

        writer = vtk.vtkSTLWriter()
        writer.SetFileName(output_filename)
        writer.SetInputConnection(reader.GetOutputPort())
        writer.Write()
