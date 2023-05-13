import vtk

from converter_interface import Converter


class VTKtoVTUConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(input_filename)
        reader.Update()

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output_filename)
        writer.SetInputConnection(reader.GetOutputPort())
        writer.Write()


class VTUtoVTKConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(input_filename)
        reader.Update()

        writer = vtk.vtkDataSetWriter()
        writer.SetFileName(output_filename)
        writer.SetInputConnection(reader.GetOutputPort())
        writer.Write()
