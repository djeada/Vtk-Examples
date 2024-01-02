import vtk
import logging
from converter_interface import Converter


class VtkToVtuConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.vtk') and output_filename.lower().endswith('.vtu')):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.vtu'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTU File
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise
class VtuToVtkConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.vtu') and output_filename.lower().endswith('.vtk')):
            raise ValueError("Invalid file extensions. Expected '.vtu' and '.vtk'.")

        try:
            # Read VTU File
            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTK File
            writer = vtk.vtkDataSetWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise
