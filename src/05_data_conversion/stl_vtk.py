import vtk
import logging
from converter_interface import Converter


class VtkToStlConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.vtk') and output_filename.lower().endswith('.stl')):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.stl'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write STL File
            writer = vtk.vtkSTLWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise

class StlToVtkConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.stl') and output_filename.lower().endswith('.vtk')):
            raise ValueError("Invalid file extensions. Expected '.stl' and '.vtk'.")

        try:
            # Read STL File
            reader = vtk.vtkSTLReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTK File
            writer = vtk.vtkPolyDataWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise
