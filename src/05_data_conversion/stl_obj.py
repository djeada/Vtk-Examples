import vtk
import logging
from converter_interface import Converter


class StlToObjConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.stl') and output_filename.lower().endswith('.obj')):
            raise ValueError("Invalid file extensions. Expected '.stl' and '.obj'.")

        try:
            # Read STL File
            reader = vtk.vtkSTLReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write OBJ File
            writer = vtk.vtkOBJWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise

class ObjToStlConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.obj') and output_filename.lower().endswith('.stl')):
            raise ValueError("Invalid file extensions. Expected '.obj' and '.stl'.")

        try:
            # Read OBJ File
            reader = vtk.vtkOBJReader()
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

