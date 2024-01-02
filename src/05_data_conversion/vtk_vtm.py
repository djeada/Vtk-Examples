import vtk
import logging
from converter_interface import Converter


class VtkToVtmConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.vtk') and output_filename.lower().endswith('.vtm')):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.vtm'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Create a multiblock dataset and add the output to it
            mbds = vtk.vtkMultiBlockDataSet()
            mbds.SetNumberOfBlocks(1)
            mbds.SetBlock(0, reader.GetOutput())

            # Write VTM File
            writer = vtk.vtkXMLMultiBlockDataWriter()
            writer.SetFileName(output_filename)
            writer.SetInputData(mbds)
            writer.Write()
            logging.info(f"Conversion successful: {input_filename} to {output_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


class VtmToVtkConverter(Converter):
    def convert(self, input_filename: str, output_filename: str):
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (input_filename.lower().endswith('.vtm') and output_filename.lower().endswith('.vtk')):
            raise ValueError("Invalid file extensions. Expected '.vtm' and '.vtk'.")

        try:
            # Read VTM File
            reader = vtk.vtkXMLMultiBlockDataReader()
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

