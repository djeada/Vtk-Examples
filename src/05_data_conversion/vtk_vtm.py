import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class VtkToVtmConverter(Converter):
    """Converter for VTK to VTM (multiblock) file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtk")
            and output_filename.lower().endswith(".vtm")
        ):
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
            logging.info(
                f"Conversion successful: {input_filename} to {output_filename}"
            )
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


class VtmToVtkConverter(Converter):
    """Converter for VTM (multiblock) to VTK file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtm")
            and output_filename.lower().endswith(".vtk")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtm' and '.vtk'.")

        try:
            # Read VTM File
            reader = vtk.vtkXMLMultiBlockDataReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Get the first block from the multiblock dataset
            mbds = reader.GetOutput()
            num_blocks = mbds.GetNumberOfBlocks()
            if num_blocks > 0:
                if num_blocks > 1:
                    logging.warning(
                        f"VTM file contains {num_blocks} blocks. "
                        "Only the first block will be converted."
                    )
                block = mbds.GetBlock(0)
                # Write VTK File
                writer = vtk.vtkDataSetWriter()
                writer.SetFileName(output_filename)
                writer.SetInputData(block)
                writer.Write()
                logging.info(
                    f"Conversion successful: {input_filename} to {output_filename}"
                )
            else:
                raise ValueError(f"No blocks found in the VTM file: {input_filename}")
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Get the directory of this script and construct paths to data files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "..", "data")

    # Example: VTK to VTM conversion
    vtk_input = os.path.join(data_dir, "vtks", "cube_polydata.vtk")
    vtm_output = os.path.join(script_dir, "cube_converted.vtm")

    vtk_to_vtm = VtkToVtmConverter()
    vtk_to_vtm.convert(vtk_input, vtm_output)
    print(f"Converted {vtk_input} to {vtm_output}")

    # Example: VTM to VTK conversion
    vtm_input = os.path.join(data_dir, "vtms", "grid_of_triangles.vtm")
    vtk_output = os.path.join(script_dir, "grid_converted.vtk")

    vtm_to_vtk = VtmToVtkConverter()
    vtm_to_vtk.convert(vtm_input, vtk_output)
    print(f"Converted {vtm_input} to {vtk_output}")
