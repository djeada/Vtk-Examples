import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class StlToVtpConverter(Converter):
    """Converter for STL to VTP (XML PolyData) file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".stl")
            and output_filename.lower().endswith(".vtp")
        ):
            raise ValueError("Invalid file extensions. Expected '.stl' and '.vtp'.")

        try:
            # Read STL File
            reader = vtk.vtkSTLReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTP File
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(
                f"Conversion successful: {input_filename} to {output_filename}"
            )
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


class VtpToStlConverter(Converter):
    """Converter for VTP (XML PolyData) to STL file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtp")
            and output_filename.lower().endswith(".stl")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtp' and '.stl'.")

        try:
            # Read VTP File
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write STL File
            writer = vtk.vtkSTLWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(
                f"Conversion successful: {input_filename} to {output_filename}"
            )
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Get the directory of this script and construct paths to data files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "..", "data")

    # Example: STL to VTP conversion
    stl_input = os.path.join(data_dir, "stls", "cube.stl")
    vtp_output = os.path.join(script_dir, "cube_converted.vtp")

    stl_to_vtp = StlToVtpConverter()
    stl_to_vtp.convert(stl_input, vtp_output)
    print(f"Converted {stl_input} to {vtp_output}")

    # Example: VTP to STL conversion
    vtp_input = os.path.join(data_dir, "vtps", "naca0012.vtp")
    stl_output = os.path.join(script_dir, "naca0012_converted.stl")

    vtp_to_stl = VtpToStlConverter()
    vtp_to_stl.convert(vtp_input, stl_output)
    print(f"Converted {vtp_input} to {stl_output}")
