import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class VtkToVtpConverter(Converter):
    """Converter for VTK to VTP (XML PolyData) file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtk")
            and output_filename.lower().endswith(".vtp")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.vtp'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
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


class VtpToVtkConverter(Converter):
    """Converter for VTP (XML PolyData) to VTK file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtp")
            and output_filename.lower().endswith(".vtk")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtp' and '.vtk'.")

        try:
            # Read VTP File
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTK File
            writer = vtk.vtkPolyDataWriter()
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

    # Example: VTK to VTP conversion
    vtk_input = os.path.join(data_dir, "vtks", "cube_polydata.vtk")
    vtp_output = os.path.join(script_dir, "cube_converted.vtp")

    vtk_to_vtp = VtkToVtpConverter()
    vtk_to_vtp.convert(vtk_input, vtp_output)
    print(f"Converted {vtk_input} to {vtp_output}")

    # Example: VTP to VTK conversion
    vtp_input = os.path.join(data_dir, "vtps", "naca0012.vtp")
    vtk_output = os.path.join(script_dir, "naca0012_converted.vtk")

    vtp_to_vtk = VtpToVtkConverter()
    vtp_to_vtk.convert(vtp_input, vtk_output)
    print(f"Converted {vtp_input} to {vtk_output}")
