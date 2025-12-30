import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class ObjToVtpConverter(Converter):
    """Converter for OBJ to VTP (XML PolyData) file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".obj")
            and output_filename.lower().endswith(".vtp")
        ):
            raise ValueError("Invalid file extensions. Expected '.obj' and '.vtp'.")

        try:
            # Read OBJ File
            reader = vtk.vtkOBJReader()
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


class VtpToObjConverter(Converter):
    """Converter for VTP (XML PolyData) to OBJ file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtp")
            and output_filename.lower().endswith(".obj")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtp' and '.obj'.")

        try:
            # Read VTP File
            reader = vtk.vtkXMLPolyDataReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write OBJ File
            writer = vtk.vtkOBJWriter()
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

    # Example: OBJ to VTP conversion
    obj_input = os.path.join(data_dir, "objs", "cube.obj")
    vtp_output = os.path.join(script_dir, "cube_converted.vtp")

    obj_to_vtp = ObjToVtpConverter()
    obj_to_vtp.convert(obj_input, vtp_output)
    print(f"Converted {obj_input} to {vtp_output}")

    # Example: VTP to OBJ conversion
    vtp_input = os.path.join(data_dir, "vtps", "naca0012.vtp")
    obj_output = os.path.join(script_dir, "naca0012_converted.obj")

    vtp_to_obj = VtpToObjConverter()
    vtp_to_obj.convert(vtp_input, obj_output)
    print(f"Converted {vtp_input} to {obj_output}")
