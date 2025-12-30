import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class PlyToStlConverter(Converter):
    """Converter for PLY to STL file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".ply")
            and output_filename.lower().endswith(".stl")
        ):
            raise ValueError("Invalid file extensions. Expected '.ply' and '.stl'.")

        try:
            # Read PLY File
            reader = vtk.vtkPLYReader()
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


class StlToPlyConverter(Converter):
    """Converter for STL to PLY file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".stl")
            and output_filename.lower().endswith(".ply")
        ):
            raise ValueError("Invalid file extensions. Expected '.stl' and '.ply'.")

        try:
            # Read STL File
            reader = vtk.vtkSTLReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write PLY File
            writer = vtk.vtkPLYWriter()
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

    # Example: STL to PLY conversion
    stl_input = os.path.join(data_dir, "stls", "cube.stl")
    ply_output = os.path.join(script_dir, "cube_converted.ply")

    stl_to_ply = StlToPlyConverter()
    stl_to_ply.convert(stl_input, ply_output)
    print(f"Converted {stl_input} to {ply_output}")

    # Note: PLY to STL conversion requires a PLY input file
    # Uncomment the following if you have a PLY file available:
    # ply_input = os.path.join(data_dir, "plys", "example.ply")
    # stl_output = os.path.join(script_dir, "example_converted.stl")
    # ply_to_stl = PlyToStlConverter()
    # ply_to_stl.convert(ply_input, stl_output)
    # print(f"Converted {ply_input} to {stl_output}")
