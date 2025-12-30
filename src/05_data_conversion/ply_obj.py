import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class PlyToObjConverter(Converter):
    """Converter for PLY to OBJ file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".ply")
            and output_filename.lower().endswith(".obj")
        ):
            raise ValueError("Invalid file extensions. Expected '.ply' and '.obj'.")

        try:
            # Read PLY File
            reader = vtk.vtkPLYReader()
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


class ObjToPlyConverter(Converter):
    """Converter for OBJ to PLY file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".obj")
            and output_filename.lower().endswith(".ply")
        ):
            raise ValueError("Invalid file extensions. Expected '.obj' and '.ply'.")

        try:
            # Read OBJ File
            reader = vtk.vtkOBJReader()
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

    # Example: OBJ to PLY conversion
    obj_input = os.path.join(data_dir, "objs", "cube.obj")
    ply_output = os.path.join(script_dir, "cube_converted.ply")

    obj_to_ply = ObjToPlyConverter()
    obj_to_ply.convert(obj_input, ply_output)
    print(f"Converted {obj_input} to {ply_output}")

    # Note: PLY to OBJ conversion requires a PLY input file
    # Uncomment the following if you have a PLY file available:
    # ply_input = os.path.join(data_dir, "plys", "example.ply")
    # obj_output = os.path.join(script_dir, "example_converted.obj")
    # ply_to_obj = PlyToObjConverter()
    # ply_to_obj.convert(ply_input, obj_output)
    # print(f"Converted {ply_input} to {obj_output}")
