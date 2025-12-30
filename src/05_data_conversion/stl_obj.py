import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class StlToObjConverter(Converter):
    """Converter for STL to OBJ file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".stl")
            and output_filename.lower().endswith(".obj")
        ):
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
            logging.info(
                f"Conversion successful: {input_filename} to {output_filename}"
            )
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


class ObjToStlConverter(Converter):
    """Converter for OBJ to STL file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".obj")
            and output_filename.lower().endswith(".stl")
        ):
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

    # Example: STL to OBJ conversion
    stl_input = os.path.join(data_dir, "stls", "cube.stl")
    obj_output = os.path.join(script_dir, "cube_converted.obj")

    stl_to_obj = StlToObjConverter()
    stl_to_obj.convert(stl_input, obj_output)
    print(f"Converted {stl_input} to {obj_output}")

    # Example: OBJ to STL conversion
    obj_input = os.path.join(data_dir, "objs", "cube.obj")
    stl_output = os.path.join(script_dir, "cube_converted.stl")

    obj_to_stl = ObjToStlConverter()
    obj_to_stl.convert(obj_input, stl_output)
    print(f"Converted {obj_input} to {stl_output}")
