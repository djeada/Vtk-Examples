import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class VtkToObjConverter(Converter):
    """Converter for VTK to OBJ file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtk")
            and output_filename.lower().endswith(".obj")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.obj'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
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


class ObjToVtkConverter(Converter):
    """Converter for OBJ to VTK file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".obj")
            and output_filename.lower().endswith(".vtk")
        ):
            raise ValueError("Invalid file extensions. Expected '.obj' and '.vtk'.")

        try:
            # Read OBJ File
            reader = vtk.vtkOBJReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTK File
            writer = vtk.vtkDataSetWriter()
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

    # Example: VTK to OBJ conversion
    vtk_input = os.path.join(data_dir, "vtks", "cube_polydata.vtk")
    obj_output = os.path.join(script_dir, "cube_converted.obj")

    vtk_to_obj = VtkToObjConverter()
    vtk_to_obj.convert(vtk_input, obj_output)
    print(f"Converted {vtk_input} to {obj_output}")

    # Example: OBJ to VTK conversion
    obj_input = os.path.join(data_dir, "objs", "cube.obj")
    vtk_output = os.path.join(script_dir, "cube_converted.vtk")

    obj_to_vtk = ObjToVtkConverter()
    obj_to_vtk.convert(obj_input, vtk_output)
    print(f"Converted {obj_input} to {vtk_output}")
