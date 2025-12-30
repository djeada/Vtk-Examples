import logging
import os

import vtk

try:
    from .converter_interface import Converter
except ImportError:
    from converter_interface import Converter


class VtkToVtuConverter(Converter):
    """Converter for VTK to VTU (XML unstructured grid) file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtk")
            and output_filename.lower().endswith(".vtu")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtk' and '.vtu'.")

        try:
            # Read VTK File
            reader = vtk.vtkGenericDataObjectReader()
            reader.SetFileName(input_filename)
            reader.Update()

            # Write VTU File
            writer = vtk.vtkXMLUnstructuredGridWriter()
            writer.SetFileName(output_filename)
            writer.SetInputConnection(reader.GetOutputPort())
            writer.Write()
            logging.info(
                f"Conversion successful: {input_filename} to {output_filename}"
            )
        except Exception as e:
            logging.error(f"Error during conversion: {e}")
            raise


class VtuToVtkConverter(Converter):
    """Converter for VTU (XML unstructured grid) to VTK file format conversion."""

    def convert(self, input_filename: str, output_filename: str) -> None:
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")

        if not (
            input_filename.lower().endswith(".vtu")
            and output_filename.lower().endswith(".vtk")
        ):
            raise ValueError("Invalid file extensions. Expected '.vtu' and '.vtk'.")

        try:
            # Read VTU File
            reader = vtk.vtkXMLUnstructuredGridReader()
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

    # Example: VTK to VTU conversion
    vtk_input = os.path.join(data_dir, "vtks", "grid_of_triangles.vtk")
    vtu_output = os.path.join(script_dir, "grid_converted.vtu")

    vtk_to_vtu = VtkToVtuConverter()
    vtk_to_vtu.convert(vtk_input, vtu_output)
    print(f"Converted {vtk_input} to {vtu_output}")

    # Example: VTU to VTK conversion
    vtu_input = os.path.join(data_dir, "vtus", "grid_of_triangles.vtu")
    vtk_output = os.path.join(script_dir, "grid_converted.vtk")

    vtu_to_vtk = VtuToVtkConverter()
    vtu_to_vtk.convert(vtu_input, vtk_output)
    print(f"Converted {vtu_input} to {vtk_output}")
