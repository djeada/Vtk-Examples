import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class PlyToStlConverter(VTKConverter):
    """Converter for PLY to STL file format conversion."""

    input_extension = ".ply"
    output_extension = ".stl"
    reader_type = vtk.vtkPLYReader
    writer_type = vtk.vtkSTLWriter
    output_type = "vtkPolyData"


class StlToPlyConverter(VTKConverter):
    """Converter for STL to PLY file format conversion."""

    input_extension = ".stl"
    output_extension = ".ply"
    reader_type = vtk.vtkSTLReader
    writer_type = vtk.vtkPLYWriter
    output_type = "vtkPolyData"


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
