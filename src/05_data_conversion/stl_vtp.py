import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class StlToVtpConverter(VTKConverter):
    """Converter for STL to VTP (XML PolyData) file format conversion."""

    input_extension = ".stl"
    output_extension = ".vtp"
    reader_type = vtk.vtkSTLReader
    writer_type = vtk.vtkXMLPolyDataWriter
    output_type = "vtkPolyData"


class VtpToStlConverter(VTKConverter):
    """Converter for VTP (XML PolyData) to STL file format conversion."""

    input_extension = ".vtp"
    output_extension = ".stl"
    reader_type = vtk.vtkXMLPolyDataReader
    writer_type = vtk.vtkSTLWriter
    output_type = "vtkPolyData"


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
