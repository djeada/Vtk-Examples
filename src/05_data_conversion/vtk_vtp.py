import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class VtkToVtpConverter(VTKConverter):
    """Converter for VTK to VTP (XML PolyData) file format conversion."""

    input_extension = ".vtk"
    output_extension = ".vtp"
    reader_type = vtk.vtkGenericDataObjectReader
    writer_type = vtk.vtkXMLPolyDataWriter
    output_type = "vtkPolyData"
    extract_surface = True


class VtpToVtkConverter(VTKConverter):
    """Converter for VTP (XML PolyData) to VTK file format conversion."""

    input_extension = ".vtp"
    output_extension = ".vtk"
    reader_type = vtk.vtkXMLPolyDataReader
    writer_type = vtk.vtkPolyDataWriter
    output_type = "vtkDataSet"


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
