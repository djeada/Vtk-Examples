import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class ObjToVtpConverter(VTKConverter):
    """Converter for OBJ to VTP (XML PolyData) file format conversion."""

    input_extension = ".obj"
    output_extension = ".vtp"
    reader_type = vtk.vtkOBJReader
    writer_type = vtk.vtkXMLPolyDataWriter
    output_type = "vtkPolyData"


class VtpToObjConverter(VTKConverter):
    """Converter for VTP (XML PolyData) to OBJ file format conversion."""

    input_extension = ".vtp"
    output_extension = ".obj"
    reader_type = vtk.vtkXMLPolyDataReader
    writer_type = vtk.vtkOBJWriter
    output_type = "vtkPolyData"


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
