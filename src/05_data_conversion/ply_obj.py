import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class PlyToObjConverter(VTKConverter):
    """Converter for PLY to OBJ file format conversion."""

    input_extension = ".ply"
    output_extension = ".obj"
    reader_type = vtk.vtkPLYReader
    writer_type = vtk.vtkOBJWriter
    output_type = "vtkPolyData"


class ObjToPlyConverter(VTKConverter):
    """Converter for OBJ to PLY file format conversion."""

    input_extension = ".obj"
    output_extension = ".ply"
    reader_type = vtk.vtkOBJReader
    writer_type = vtk.vtkPLYWriter
    output_type = "vtkPolyData"


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
