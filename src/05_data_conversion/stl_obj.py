import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class StlToObjConverter(VTKConverter):
    """Converter for STL to OBJ file format conversion."""

    input_extension = ".stl"
    output_extension = ".obj"
    reader_type = vtk.vtkSTLReader
    writer_type = vtk.vtkOBJWriter
    output_type = "vtkPolyData"


class ObjToStlConverter(VTKConverter):
    """Converter for OBJ to STL file format conversion."""

    input_extension = ".obj"
    output_extension = ".stl"
    reader_type = vtk.vtkOBJReader
    writer_type = vtk.vtkSTLWriter
    output_type = "vtkPolyData"


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
