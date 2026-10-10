import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class VtkToObjConverter(VTKConverter):
    """Converter for VTK to OBJ file format conversion."""

    input_extension = ".vtk"
    output_extension = ".obj"
    reader_type = vtk.vtkGenericDataObjectReader
    writer_type = vtk.vtkOBJWriter
    output_type = "vtkPolyData"
    extract_surface = True


class ObjToVtkConverter(VTKConverter):
    """Converter for OBJ to VTK file format conversion."""

    input_extension = ".obj"
    output_extension = ".vtk"
    reader_type = vtk.vtkOBJReader
    writer_type = vtk.vtkDataSetWriter
    output_type = "vtkDataSet"


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
