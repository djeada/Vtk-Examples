import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class VtkToStlConverter(VTKConverter):
    """Converter for VTK to STL file format conversion."""

    input_extension = ".vtk"
    output_extension = ".stl"
    reader_type = vtk.vtkGenericDataObjectReader
    writer_type = vtk.vtkSTLWriter
    output_type = "vtkPolyData"
    extract_surface = True


class StlToVtkConverter(VTKConverter):
    """Converter for STL to VTK file format conversion."""

    input_extension = ".stl"
    output_extension = ".vtk"
    reader_type = vtk.vtkSTLReader
    writer_type = vtk.vtkPolyDataWriter
    output_type = "vtkDataSet"


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Get the directory of this script and construct paths to data files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "..", "data")

    # Example: VTK to STL conversion
    vtk_input = os.path.join(data_dir, "vtks", "cube_polydata.vtk")
    stl_output = os.path.join(script_dir, "cube_converted.stl")

    vtk_to_stl = VtkToStlConverter()
    vtk_to_stl.convert(vtk_input, stl_output)
    print(f"Converted {vtk_input} to {stl_output}")

    # Example: STL to VTK conversion
    stl_input = os.path.join(data_dir, "stls", "cube.stl")
    vtk_output = os.path.join(script_dir, "cube_converted.vtk")

    stl_to_vtk = StlToVtkConverter()
    stl_to_vtk.convert(stl_input, vtk_output)
    print(f"Converted {stl_input} to {vtk_output}")
