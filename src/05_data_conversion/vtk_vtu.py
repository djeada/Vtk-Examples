import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class VtkToVtuConverter(VTKConverter):
    """Converter for VTK to VTU (XML unstructured grid) file format conversion."""

    input_extension = ".vtk"
    output_extension = ".vtu"
    reader_type = vtk.vtkGenericDataObjectReader
    writer_type = vtk.vtkXMLUnstructuredGridWriter
    output_type = "vtkUnstructuredGrid"


class VtuToVtkConverter(VTKConverter):
    """Converter for VTU (XML unstructured grid) to VTK file format conversion."""

    input_extension = ".vtu"
    output_extension = ".vtk"
    reader_type = vtk.vtkXMLUnstructuredGridReader
    writer_type = vtk.vtkDataSetWriter
    output_type = "vtkDataSet"


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
