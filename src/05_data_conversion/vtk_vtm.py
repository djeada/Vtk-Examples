import logging
import os

import vtk

try:
    from .converter_interface import VTKConverter
except ImportError:
    from converter_interface import VTKConverter


class VtkToVtmConverter(VTKConverter):
    """Wrap a legacy VTK dataset in a single-block VTM file."""

    input_extension = ".vtk"
    output_extension = ".vtm"
    reader_type = vtk.vtkGenericDataObjectReader
    writer_type = vtk.vtkXMLMultiBlockDataWriter
    output_type = "vtkMultiBlockDataSet"

    def prepare_data(self, data):
        if not data.IsA("vtkDataSet") or data.GetNumberOfPoints() == 0:
            raise ValueError("VTM conversion requires a nonempty VTK dataset.")
        multiblock = vtk.vtkMultiBlockDataSet()
        multiblock.SetNumberOfBlocks(1)
        multiblock.SetBlock(0, data)
        return multiblock


class VtmToVtkConverter(VTKConverter):
    """Export the first block of a VTM file as a legacy VTK dataset.

    Legacy dataset writers cannot preserve multiple blocks. Warn when other
    blocks are omitted and reject a missing or nested first block explicitly.
    """

    input_extension = ".vtm"
    output_extension = ".vtk"
    reader_type = vtk.vtkXMLMultiBlockDataReader
    writer_type = vtk.vtkDataSetWriter
    output_type = "vtkDataSet"

    def prepare_data(self, data):
        count = data.GetNumberOfBlocks()
        if count == 0:
            raise ValueError("No blocks found in the VTM file.")
        if count > 1:
            logging.warning(
                "VTM file contains %s blocks. Only the first block will be converted.",
                count,
            )
        block = data.GetBlock(0)
        if block is None:
            raise ValueError("The first VTM block is empty.")
        return block


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Get the directory of this script and construct paths to data files
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, "..", "..", "data")

    # Example: VTK to VTM conversion
    vtk_input = os.path.join(data_dir, "vtks", "cube_polydata.vtk")
    vtm_output = os.path.join(script_dir, "cube_converted.vtm")

    vtk_to_vtm = VtkToVtmConverter()
    vtk_to_vtm.convert(vtk_input, vtm_output)
    print(f"Converted {vtk_input} to {vtm_output}")

    # Example: VTM to VTK conversion
    vtm_input = os.path.join(data_dir, "vtms", "grid_of_triangles.vtm")
    vtk_output = os.path.join(script_dir, "grid_converted.vtk")

    vtm_to_vtk = VtmToVtkConverter()
    vtm_to_vtk.convert(vtm_input, vtk_output)
    print(f"Converted {vtm_input} to {vtk_output}")
