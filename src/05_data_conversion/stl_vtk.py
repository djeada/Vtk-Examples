import vtk
from vtkmodules.vtkIOLegacy import vtkPolyDataWriter, vtkPolyDataReader
from vtkmodules.vtkRenderingCore import vtkDataSetMapper

from src.io.read_stl import read_stl
from src.io.read_vtk import read_vtk
from src.io.write_vtk import write_vtk
from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

input_file_name = "../../data/cube.stl"
output_file_name = "example.vtk"


def convert_stl_to_vtk(input_file_name, output_file_name):
    print(f"Attempting to convert {input_file_name} to {output_file_name}")
    data = read_stl(input_file_name)
    write_vtk(output_file_name, data, vtkPolyDataWriter)
    print(f"Conversion successful")


if __name__ == "__main__":
    convert_stl_to_vtk(input_file_name, output_file_name)

    # check the conversion result
    output = read_vtk(output_file_name, vtkPolyDataReader)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
