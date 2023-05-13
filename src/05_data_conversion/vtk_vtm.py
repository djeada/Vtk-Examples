from vtkmodules.vtkRenderingCore import vtkDataSetMapper

from src.io.read_vtk import read_vtk
from src.io.read_vtm import read_vtm
from src.io.write_vtm import write_vtm
from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

input_file_name = "../../data/unstructured_grid.vtk"
output_file_name = "example.vtm"


def convert_vtk_to_vtm(input_file_name, output_file_name):
    print(f"Attempting to convert {input_file_name} to {output_file_name}")
    data = read_vtk(input_file_name)
    write_vtm(output_file_name, data)
    print(f"Conversion successful")


if __name__ == "__main__":
    convert_vtk_to_vtm(input_file_name, output_file_name)

    # check the conversion result
    output = read_vtm(output_file_name)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
