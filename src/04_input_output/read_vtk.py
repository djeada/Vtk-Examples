from vtkmodules.vtkIOLegacy import vtkUnstructuredGridReader
from vtk import (
    vtkDataSetMapper,
    vtkActor,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
)

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "../../data/unstructured_grid.vtk"


def read_vtk(file_name, reader_type=vtkUnstructuredGridReader):
    reader = reader_type()
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    return output


if __name__ == "__main__":

    output = read_vtk(file_name)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
