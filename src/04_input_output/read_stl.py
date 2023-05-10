import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "../../data/cube.stl"


def read_stl(file_name):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()
    return output


if __name__ == "__main__":

    output = read_stl(file_name)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(output)

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
