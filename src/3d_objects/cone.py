import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

if __name__ == "__main__":
    cone = vtk.vtkConeSource()
    cone.SetHeight(5.0)
    cone.SetRadius(2.0)
    cone.SetResolution(5)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cone.GetOutputPort())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
