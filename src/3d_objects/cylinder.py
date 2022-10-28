import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

if __name__ == "__main__":

    cylinder = vtk.vtkCylinderSource()
    cylinder.SetResolution(8)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cylinder.GetOutputPort())
    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
