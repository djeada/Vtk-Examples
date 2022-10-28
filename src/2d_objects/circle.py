import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline


if __name__ == "__main__":
    # create a circle
    circle = vtk.vtkRegularPolygonSource()
    circle.SetNumberOfSides(100)
    circle.SetRadius(1)
    circle.SetCenter(0, 0, 0)
    circle.SetNormal(0, 0, 1)

    # create a mapper
    circleMapper = vtk.vtkPolyDataMapper()
    circleMapper.SetInputConnection(circle.GetOutputPort())

    # display the circle
    pipeline = VisualisationPipeline([circleMapper])
    pipeline.run()
