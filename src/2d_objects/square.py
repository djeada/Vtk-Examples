import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline


if __name__ == "__main__":
    # create a square
    square = vtk.vtkCubeSource()
    square.SetXLength(1)
    square.SetYLength(1)
    square.SetZLength(0)

    # create a mapper
    squareMapper = vtk.vtkPolyDataMapper()
    squareMapper.SetInputConnection(square.GetOutputPort())

    # display the square
    pipeline = VisualisationPipeline([squareMapper])
    pipeline.run()
