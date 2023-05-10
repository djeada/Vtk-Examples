import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline


def setup_cone(height, radius, resolution):
    cone = vtk.vtkConeSource()
    cone.SetHeight(height)
    cone.SetRadius(radius)
    cone.SetResolution(resolution)
    return cone


if __name__ == "__main__":

    # define sources
    cone1 = setup_cone(5, 2, 10)
    cone2 = setup_cone(3, 2, 5)

    # position and orientation
    cone2.SetCenter(-2, 0, 0)
    cone1.SetCenter(2, 0, 0)
    cone1.SetDirection((0.0, 1.0, 0.0))
    cone2.SetDirection((0.0, 1.0, 0.0))

    # define mapper
    mappers = []
    for cone in (cone1, cone2):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cone.GetOutputPort())
        mappers.append(mapper)

    pipeline = VisualisationPipeline(mappers)
    pipeline.run()
