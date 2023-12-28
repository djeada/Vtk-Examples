import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_circle_mapper(radius, center, resolution=100):
    circle = vtk.vtkRegularPolygonSource()
    circle.SetNumberOfSides(resolution)
    circle.SetRadius(radius)
    circle.SetCenter(center)
    circle.SetNormal(0, 0, 1)

    circle_mapper = vtk.vtkPolyDataMapper()
    circle_mapper.SetInputConnection(circle.GetOutputPort())

    return circle_mapper


if __name__ == "__main__":
    # Create a circle with different radius and center
    circle_mapper1 = create_circle_mapper(radius=1, center=(0, 0, 0))
    circle_mapper2 = create_circle_mapper(radius=2, center=(4, 0, 0))

    # Display the circles
    pipeline = VisualisationPipeline(mappers=[circle_mapper1, circle_mapper2])
    pipeline.run()
