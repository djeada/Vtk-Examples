import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_cone_mapper(radius, height, center, resolution=100):
    cone = vtk.vtkConeSource()
    cone.SetRadius(radius)
    cone.SetHeight(height)
    cone.SetCenter(center)
    cone.SetDirection(0, 0, 1)
    cone.SetResolution(resolution)

    cone_mapper = vtk.vtkPolyDataMapper()
    cone_mapper.SetInputConnection(cone.GetOutputPort())

    return cone_mapper


if __name__ == "__main__":
    # Create cones with different radius, height, and center
    cone_mapper1 = create_cone_mapper(radius=1, height=2, center=(0, 0, 0))
    cone_mapper2 = create_cone_mapper(radius=1.5, height=3, center=(4, 0, 0))

    # Display the cones
    pipeline = VisualisationPipeline(mappers=[cone_mapper1, cone_mapper2])
    pipeline.run()
