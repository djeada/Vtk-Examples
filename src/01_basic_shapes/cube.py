import vtk
from src.simple_pipeline import VisualisationPipeline


def create_cube_mapper(bounds):
    cube = vtk.vtkCubeSource()
    cube.SetBounds(bounds)

    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube.GetOutputPort())

    return cube_mapper


if __name__ == "__main__":
    # Create cubes with different bounds
    cube_mapper1 = create_cube_mapper(bounds=[0, 1, 0, 1, 0, 1])
    cube_mapper2 = create_cube_mapper(bounds=[2, 4, 0, 1, 0, 1])

    # Display the cubes
    pipeline = VisualisationPipeline(mappers=[cube_mapper1, cube_mapper2])
    pipeline.run()
