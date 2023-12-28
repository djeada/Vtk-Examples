from typing import List

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_cube_mapper(bounds: List[float]) -> vtk.vtkPolyDataMapper:
    cube = vtk.vtkCubeSource()
    cube.SetBounds(bounds)

    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube.GetOutputPort())

    return cube_mapper


if __name__ == "__main__":
    # Create cubes with different bounds
    cube_mapper1 = create_cube_mapper(bounds=[0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
    cube_mapper2 = create_cube_mapper(bounds=[2.0, 4.0, 0.0, 1.0, 0.0, 1.0])

    # Display the cubes
    pipeline = VisualisationPipeline(mappers=[cube_mapper1, cube_mapper2])
    pipeline.run()
