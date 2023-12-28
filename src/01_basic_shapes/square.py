from typing import Tuple

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_square_mapper(
    origin: Tuple[float, float, float],
    point1: Tuple[float, float, float],
    point2: Tuple[float, float, float],
) -> vtk.vtkPolyDataMapper:
    square = vtk.vtkPlaneSource()
    square.SetOrigin(origin)
    square.SetPoint1(point1)
    square.SetPoint2(point2)

    square_mapper = vtk.vtkPolyDataMapper()
    square_mapper.SetInputConnection(square.GetOutputPort())

    return square_mapper


if __name__ == "__main__":
    # Create squares with different positions and sizes
    square_mapper1 = create_square_mapper(
        origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 1.0, 0.0)
    )
    square_mapper2 = create_square_mapper(
        origin=(3.0, 0.0, 0.0), point1=(4.0, 0.0, 0.0), point2=(3.0, 1.0, 0.0)
    )

    # Display the squares
    pipeline = VisualisationPipeline(mappers=[square_mapper1, square_mapper2])
    pipeline.run()
