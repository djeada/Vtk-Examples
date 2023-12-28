from typing import Tuple

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_cylinder_mapper(
    radius: float,
    height: float,
    center: Tuple[float, float, float],
    resolution: int = 100,
) -> vtk.vtkPolyDataMapper:
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(radius)
    cylinder.SetHeight(height)
    cylinder.SetCenter(center)
    cylinder.SetResolution(resolution)

    cylinder_mapper = vtk.vtkPolyDataMapper()
    cylinder_mapper.SetInputConnection(cylinder.GetOutputPort())

    return cylinder_mapper


if __name__ == "__main__":
    # Create cylinders with different radius, height, and center
    cylinder_mapper1 = create_cylinder_mapper(
        radius=1.0, height=2.0, center=(0.0, 0.0, 0.0)
    )
    cylinder_mapper2 = create_cylinder_mapper(
        radius=1.5, height=3.0, center=(4.0, 0.0, 0.0)
    )

    # Display the cylinders
    pipeline = VisualisationPipeline(mappers=[cylinder_mapper1, cylinder_mapper2])
    pipeline.run()
