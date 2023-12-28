import vtk

from src.simple_pipeline import VisualisationPipeline


def create_cylinder_mapper(radius, height, center, resolution=100):
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
    cylinder_mapper1 = create_cylinder_mapper(radius=1, height=2, center=(0, 0, 0))
    cylinder_mapper2 = create_cylinder_mapper(radius=1.5, height=3, center=(4, 0, 0))

    # Display the cylinders
    pipeline = VisualisationPipeline(mappers=[cylinder_mapper1, cylinder_mapper2])
    pipeline.run()
