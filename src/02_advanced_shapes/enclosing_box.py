import vtk
import random
from src.simple_pipeline import VisualisationPipeline


def create_random_spheres_mapper(n_spheres, box_size):
    append_filter = vtk.vtkAppendPolyData()

    for _ in range(n_spheres):
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(
            random.uniform(-box_size / 2, box_size / 2),
            random.uniform(-box_size / 2, box_size / 2),
            random.uniform(-box_size / 2, box_size / 2),
        )
        sphere.SetRadius(0.2)
        append_filter.AddInputConnection(sphere.GetOutputPort())

    append_filter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(append_filter.GetOutputPort())

    return mapper


def create_enclosing_box_mapper(input_mapper):
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(input_mapper.GetInputConnection(0, 0))

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(outline.GetOutputPort())

    return mapper


if __name__ == "__main__":
    n_spheres = 50
    box_size = 5.0

    spheres_mapper = create_random_spheres_mapper(n_spheres, box_size)
    box_mapper = create_enclosing_box_mapper(spheres_mapper)

    pipeline = VisualisationPipeline(mappers=[spheres_mapper, box_mapper])
    pipeline.run()
