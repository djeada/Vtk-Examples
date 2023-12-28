"""
This module demonstrates the workflow of creating and visualizing a set of randomly positioned spheres within an enclosing box using the Visualization Toolkit (VTK). It highlights the use of VTK's capabilities for generating and displaying 3D objects in a defined space, showcasing fundamental concepts in computational geometry and 3D visualization.

Workflow Overview:

1. Sphere Creation (create_random_sphere):
   - Generates individual spheres with specified centers and radii.
   - Illustrates how to create basic 3D objects (spheres) in VTK.

2. Random Spheres Mapper (create_random_spheres_mapper):
   - Constructs multiple spheres, each positioned randomly within a specified box size.
   - Utilizes vtkAppendPolyData to combine multiple vtkSphereSources into a single data structure.
   - Creates a vtkPolyDataMapper for the set of spheres, preparing them for visualization.

3. Enclosing Box Mapper (create_enclosing_box_mapper):
   - Generates an enclosing box (using vtkCubeSource) that defines the boundary for the random spheres.
   - Applies vtkOutlineFilter to create an outline of the box, enhancing spatial understanding.
   - Produces a vtkPolyDataMapper for the box, allowing it to be rendered alongside the spheres.

4. Visualization Pipeline:
   - The module uses a predefined 'VisualisationPipeline' class to manage the rendering process.
   - Both the random spheres and the enclosing box are visualized together, providing context and enhancing the perception of the 3D space.
"""

import random

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_random_sphere(center, radius):
    """
    Create a sphere source with specified center and radius.

    Args:
    center (tuple): The center of the sphere.
    radius (float): The radius of the sphere.

    Returns:
    vtkSphereSource: A VTK sphere source.
    """
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(*center)
    sphere.SetRadius(radius)
    sphere.Update()
    return sphere


def create_random_spheres_mapper(n_spheres, box_size, sphere_radius=0.2):
    """
    Create a mapper for a set of randomly positioned spheres within a box.

    Args:
    n_spheres (int): Number of spheres to create.
    box_size (float): Size of the box in which spheres are contained.
    sphere_radius (float): Radius of each sphere.

    Returns:
    vtkPolyDataMapper: A mapper for the random spheres.
    """
    append_filter = vtk.vtkAppendPolyData()

    for _ in range(n_spheres):
        center = (
            random.uniform(-box_size / 2, box_size / 2),
            random.uniform(-box_size / 2, box_size / 2),
            random.uniform(-box_size / 2, box_size / 2),
        )
        sphere = create_random_sphere(center, sphere_radius)
        append_filter.AddInputConnection(sphere.GetOutputPort())

    append_filter.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(append_filter.GetOutputPort())
    return mapper


def create_enclosing_box_mapper(box_size):
    """
    Create a mapper for an enclosing box.

    Args:
    box_size (float): Size of the box.

    Returns:
    vtkPolyDataMapper: A mapper for the enclosing box.
    """
    box = vtk.vtkCubeSource()
    box.SetXLength(box_size)
    box.SetYLength(box_size)
    box.SetZLength(box_size)
    box.Update()

    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(box.GetOutputPort())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(outline.GetOutputPort())
    return mapper


if __name__ == "__main__":
    n_spheres = 50
    box_size = 5.0

    spheres_mapper = create_random_spheres_mapper(n_spheres, box_size)
    box_mapper = create_enclosing_box_mapper(box_size)

    pipeline = VisualisationPipeline(mappers=[spheres_mapper, box_mapper])
    pipeline.run()
