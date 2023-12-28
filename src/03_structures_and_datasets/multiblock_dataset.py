"""
A multiblock dataset in VTK is a composite dataset that can store other datasets including other composite datasets. 
  - It is used to represent complex data where different regions can be meshed differently. 
  - It is also used to group multiple datasets into one.

Workflow Overview:

1. Geometric Source Creation:
   - The process begins with the dynamic creation of various geometric shapes such as spheres, cylinders, cones, and cubes. 
   - Each shape is generated based on specific parameters (like radius, height, and dimensions), illustrating how VTK can be used to create a wide range of 3D objects.

2. Multiblock Dataset Assembly:
   - Once individual geometric shapes are created, they are assembled into a vtkMultiBlockDataSet. 
   - This dataset serves as a container that can hold multiple disparate data sources, showcasing its capability to manage complex collections of 3D objects.

3. Conversion to PolyData:
   - The vtkMultiBlockDataSet is then converted into vtkPolyData, a versatile format in VTK used for rendering and visualizing geometric data.
   - This step is crucial as it transforms the multiblock dataset into a format that can be easily rendered by VTK's visualization pipeline.

4. Visualization Pipeline:
   - Finally, the polydata is fed into a pre-defined visualization pipeline, which is responsible for rendering the 3D objects on the screen.
   - This pipeline not only displays the assembled 3D shapes but also allows for interactive exploration of the scene, such as zooming and rotating the objects.
"""

import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_geometric_source(source_type, **params):
    """
    Create and return a geometric source based on the specified type and parameters.
    Supported types: 'Sphere', 'Cylinder', 'Cone', 'Cube'.
    """
    if source_type == "Sphere":
        source = vtk.vtkSphereSource()
        source.SetRadius(params.get("radius", 1.0))
        source.SetCenter(params.get("center", [0.0, 0.0, 0.0]))
    elif source_type == "Cylinder":
        source = vtk.vtkCylinderSource()
        source.SetHeight(params.get("height", 2.0))
        source.SetRadius(params.get("radius", 0.5))
    elif source_type == "Cone":
        source = vtk.vtkConeSource()
        source.SetHeight(params.get("height", 2.0))
        source.SetRadius(params.get("radius", 0.5))
    elif source_type == "Cube":
        source = vtk.vtkCubeSource()
        source.SetXLength(params.get("length", 1.0))
        source.SetYLength(params.get("width", 1.0))
        source.SetZLength(params.get("height", 1.0))
    else:
        raise ValueError("Unsupported source type")

    source.Update()
    return source.GetOutput()


def create_multiblock_dataset(sources):
    """
    Create a multiblock dataset from a list of VTK data sources.
    """
    multiBlock = vtk.vtkMultiBlockDataSet()
    multiBlock.SetNumberOfBlocks(len(sources))
    for i, source in enumerate(sources):
        multiBlock.SetBlock(i, source)
    return multiBlock


def convert_to_polydata(multiBlock):
    """
    Convert the multiblock dataset to polydata.
    """
    geometryFilter = vtk.vtkCompositeDataGeometryFilter()
    geometryFilter.SetInputDataObject(multiBlock)
    geometryFilter.Update()
    return geometryFilter.GetOutput()


def main():
    # Create various geometric sources
    sphere = create_geometric_source("Sphere", radius=1.0)
    cylinder = create_geometric_source("Cylinder", height=2.0, radius=0.5)
    cone = create_geometric_source("Cone", height=2.0, radius=0.5)
    cube = create_geometric_source("Cube", length=1.0, width=1.0, height=1.0)

    # Create a multiblock dataset containing all sources
    sources = [sphere, cylinder, cone, cube]
    multiBlock = create_multiblock_dataset(sources)
    polydata = convert_to_polydata(multiBlock)

    # Create a mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    # Run the visualization pipeline
    pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
    pipeline.run()


if __name__ == "__main__":
    main()
