"""
A multiblock dataset in VTK is a composite dataset that can store other datasets including other composite datasets. It is used to represent complex data where different regions can be meshed differently. It is also used to group multiple datasets into one.
"""

import vtk

from src.simple_pipeline import VisualisationPipeline

# Create a sphere source
sphereSource = vtk.vtkSphereSource()
sphereSource.Update()

# Create a cylinder source
cylinderSource = vtk.vtkCylinderSource()
cylinderSource.Update()

# Create a multiblock dataset and add the sphere and the cylinder
multiBlock = vtk.vtkMultiBlockDataSet()
multiBlock.SetNumberOfBlocks(2)
multiBlock.SetBlock(0, sphereSource.GetOutput())
multiBlock.SetBlock(1, cylinderSource.GetOutput())

# Use vtkCompositeDataGeometryFilter to convert the multiblock dataset to polydata
geometryFilter = vtk.vtkCompositeDataGeometryFilter()
geometryFilter.SetInputDataObject(multiBlock)
geometryFilter.Update()

# Create a mapper and actor for the polydata
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(geometryFilter.GetOutputPort())

pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
pipeline.run()
