"""
Fields in VTK can be attached to both points and cells. In this example, I will demonstrate how to attach a scalar field to points and a vector field to cells. We can create a simple situation where we assign random scalar values to points and random vector values to cells.
in VTK, you can store physical properties such as pressure, temperature, or any scalar or vector quantity associated with the points or cells in your dataset. These are often referred to as point data or cell data, respectively.

    Point Data: These are data attributes that are associated with the points of a dataset. For example, in a computational fluid dynamics (CFD) simulation, you might store the velocity vector at each point in the mesh as point data.

    Cell Data: These are data attributes that are associated with the cells of a dataset. For example, you might store the pressure or temperature within each cell of the mesh as cell data.

Workflow Overview:

1. Point Creation (create_points):
   - Initializes a set of points in 3D space, laying the groundwork for constructing geometric shapes.
   - In this case, points defining the corners of a quad are created.

2. Quad Creation (create_quad):
   - Utilizes the previously defined points to create a quad cell, a basic polygonal element.
   - Demonstrates how cells are constructed from points to form geometric shapes.

3. Scalar Field Attachment (attach_scalar_field):
   - Attaches a scalar field to the points, assigning a random scalar value to each point.
   - This highlights how to associate data with geometric entities in VTK.

4. Vector Field Attachment (attach_vector_field):
   - Similar to the scalar field, a random vector field is attached to the cells of the polydata.
   - Shows the capability of VTK to handle complex data associations, including vectorial data.

5. Polydata Assembly (create_polydata):
   - Integrates the points and the quad cell into a vtkPolyData object, a versatile format for representing and manipulating geometric data in VTK.

6. Visualization Process:
   - The vtkPolyDataMapper maps the polydata to graphical primitives for rendering.
   - A visualization pipeline, provided by the 'VisualisationPipeline' class, is employed to render the polydata, showcasing the scalar and vector fields attached to it.
"""
import vtk
import numpy as np

from src.simple_pipeline import VisualisationPipeline

def create_points():
    """
    Create a set of points for the polydata.
    """
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(1, 0, 0)
    points.InsertNextPoint(0, 1, 0)
    points.InsertNextPoint(1, 1, 0)
    return points

def create_quad(points):
    """
    Create a quad cell using the provided points.
    """
    quad = vtk.vtkQuad()
    # Assign points to the quad corners
    quad.GetPointIds().SetId(0, 0)
    quad.GetPointIds().SetId(1, 1)
    quad.GetPointIds().SetId(2, 3)
    quad.GetPointIds().SetId(3, 2)
    return quad

def attach_scalar_field(points):
    """
    Attach a random scalar field to the points.
    """
    scalars = vtk.vtkFloatArray()
    scalars.SetName("Scalars")
    for i in range(points.GetNumberOfPoints()):
        scalars.InsertNextValue(np.random.rand())
    return scalars

def attach_vector_field(cells):
    """
    Attach a random vector field to the cells.
    """
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetName("Vectors")
    for i in range(cells.GetNumberOfCells()):
        vectors.InsertNextTuple3(np.random.rand(), np.random.rand(), np.random.rand())
    return vectors

def create_polydata(points, quad):
    """
    Create polydata to hold the points and cells.
    """
    cells = vtk.vtkCellArray()
    cells.InsertNextCell(quad)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)
    return polydata, cells

def main():
    points = create_points()
    quad = create_quad(points)
    polydata, cells = create_polydata(points, quad)

    polydata.GetPointData().SetScalars(attach_scalar_field(points))
    polydata.GetCellData().SetVectors(attach_vector_field(cells))

    # Visualizing the fields
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.SetScalarModeToUsePointData()
    mapper.SetColorModeToMapScalars()
    mapper.SelectColorArray("Scalars")

    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()

if __name__ == "__main__":
    main()
