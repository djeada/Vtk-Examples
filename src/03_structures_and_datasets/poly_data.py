"""
This module demonstrates the creation and visualization of a simple triangle using the Visualization Toolkit (VTK). It is structured into several distinct steps, each utilizing specific VTK components for 3D graphics and visualization:

1. Point Creation (vtkPoints):
   - Responsible for defining and storing the coordinates of vertices in a 3D space.
   - In this example, three points are created to represent the vertices of a triangle.

2. Triangle Geometry (vtkTriangle):
   - Defines the geometry of a triangle by specifying the connectivity between the points.
   - It uses the points' indices to create a triangular cell.

3. Cell Array (vtkCellArray):
   - A collection of cells (in this case, just one triangle).
   - Used to store the triangle geometry and facilitates the handling of multiple cells.

4. PolyData Creation (vtkPolyData):
   - A data object that holds geometric and topological data.
   - In this example, it stores the points and the triangle cell, representing the complete geometry of the triangle.

5. Data Mapping (vtkPolyDataMapper):
   - Maps the polydata (geometric data) to graphical primitives.
   - Essential for rendering the data in a window.

6. Actor Creation (vtkActor):
   - Represents an entity in the rendering process.
   - It is linked with the data mapper and controls the rendering of the polydata.

7. Renderer and Render Window:
   - The renderer (vtkRenderer) adds the actor to the scene.
   - The render window (vtkRenderWindow) displays the rendered graphics.

8. Interaction:
   - The render window interactor (vtkRenderWindowInteractor) allows user interaction with the visualization, such as rotating and zooming.

The module serves as an educational example to illustrate the basic components and workflow in VTK for creating and displaying simple geometric shapes in 3D. It highlights the use of polydata as a versatile structure for representing and manipulating geometric data in VTK, which is a cornerstone for more complex 3D visualizations in scientific and engineering applications.
"""
import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_points():
    """
    Create and return a vtkPoints object with predefined points.
    """
    points = vtk.vtkPoints()
    points.InsertNextPoint([0.0, 0.0, 0.0])  # Point 1
    points.InsertNextPoint([1.0, 0.0, 0.0])  # Point 2
    points.InsertNextPoint([0.5, 0.866, 0.0])  # Point 3
    return points


def create_triangle():
    """
    Create and return a vtkTriangle.
    """
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 0)  # Vertex 1
    triangle.GetPointIds().SetId(1, 1)  # Vertex 2
    triangle.GetPointIds().SetId(2, 2)  # Vertex 3
    return triangle


def create_polydata(points, triangle):
    """
    Create and return a vtkPolyData object containing the points and triangle.
    """
    triangles = vtk.vtkCellArray()
    triangles.InsertNextCell(triangle)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(triangles)
    return polydata


def main():
    points = create_points()
    triangle = create_triangle()
    polydata = create_polydata(points, triangle)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    pipeline = VisualisationPipeline(mappers=[mapper])
    pipeline.run()


if __name__ == "__main__":
    main()
