"""
This module demonstrates the process of 3D triangulation and its visualization using the Visualization Toolkit (VTK). Triangulation is a fundamental operation in computational geometry, where a surface or volume is divided into simpler units, specifically triangles in 2D or tetrahedrons in 3D. These simple shapes facilitate easier rendering and enable complex geometrical shapes to be approximated with high accuracy.

Workflow:
1. Generate a set of random points in 3D space.
2. Use the vtkDelaunay3D filter to perform Delaunay triangulation on these points.
3. Convert the output vtkUnstructuredGrid to vtkPolyData using vtkGeometryFilter.
4. Pass this vtkPolyData to vtkTriangleFilter to process and prepare it for visualization.
5. Create a vtkPolyDataMapper and vtkActor to render the triangulated geometry.
6. Set up a VTK renderer, render window, and interactor to display the result.

Mathematics:
- Delaunay Triangulation: For a given set of points in a space (R^3), Delaunay triangulation generates a mesh of tetrahedra such that no point is inside the circum-hypersphere of any tetrahedra. It maximizes the minimum angle of all the angles of the triangles in the triangulation, helping avoid skinny triangles.
- 3D Space Representation: Points are randomly generated in a 3D space, each defined by (x, y, z) coordinates.
- vtkDelaunay3D: This VTK filter computes the Delaunay triangulation in 3D. It outputs a vtkUnstructuredGrid representing the tetrahedral mesh.
- vtkGeometryFilter: Converts the vtkUnstructuredGrid (tetrahedral mesh) to vtkPolyData, which is a suitable format for further processing and rendering.
- vtkTriangleFilter: Processes the vtkPolyData to ensure that the geometry is represented purely with triangles, a requirement for many rendering and processing operations.

"""

import random

import vtk


class Triangulation3D:
    def __init__(self, num_points=25):
        self.num_points = num_points
        self.points = vtk.vtkPoints()
        self.generate_random_points()
        self.triangulate()

    def generate_random_points(self):
        for _ in range(self.num_points):
            x, y, z = (
                random.uniform(-1, 1),
                random.uniform(-1, 1),
                random.uniform(-1, 1),
            )
            self.points.InsertNextPoint(x, y, z)

    def triangulate(self):
        input_polydata = vtk.vtkPolyData()
        input_polydata.SetPoints(self.points)

        delaunay = vtk.vtkDelaunay3D()
        delaunay.SetInputData(input_polydata)
        delaunay.Update()

        # Convert vtkUnstructuredGrid to vtkPolyData
        geometry_filter = vtk.vtkGeometryFilter()
        geometry_filter.SetInputConnection(delaunay.GetOutputPort())
        geometry_filter.Update()

        self.triangle_filter = vtk.vtkTriangleFilter()
        self.triangle_filter.SetInputConnection(geometry_filter.GetOutputPort())
        self.triangle_filter.PassLinesOn()
        self.triangle_filter.PassVertsOn()
        self.triangle_filter.Update()

    def get_actor(self):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.triangle_filter.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetRepresentationToWireframe()

        return actor


def main():
    triangulation = Triangulation3D()

    renderer = vtk.vtkRenderer()
    renderer.AddActor(triangulation.get_actor())
    renderer.SetBackground(0, 0, 0)

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    interactor.Initialize()
    interactor.Start()


if __name__ == "__main__":
    main()
