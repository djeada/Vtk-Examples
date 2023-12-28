import vtk

from src.simple_pipeline import VisualisationPipeline


def create_triangle_mapper(points):
    triangle = vtk.vtkTriangle()
    triangle.GetPointIds().SetId(0, 0)
    triangle.GetPointIds().SetId(1, 1)
    triangle.GetPointIds().SetId(2, 2)

    points_array = vtk.vtkPoints()
    for point in points:
        points_array.InsertNextPoint(point)

    triangles = vtk.vtkCellArray()
    triangles.InsertNextCell(triangle)

    triangle_polydata = vtk.vtkPolyData()
    triangle_polydata.SetPoints(points_array)
    triangle_polydata.SetPolys(triangles)

    triangle_mapper = vtk.vtkPolyDataMapper()
    triangle_mapper.SetInputData(triangle_polydata)

    return triangle_mapper


if __name__ == "__main__":
    # Create triangles with different vertex positions
    triangle_mapper1 = create_triangle_mapper([(0, 0, 0), (1, 0, 0), (0.5, 1, 0)])
    triangle_mapper2 = create_triangle_mapper([(3, 0, 0), (4, 0, 0), (3.5, 1, 0)])

    # Display the triangles
    pipeline = VisualisationPipeline(mappers=[triangle_mapper1, triangle_mapper2])
    pipeline.run()
