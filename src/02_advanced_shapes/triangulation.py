import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

# can be read from stl file


def generate_data():
    polydata = vtk.vtkSphereSource()
    polydata.SetRadius(10)
    polydata.SetThetaResolution(5)
    polydata.SetPhiResolution(3)
    polydata.Update()
    return polydata


def remove_double_vertices(polydata):
    clean_polydata = vtk.vtkCleanPolyData()
    clean_polydata.SetInputData(polydata.GetOutput())
    clean_polydata.Update()
    return clean_polydata


def create_normals(polydata):
    normals = vtk.vtkPolyDataNormals()
    normals.SetComputeCellNormals(1)
    normals.SetInputData(polydata.GetOutput())
    normals.SplittingOff()
    normals.Update()
    return normals


def triangulate_new_connecting_faces(normals):
    polydata = vtk.vtkPolyData()
    polydata.DeepCopy(normals.GetOutput())
    triangle_filter = vtk.vtkTriangleFilter()
    triangle_filter.SetInputData(polydata)
    triangle_filter.Update()
    return triangle_filter


if __name__ == "__main__":
    polydata = generate_data()
    polydata = remove_double_vertices(polydata)
    normals = create_normals(polydata)
    triangle_filter = triangulate_new_connecting_faces(normals)
    data = triangle_filter.GetOutput()
    # visualize
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(data)

    pipeline = VisualisationPipeline([mapper], edges_visible=True)
    pipeline.run()
