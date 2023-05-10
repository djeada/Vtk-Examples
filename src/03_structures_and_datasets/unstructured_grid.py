from pyvtk import VtkData, UnstructuredGrid, PointData, Scalars, Vectors
from vtk.vtkIOLegacy import vtkPolyDataReader
from vtk import (
    vtkDataSetMapper,
    vtkActor,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkUnstructuredGridReader,
)

from src.io.read_vtk import read_vtk
from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "example.vtk"


def generate_data():

    points = [
        [0, 0, 0],
        [1, 0, 0],
        [2, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [2, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [2, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
        [2, 1, 1],
        [0, 1, 2],
        [1, 1, 2],
        [2, 1, 2],
        [0, 1, 3],
        [1, 1, 3],
        [2, 1, 3],
        [0, 1, 4],
        [1, 1, 4],
        [2, 1, 4],
        [0, 1, 5],
        [1, 1, 5],
        [2, 1, 5],
        [0, 1, 6],
        [1, 1, 6],
        [2, 1, 6],
    ]
    vectors = [
        [1, 0, 0],
        [1, 1, 0],
        [0, 2, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 2, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 2, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 2, 0],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
        [0, 0, 1],
    ]
    vtk_data = VtkData(
        UnstructuredGrid(
            points,
            hexahedron=[
                [0, 1, 4, 3, 6, 7, 10, 9],
                [1, 2, 5, 4, 7, 8, 11, 10],
                [15, 16, 17, 19, 18, 14, 13, 12],
            ],
            tetra=[[6, 10, 9, 12], [5, 11, 10, 14]],
            polygon=[15, 16, 17, 14, 13, 12],
            triangle_strip=[18, 15, 19, 16, 20, 17],
            quad=[22, 23, 20, 19],
            triangle=[[21, 22, 18], [22, 19, 18]],
            line=[[26, 25], [3, 20]],
            vertex=[24],
        ),
        PointData(Vectors(vectors), Scalars(range(27))),
        "Unstructured Grid",
    )

    return vtk_data


if __name__ == "__main__":

    vtk_data = generate_data()
    vtk_data.tofile(file_name, "ascii")

    output = read_vtk(file_name)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
