from pyvtk import *
from vtkmodules.vtkIOLegacy import vtkPolyDataReader
from vtkmodules.vtkRenderingCore import vtkDataSetMapper

from src.io.read_vtk import read_vtk
from src.visualization_pipeline.simple_pipeline import VisualisationPipeline

file_name = "example.vtk"


def generate_data():
    structure = PolyData(
        points=[
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        polygons=[
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [2, 3, 7, 6],
            [0, 4, 7, 3],
            [1, 2, 6, 5],
        ],
    )

    point_data = PointData(Scalars([0, 1, 2, 3, 4, 5, 6, 7], name="sample_scalars"))

    cell_data = CellData(
        Scalars([0, 1, 2, 3, 4, 5], name="cell_scalars"),
        Normals(
            [[0, 0, 1], [0, 0, 1], [0, -1, 0], [0, 1, 0], [-1, 0, 0], [1, 0, 0]],
            name="cell_normals",
        ),
        Field(
            "FieldData",
            cellIds=[[0], [1], [2], [3], [4], [5]],
            faceAttributes=[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6]],
        ),
    )

    vtk_data = VtkData(structure, point_data, cell_data)
    return vtk_data


if __name__ == "__main__":
    vtk_data = generate_data()
    vtk_data.tofile(file_name, "ascii")
    # vtk_data.tofile(f"{file_name}b", "binary")

    output = read_vtk(file_name, vtkPolyDataReader)

    mapper = vtkDataSetMapper()
    mapper.SetInputData(output)
    mapper.SetScalarRange(output.GetScalarRange())

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
