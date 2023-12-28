import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_arrow_glyphs(vector_field):
    # Create arrow source
    arrow_source = vtk.vtkArrowSource()

    # Create glyph points and vectors
    glyph_points = vtk.vtkPoints()
    glyph_vectors = vtk.vtkDoubleArray()
    glyph_vectors.SetNumberOfComponents(3)
    glyph_vectors.SetNumberOfTuples(vector_field.shape[0] * vector_field.shape[1])

    idx = 0
    for i in range(vector_field.shape[0]):
        for j in range(vector_field.shape[1]):
            x, y, z = i, j, 0
            u, v, w = vector_field[i, j]
            glyph_points.InsertPoint(idx, x, y, z)
            glyph_vectors.SetTuple(idx, (u, v, w))
            idx += 1

    # Create glyph polydata
    glyph_polydata = vtk.vtkPolyData()
    glyph_polydata.SetPoints(glyph_points)
    glyph_polydata.GetPointData().SetVectors(glyph_vectors)

    # Create the glyph3D filter
    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetSourceConnection(arrow_source.GetOutputPort())
    glyph3D.SetInputData(glyph_polydata)
    glyph3D.SetVectorModeToUseVector()
    glyph3D.SetScaleModeToScaleByVector()
    glyph3D.SetScaleFactor(0.5)

    glyph_mapper = vtk.vtkPolyDataMapper()
    glyph_mapper.SetInputConnection(glyph3D.GetOutputPort())

    return glyph_mapper


def generate_vector_field(size, scale):
    import numpy as np

    x, y = np.meshgrid(
        np.linspace(-scale, scale, size), np.linspace(-scale, scale, size)
    )
    u = -y
    v = x
    vector_field = np.stack((u, v, np.zeros_like(u)), axis=-1)
    return vector_field


if __name__ == "__main__":
    vector_field = generate_vector_field(size=11, scale=5)
    glyph_mapper = create_arrow_glyphs(vector_field)

    # Display the glyphs
    pipeline = VisualisationPipeline(mappers=[glyph_mapper])
    pipeline.run()
