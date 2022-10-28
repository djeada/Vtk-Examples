import vtk

from src.visualization_pipeline.simple_pipeline import VisualisationPipeline


if __name__ == "__main__":
    ## glyph
    glyph = vtk.vtkGlyph2D()

    ## locations
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(3, 3, 0)
    points.InsertNextPoint(6, 6, 0)
    poly_data = vtk.vtkPolyData()
    poly_data.SetPoints(points)
    glyph.SetInputData(poly_data)

    ## glyph source
    point_sources = vtk.vtkGlyphSource2D()
    point_sources.SetGlyphTypeToCircle()
    glyph.SetSourceConnection(point_sources.GetOutputPort())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.Update()

    pipeline = VisualisationPipeline([mapper])
    pipeline.run()
