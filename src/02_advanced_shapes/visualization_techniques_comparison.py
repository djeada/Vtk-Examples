import random
import vtk


def create_render_window():
    render_window = vtk.vtkRenderWindow()
    render_window.SetSize(800, 600)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    return render_window, interactor


def create_data_source():
    volume_source = vtk.vtkRTAnalyticSource()

    vector_source = vtk.vtkPointSource()
    vector_source.SetNumberOfPoints(10)
    vector_source.Update()

    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetNumberOfTuples(vector_source.GetOutput().GetNumberOfPoints())

    for i in range(vectors.GetNumberOfTuples()):
        vectors.SetTuple3(i, 1, 2, 3)

    vector_source.GetOutput().GetPointData().SetVectors(vectors)
    return volume_source, vector_source


def create_structured_grid():
    grid_dimensions = [10, 10, 10]

    structured_grid = vtk.vtkStructuredGrid()
    structured_grid.SetDimensions(grid_dimensions)

    points = vtk.vtkPoints()

    for z in range(grid_dimensions[2]):
        for y in range(grid_dimensions[1]):
            for x in range(grid_dimensions[0]):
                points.InsertNextPoint(x, y, z)

    structured_grid.SetPoints(points)
    return structured_grid


def create_point_source():
    point_source = vtk.vtkPointSource()
    point_source.SetNumberOfPoints(50)
    point_source.SetRadius(3)
    point_source.Update()
    return point_source


def create_visualizations(
    render_window,
    xmins,
    ymins,
    xmaxs,
    ymaxs,
    volume_source,
    structured_grid,
    point_source,
    vector_source,
):
    for i in range(4):
        renderer = vtk.vtkRenderer()
        render_window.AddRenderer(renderer)
        renderer.SetViewport(xmins[i], ymins[i], xmaxs[i], ymaxs[i])

        if i == 0:  # Volume rendering
            create_volume_rendering(renderer, volume_source)
        elif i == 1:  # Streamlines
            create_streamlines(renderer, structured_grid, point_source)
        elif i == 2:  # Glyphs
            create_glyphs(renderer, vector_source)
        elif i == 3:  # Contouring
            create_contouring(renderer, volume_source)


def create_volume_rendering(renderer, volume_source):
    volume_mapper = vtk.vtkSmartVolumeMapper()
    volume_mapper.SetInputConnection(volume_source.GetOutputPort())
    volume = vtk.vtkVolume()
    volume.SetMapper(volume_mapper)
    renderer.AddVolume(volume)

    # Add text annotation
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("Volume Rendering")
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.GetTextProperty().SetColor(1, 1, 1)  # White text
    renderer.AddActor2D(text_actor)


def create_streamlines(renderer, structured_grid, point_source):
    vector_data_array = vtk.vtkDoubleArray()
    vector_data_array.SetNumberOfComponents(3)
    vector_data_array.SetNumberOfTuples(structured_grid.GetNumberOfPoints())

    for j in range(structured_grid.GetNumberOfPoints()):
        vector_data_array.SetTuple3(
            j, random.random(), random.random(), random.random()
        )

    vector_data_array.SetName("VectorField")
    structured_grid.GetPointData().SetVectors(vector_data_array)

    stream_tracer = vtk.vtkStreamTracer()
    stream_tracer.SetInputData(structured_grid)
    stream_tracer.SetSourceConnection(point_source.GetOutputPort())
    stream_tracer.SetMaximumPropagation(500)
    stream_tracer.SetInitialIntegrationStep(0.1)
    stream_tracer.SetIntegratorType(2)

    tube_filter = vtk.vtkTubeFilter()
    tube_filter.SetInputConnection(stream_tracer.GetOutputPort())
    tube_filter.SetRadius(0.1)
    tube_filter.SetNumberOfSides(6)

    streamline_mapper = vtk.vtkPolyDataMapper()
    streamline_mapper.SetInputConnection(tube_filter.GetOutputPort())
    streamline_actor = vtk.vtkActor()
    streamline_actor.SetMapper(streamline_mapper)

    renderer.AddActor(streamline_actor)

    # Add text annotation
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("Streamlines")
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.GetTextProperty().SetColor(1, 1, 1)  # White text
    renderer.AddActor2D(text_actor)


def create_glyphs(renderer, vector_source):
    # Create a glyph source and mapper
    glyph_source = vtk.vtkArrowSource()
    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetSourceConnection(glyph_source.GetOutputPort())
    glyph3D.SetInputConnection(vector_source.GetOutputPort())
    glyph3D.SetVectorModeToUseVector()
    glyph3D.SetScaleModeToScaleByVector()
    glyph3D.SetScaleFactor(0.1)
    glyph_mapper = vtk.vtkPolyDataMapper()
    glyph_mapper.SetInputConnection(glyph3D.GetOutputPort())
    # Create an actor
    glyph_actor = vtk.vtkActor()
    glyph_actor.SetMapper(glyph_mapper)
    renderer.AddActor(glyph_actor)

    # Add text annotation
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("Glyphs")
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.GetTextProperty().SetColor(1, 1, 1)  # White text
    renderer.AddActor2D(text_actor)


def create_contouring(renderer, volume_source):
    contour_filter = vtk.vtkContourFilter()
    contour_filter.SetInputConnection(volume_source.GetOutputPort())
    contour_filter.SetValue(0, 150)

    contour_mapper = vtk.vtkPolyDataMapper()
    contour_mapper.SetInputConnection(contour_filter.GetOutputPort())

    contour_actor = vtk.vtkActor()
    contour_actor.SetMapper(contour_mapper)
    renderer.AddActor(contour_actor)

    # Add text annotation
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("Contouring")
    text_actor.GetTextProperty().SetFontSize(24)
    text_actor.GetTextProperty().SetColor(1, 1, 1)  # White text
    renderer.AddActor2D(text_actor)


def main():
    xmins = [0, 0.5, 0, 0.5]
    xmaxs = [0.5, 1, 0.5, 1]
    ymins = [0, 0, 0.5, 0.5]
    ymaxs = [0.5, 0.5, 1, 1]

    render_window, interactor = create_render_window()
    volume_source, vector_source = create_data_source()
    structured_grid = create_structured_grid()
    point_source = create_point_source()

    create_visualizations(
        render_window,
        xmins,
        ymins,
        xmaxs,
        ymaxs,
        volume_source,
        structured_grid,
        point_source,
        vector_source,
    )

    interactor.Initialize()
    interactor.Start()


if __name__ == "__main__":
    main()
