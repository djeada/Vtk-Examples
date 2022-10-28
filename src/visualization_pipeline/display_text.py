# create a cube and show it with text label below

import vtk


def setup_data():
    cube = vtk.vtkCubeSource()
    cube.SetXLength(1)
    cube.SetYLength(1)
    cube.SetZLength(1)
    return cube


def setup_mapper(cube):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cube.GetOutputPort())
    return mapper


def setup_actor(mapper):
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(0, 0, 0)
    return actor


def setup_renderer(actor):
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    return renderer


def setup_render_window(renderer):
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    return render_window


def setup_interactor(render_window):
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    return interactor


# add the actor to the scene
def add_actor(renderer, actor):
    renderer.AddActor(actor)
    renderer.SetBackground(0.1, 0.2, 0.4)  # Background color dark blue


# create a text property
def create_text_property():
    textProperty = vtk.vtkTextProperty()
    textProperty.SetFontSize(18)
    textProperty.SetJustificationToCentered()
    textProperty.BoldOn()
    textProperty.ItalicOn()
    textProperty.ShadowOn()
    return textProperty


def create_text_mapper(textProperty, text="Hello World"):
    textMapper = vtk.vtkTextMapper()
    textMapper.SetInput(text)
    textMapper.SetTextProperty(textProperty)
    return textMapper


# create an actor2D
def create_text_actor(textMapper):
    textActor = vtk.vtkActor2D()
    textActor.SetMapper(textMapper)
    textActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
    textActor.GetPositionCoordinate().SetValue(0.5, 0.1)
    return textActor


# add the actor to the scene
def add_text_actor(renderer, textActor):
    renderer.AddActor2D(textActor)


if __name__ == "__main__":
    cube = setup_data()
    mapper = setup_mapper(cube)
    actor = setup_actor(mapper)
    renderer = setup_renderer(actor)
    render_window = setup_render_window(renderer)
    interactor = setup_interactor(render_window)

    textProperty = create_text_property()
    textMapper = create_text_mapper(textProperty)
    textActor = create_text_actor(textMapper)

    add_actor(renderer, actor)
    add_text_actor(renderer, textActor)

    # make the cube fit half of the window
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Zoom(0.5)

    # render and interact
    render_window.SetSize(500, 500)
    render_window.Render()
    interactor.Start()
