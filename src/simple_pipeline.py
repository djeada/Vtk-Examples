import vtk


class VisualisationPipeline:
    def __init__(self, mappers, point_size=1, edges_visible=False):
        self.mappers = mappers
        self.point_size = point_size
        self.edges_visible = edges_visible

    def create_actors(self):
        actors = []
        for mapper in self.mappers:
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetPointSize(self.point_size)
            if self.edges_visible:
                actor.GetProperty().SetEdgeVisibility(1)
                actor.GetProperty().SetOpacity(0.5)
                actor.GetProperty().BackfaceCullingOn()
            actors.append(actor)

        return actors

    def create_renderer(self, actors):
        renderer = vtk.vtkRenderer()
        for actor in actors:
            renderer.AddActor(actor)
        renderer.SetBackground(0, 0, 0)
        return renderer

    def create_window(self, renderer):
        window = vtk.vtkRenderWindow()
        window.AddRenderer(renderer)
        window.SetSize(800, 600)
        return window

    def create_interactor(self, window):
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(window)
        interactor.Initialize()
        return interactor

    def run(self):
        actors = self.create_actors()
        renderer = self.create_renderer(actors)
        window = self.create_window(renderer)
        interactor = self.create_interactor(window)
        interactor.Start()
