import vtk


class VisualisationPipeline:
    def __init__(
        self,
        mappers,
        point_size=1,
        edges_visible=False,
        background_color=(0, 0, 0),
        window_title="VTK Visualization",
    ):
        """
        Initialize the visualization pipeline with specified parameters.

        Args:
        mappers (list of vtkPolyDataMapper): The data mappers for rendering.
        point_size (int): Size of points in the rendered actors.
        edges_visible (bool): Flag to render edges of the actors.
        background_color (tuple): Background color of the renderer.
        window_title (str): Title of the visualization window.
        """
        self.mappers = mappers
        self.point_size = point_size
        self.edges_visible = edges_visible
        self.background_color = background_color
        self.window_title = window_title

    def create_actors(self):
        """
        Create actors based on the mappers provided.

        Returns:
        list of vtkActor: List of actors for rendering.
        """
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
        """
        Create a renderer and add actors to it.

        Args:
        actors (list of vtkActor): Actors to be rendered.

        Returns:
        vtkRenderer: Configured renderer.
        """
        renderer = vtk.vtkRenderer()
        for actor in actors:
            renderer.AddActor(actor)
        renderer.SetBackground(*self.background_color)
        return renderer

    def create_window(self, renderer):
        """
        Create a rendering window.

        Args:
        renderer (vtkRenderer): Renderer to be displayed in the window.

        Returns:
        vtkRenderWindow: Configured rendering window.
        """
        window = vtk.vtkRenderWindow()
        window.AddRenderer(renderer)
        window.SetSize(800, 600)
        # window.SetTitle(self.window_title)
        return window

    def create_interactor(self, window):
        """
        Create an interactor for the rendering window.

        Args:
        window (vtkRenderWindow): Rendering window.

        Returns:
        vtkRenderWindowInteractor: Configured window interactor.
        """
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(window)
        return interactor

    def run(self):
        """
        Run the visualization pipeline.
        """
        actors = self.create_actors()
        renderer = self.create_renderer(actors)
        window = self.create_window(renderer)
        interactor = self.create_interactor(window)
        interactor.Start()
