import unittest
from unittest.mock import Mock, patch

import vtk

from src.common.simple_pipeline import VisualisationPipeline


class PipelineTests(unittest.TestCase):
    def make_pipeline(self):
        source = vtk.vtkSphereSource()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        return VisualisationPipeline([mapper], window_title="Example sphere")

    def test_window_title_and_interaction_style(self):
        pipeline = self.make_pipeline()
        renderer = pipeline.create_renderer(pipeline.create_actors())
        window = pipeline.create_window(renderer)
        self.addCleanup(window.Finalize)
        self.assertEqual(window.GetWindowName(), "Example sphere")
        self.assertEqual(window.GetSize(), (800, 600))
        interactor = pipeline.create_interactor(window)
        self.assertIsInstance(
            interactor.GetInteractorStyle(), vtk.vtkInteractorStyleTrackballCamera
        )

    def test_run_frames_objects_and_initializes_before_starting(self):
        pipeline = self.make_pipeline()
        events = []
        window = Mock()
        interactor = Mock()
        window.Render.side_effect = lambda: events.append("render")
        window.Finalize.side_effect = lambda: events.append("finalize")
        interactor.Initialize.side_effect = lambda: events.append("initialize")
        interactor.Start.side_effect = lambda: events.append("start")
        with (
            patch.object(
                pipeline, "create_window", return_value=window
            ) as create_window,
            patch.object(pipeline, "create_interactor", return_value=interactor),
        ):
            pipeline.run()
        renderer = create_window.call_args.args[0]
        # Default camera position is (0, 0, 1); framing moves it to fit the sphere.
        self.assertGreater(renderer.GetActiveCamera().GetPosition()[2], 1)
        self.assertEqual(events, ["initialize", "render", "start", "finalize"])

    def test_run_finalizes_window_if_rendering_or_interaction_fails(self):
        for operation in ("Initialize", "Render", "Start"):
            with self.subTest(operation=operation):
                pipeline = self.make_pipeline()
                window = Mock()
                interactor = Mock()
                owner = window if operation == "Render" else interactor
                getattr(owner, operation).side_effect = RuntimeError("VTK failed")
                with (
                    patch.object(pipeline, "create_window", return_value=window),
                    patch.object(
                        pipeline, "create_interactor", return_value=interactor
                    ),
                ):
                    with self.assertRaisesRegex(RuntimeError, "VTK failed"):
                        pipeline.run()
                window.Finalize.assert_called_once()


if __name__ == "__main__":
    unittest.main()
