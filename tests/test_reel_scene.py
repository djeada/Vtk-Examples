"""Verify that reel displays preserve scientific arrays and original viewports."""

import os
from pathlib import Path
import sys
import unittest

import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

from src.reels.models import Model
from src.reels.scene import SceneRenderer, capture_example, color_samples, lookup_table


class ReelColorTests(unittest.TestCase):
    def test_color_range_covers_all_snapshots_and_ignores_solid_nan(self):
        data = np.array([[[0.1, np.nan]], [[0.5, np.nan]]])
        table = lookup_table(data, color_samples("viridis"))
        self.assertEqual(table.GetRange(), (0.1, 0.5))
        self.assertEqual(table.GetNumberOfTableValues(), 256)


@unittest.skipIf(
    sys.platform.startswith("linux") and not os.environ.get("DISPLAY"),
    "VTK render checks require a display or xvfb-run",
)
class ReelSceneTests(unittest.TestCase):
    def model(self):
        x, y = np.arange(10.0), np.arange(8.0)
        vectors = np.zeros((2, 2, 8, 10))
        vectors[0, 0] = 0.1
        vectors[1, 1] = 0.2
        solid = np.zeros((8, 10), bool)
        solid[3:5, 3:5] = True
        fields = np.hypot(vectors[:, 0], vectors[:, 1])
        fields[:, solid] = np.nan
        return Model(
            "synthetic",
            "flow",
            x,
            y,
            fields,
            "lattice units",
            "SPEED",
            {},
            "test",
            [],
            vectors=vectors,
            solid=solid,
            steps=np.array([20, 40]),
        )

    def test_playback_interpolates_vectors_before_speed_and_preserves_orientation(self):
        model = self.model()
        scene = SceneRenderer(model, 320, 240)
        try:
            image = scene.render(0.5)
            displayed = vtk_to_numpy(scene.scalars).reshape(9, 10)
            np.testing.assert_array_equal(displayed[-1], displayed[0])
            values = displayed[:-1]
            np.testing.assert_allclose(values[~model.solid], np.hypot(0.05, 0.1))
            self.assertTrue(np.isnan(values[model.solid]).all())
            self.assertEqual(image.size, (320, 240))
            self.assertEqual(scene.mesh.GetPoint(2 + 6 * 10), (2.0, 6.0, 0.0))
            self.assertEqual(scene.legend[0], (0.1, 0.2))
        finally:
            scene.close()

    def test_preview_then_restart_resets_tracers_outside_solid(self):
        model = self.model()
        scene = SceneRenderer(model, 320, 240)
        try:
            initial = scene.initial_particles.copy()
            ix, iy = (
                np.rint(initial[:, 0]).astype(int),
                np.rint(initial[:, 1]).astype(int) % len(model.y),
            )
            self.assertFalse(model.solid[iy, ix].any())
            scene.render(0)
            scene.render(0.9)
            self.assertFalse(np.array_equal(scene.particles, initial))
            scene.render(0)
            np.testing.assert_array_equal(scene.particles, initial)
            np.testing.assert_allclose(
                vtk_to_numpy(scene.scalars).reshape(9, 10)[:-1][~model.solid], 0.1
            )
        finally:
            scene.close()

    def test_periodic_tracers_interpolate_across_seam_and_use_full_row_period(self):
        model = self.model()
        scene = SceneRenderer(model, 320, 240)
        try:
            velocity = np.zeros((2, 8, 10))
            velocity[0, -1] = 4
            np.testing.assert_allclose(
                scene._velocity_at(np.array([[1.0, 7.5]]), velocity), [[2.0, 0.0]]
            )
            velocity[:] = 0
            velocity[1] = 1
            scene.particles[:] = [1.0, 7.9]
            scene._update_particles(velocity, 0.2)
            np.testing.assert_allclose(scene.particles[:, 0], 1.0)
            np.testing.assert_allclose(scene.particles[:, 1], 0.1, atol=1e-12)
        finally:
            scene.close()

    def test_generic_capture_keeps_all_views_and_restores_script_arguments(self):
        arguments = sys.argv[:]
        model, renderers = capture_example(
            Path("src/02_advanced_shapes/multiple_independent_objects.py"), 320, 240
        )
        self.assertEqual(sys.argv, arguments)
        self.assertEqual(model.params["viewports"], 3)
        self.assertEqual(model.metrics[0], ("Scene actors", "3"))
        self.assertGreater(int(model.metrics[1][1].replace(",", "")), 0)
        scene = SceneRenderer(model, 320, 240, renderers)
        try:
            self.assertEqual(scene.window.GetRenderers().GetNumberOfItems(), 3)
            self.assertIsNone(scene.legend)
            self.assertEqual(scene.render(0.5).size, (320, 240))
        finally:
            scene.close()


if __name__ == "__main__":
    unittest.main()
