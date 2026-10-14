"""Headless geometry checks for the standalone procedural VTK scenes."""
import importlib
import math
import unittest

import numpy as np


ROOT = "src.02_advanced_shapes"


def scene(name):
    return importlib.import_module(f"{ROOT}.{name}")


class ProceduralSceneTests(unittest.TestCase):
    def test_lava_planet_has_displaced_surface_and_scalars(self):
        mesh = scene("lava_planet").make_planet()
        self.assertGreater(mesh.GetNumberOfPolys(), 100)
        self.assertEqual(mesh.GetPointData().GetScalars().GetName(), "elevation")
        lo, hi = mesh.GetPointData().GetScalars().GetRange()
        self.assertLess(lo, 0)
        self.assertGreater(hi, 0)

    def test_metaballs_field_updates_and_contour_exists(self):
        vtk = importlib.import_module("vtk")
        image, update = scene("metaballs").make_field(30)
        scalars = image.GetPointData().GetScalars()
        before = np.array([
            scalars.GetTuple1(i) for i in range(scalars.GetNumberOfTuples())
        ])
        update(1.2)
        after = np.array([
            scalars.GetTuple1(i) for i in range(scalars.GetNumberOfTuples())
        ])
        self.assertFalse(np.array_equal(before, after))
        contour = vtk.vtkFlyingEdges3D()
        contour.SetInputData(image)
        contour.SetValue(0, 0.52)
        contour.Update()
        self.assertGreater(contour.GetOutput().GetNumberOfCells(), 0)

    def test_gravity_well_sags_at_center(self):
        mod = scene("gravity_well")
        self.assertLess(mod.well_height(0, 0), mod.well_height(3, 0))
        self.assertGreater(mod.make_well().GetNumberOfCells(), 100)

    def test_tornado_has_separate_strands(self):
        mesh = scene("tornado_ribbons").make_tornado(4, 12)
        self.assertEqual(mesh.GetNumberOfLines(), 4)
        self.assertEqual(mesh.GetNumberOfPoints(), 48)

    def test_slicer_can_advance(self):
        renderer, step = scene("slice_sphere").build_scene()
        self.assertEqual(renderer.GetActors().GetNumberOfItems(), 2)
        step()
        step()

    def test_explosion_displaces_independent_triangles(self):
        mesh, update = scene("exploding_mesh").make_fragments()
        self.assertGreater(mesh.GetNumberOfPolys(), 100)
        self.assertEqual(mesh.GetNumberOfPoints(), mesh.GetNumberOfPolys() * 3)
        before = mesh.GetPoint(0)
        update(1)
        self.assertNotEqual(before, mesh.GetPoint(0))
        update(0)
        for a, b in zip(before, mesh.GetPoint(0)):
            self.assertTrue(math.isclose(a, b, abs_tol=1e-6))

    def test_terrain_is_not_flat(self):
        mesh = scene("procedural_terrain").make_terrain()
        self.assertGreater(mesh.GetNumberOfPoints(), 1000)
        bounds = mesh.GetBounds()
        self.assertGreater(bounds[5] - bounds[4], 1)

    def test_synthetic_volume_has_multiple_intensities(self):
        image = scene("synthetic_volume_scanner").make_volume(25)
        self.assertEqual(image.GetDimensions(), (25, 25, 25))
        low, high = image.GetPointData().GetScalars().GetRange()
        self.assertEqual(low, 0)
        self.assertGreater(high, 100)

    def test_implicit_boolean_has_surface(self):
        mesh = scene("boolean_geometry").make_difference(35)
        self.assertGreater(mesh.GetNumberOfCells(), 100)
        self.assertGreater(mesh.GetNumberOfPoints(), 100)

    def test_recursive_sculpture_grows_with_depth(self):
        mod = scene("fractal_sculpture")
        small = mod.make_sculpture(0)
        large = mod.make_sculpture(2)
        self.assertGreater(large.GetNumberOfCells(), small.GetNumberOfCells())

    def test_scene_builders_can_step_without_a_window(self):
        for name in (
            "lava_planet",
            "metaballs",
            "gravity_well",
            "tornado_ribbons",
            "slice_sphere",
            "exploding_mesh",
            "procedural_terrain",
            "synthetic_volume_scanner",
            "boolean_geometry",
            "fractal_sculpture",
        ):
            with self.subTest(name=name):
                renderer, advance = scene(name).build_scene()
                self.assertGreater(renderer.GetViewProps().GetNumberOfItems(), 0)
                advance()


if __name__ == "__main__":
    unittest.main()
