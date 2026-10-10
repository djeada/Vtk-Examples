import importlib
import unittest

import numpy as np


create_surface = importlib.import_module(
    "src.08_integration_with_ui.matplotlib_surface_plot"
).create_surface_from_grid


class SurfaceTests(unittest.TestCase):
    def test_rectangular_mesh_has_valid_coordinates_and_cells(self):
        x, y = np.meshgrid(np.arange(5), np.arange(3))
        z = x + y
        surface = create_surface(x, y, z)
        self.assertEqual(surface.GetNumberOfPoints(), 15)
        self.assertEqual(surface.GetNumberOfPolys(), 8)
        for index, coordinate in enumerate(zip(x.flat, y.flat, z.flat)):
            self.assertEqual(surface.GetPoint(index), coordinate)
        for index in range(surface.GetNumberOfCells()):
            cell = surface.GetCell(index)
            self.assertEqual(cell.GetNumberOfPoints(), 4)
            points = [cell.GetPoints().GetPoint(j) for j in range(4)]
            self.assertEqual(max(p[0] for p in points) - min(p[0] for p in points), 1)
            self.assertEqual(max(p[1] for p in points) - min(p[1] for p in points), 1)

    def test_incompatible_coordinate_arrays_are_rejected(self):
        for arrays in (
            (np.zeros((2, 3)), np.zeros((3, 2)), np.zeros((2, 3))),
            (np.zeros(3), np.zeros(3), np.zeros(3)),
        ):
            with self.subTest(shapes=[array.shape for array in arrays]):
                with self.assertRaisesRegex(ValueError, "same 2D shape"):
                    create_surface(*arrays)


if __name__ == "__main__":
    unittest.main()
