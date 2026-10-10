"""Numerical and cache contracts for the educational reel model adapters."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from src.reels import models
from src.reels.models import Model, prepare_model
from src.reels.stories import STORIES


class ReelModelAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.heat_models = {
            key: prepare_model(key, cache_dir=None)
            for key in (
                "fixed_end_heat_transfer_1d",
                "convective_end_heat_transfer_1d",
                "enhanced_heat_transfer_solver_1d",
                "enhanced_heat_transfer_solver_2d",
                "heat_convection_solver_1d",
                "heat_convection_solver_2d",
                "heat_conduction_solver_2d",
            )
        }

    def test_catalog_matches_nine_cfd_examples(self):
        expected = {
            path.stem
            for path in Path("src/09_cfd").glob("*.py")
            if path.stem != "__init__"
        }
        self.assertEqual(set(STORIES), expected)
        self.assertEqual(set(STORIES), {story.key for story in STORIES.values()})
        for story in STORIES.values():
            self.assertEqual(len(story.beats), 6)

    def test_fixed_end_adapter_conserves_generation_and_boundaries(self):
        model = self.heat_models["fixed_end_heat_transfer_1d"]
        temperature = model.fields[0, 0]
        self.assertAlmostEqual(temperature[0], 100.0, places=8)
        self.assertAlmostEqual(temperature[-1], 200.0, places=8)
        self.assertGreater(temperature.max(), 200.0)
        self.assertAlmostEqual(
            float(model.metrics[1][1].split()[0]), 100000.0, delta=1.0
        )

    def test_robin_adapter_matches_surface_and_analytic_profile(self):
        model = self.heat_models["convective_end_heat_transfer_1d"]
        numerical, analytical = model.chart["series"]
        self.assertAlmostEqual(numerical[0], 140.0, places=8)
        self.assertAlmostEqual(numerical[-1], 220.0, places=8)
        np.testing.assert_allclose(numerical, analytical, atol=1e-8, rtol=1e-10)

    def test_channel_adapter_uses_the_tested_imposed_flux_case(self):
        model = self.heat_models["heat_conduction_solver_2d"]
        self.assertEqual(model.params["length"], 30.0)
        self.assertEqual(model.params["hf"], 1.0)
        self.assertAlmostEqual(float(model.metrics[0][1]), 140 / 17, delta=0.02)

    def test_heat_adapters_return_consistent_finite_fields_and_charts(self):
        expected_units = {
            "fixed_end_heat_transfer_1d": "°C",
            "convective_end_heat_transfer_1d": "°C",
            "enhanced_heat_transfer_solver_1d": "°C",
            "enhanced_heat_transfer_solver_2d": "°C",
            "heat_convection_solver_1d": "scalar",
            "heat_convection_solver_2d": "scalar",
            "heat_conduction_solver_2d": "°C",
        }
        for key, model in self.heat_models.items():
            with self.subTest(key=key):
                self.assertEqual(model.unit, expected_units[key])
                self.assertEqual(model.fields.shape[0], 1)
                self.assertEqual(model.fields.shape[1:], (model.y.size, model.x.size))
                self.assertTrue(np.isfinite(model.fields).all())
                self.assertTrue(np.isfinite(model.x).all())
                self.assertTrue(np.isfinite(model.y).all())
                self.assertIn("series", model.chart)
                self.assertTrue(all(np.isfinite(series).all() for series in model.chart["series"]))

    def test_flow_adapters_preserve_grid_shapes_units_and_solid_mask(self):
        cavity = models.example_module("fluid_flow_simulator")
        x = np.linspace(0.0, 1.0, 7)
        y = np.linspace(0.0, 1.0, 5)
        u = np.full((y.size, x.size), 0.2)
        v = np.full_like(u, -0.1)
        pressure = np.zeros_like(u)
        with (
            patch.object(cavity, "run_simulation", return_value=(x, y, u, v, pressure)),
            patch.object(cavity, "compute_divergence", return_value=0.0),
        ):
            cavity_model = prepare_model("fluid_flow_simulator", cache_dir=None)

        self.assertEqual(cavity_model.unit, "m/s")
        self.assertEqual(cavity_model.fields.shape, (1, y.size, x.size))
        self.assertEqual(cavity_model.vectors.shape, (1, 2, y.size, x.size))
        self.assertTrue(np.isfinite(cavity_model.fields).all())
        self.assertTrue(np.isfinite(cavity_model.vectors).all())

        obstacle = models.example_module("obstacle_flow_simulation")
        # Match the adapter's fast-grid contract while replacing the expensive
        # time integration with finite snapshots of the same expected shape.
        nx, ny, max_iter, interval = 180, 76, 1600, 20
        snapshot_count = max_iter // interval
        velocity_snapshots = np.zeros((snapshot_count, 2, ny, nx), dtype=float)
        velocity_snapshots[:, 0] = 0.05
        solid = np.zeros((ny, nx), dtype=bool)
        solid[ny // 2 - 2 : ny // 2 + 2, nx // 4 - 3 : nx // 4 + 3] = True
        with patch.object(
            obstacle,
            "simulate_flow",
            return_value=(velocity_snapshots, solid, np.empty((0,)), np.empty((0,))),
        ) as simulation:
            obstacle_model = prepare_model("obstacle_flow_simulation", cache_dir=None)
            callback = simulation.call_args.kwargs["progress_callback"]
            callback(1, 2)

        self.assertEqual(obstacle_model.unit, "lattice units")
        self.assertEqual(obstacle_model.fields.shape, (snapshot_count, ny, nx))
        self.assertEqual(obstacle_model.vectors.shape, (snapshot_count, 2, ny, nx))
        np.testing.assert_array_equal(obstacle_model.solid, solid)
        self.assertTrue(np.isnan(obstacle_model.fields[:, solid]).all())
        self.assertTrue(np.isfinite(obstacle_model.fields[:, ~solid]).all())
        self.assertTrue(np.isfinite(obstacle_model.vectors).all())
        np.testing.assert_array_equal(
            obstacle_model.steps, np.arange(interval, max_iter + 1, interval)
        )


class ReelModelCacheTests(unittest.TestCase):
    @staticmethod
    def fake_model(key="fixed_end_heat_transfer_1d"):
        x = np.array([0.0, 0.5, 1.0])
        field = np.array([[[100.0, 150.0, 200.0]]])
        return Model(
            key=key,
            kind="rod",
            x=x,
            y=np.array([0.0]),
            fields=field,
            unit="°C",
            field_name="TEMPERATURE",
            params={"k": 50.0, "count": 3},
            note="steady test field",
            metrics=[("Peak", "200 °C")],
            chart={
                "x": x.copy(),
                "series": [field[0, 0].copy(), np.array([100.0, 145.0, 200.0])],
                "labels": ["Numerical", "Analytical"],
                "xlabel": "x (m)",
                "ylabel": "T (°C)",
            },
        )

    def test_cache_round_trip_preserves_arrays_chart_and_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(models, "source_fingerprint", return_value="a" * 64):
                with patch.object(models, "_solve", return_value=self.fake_model()) as solve:
                    first = prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)
                    second = prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)

        solve.assert_called_once()
        self.assertFalse(first.cache_hit)
        self.assertTrue(second.cache_hit)
        self.assertEqual(second.params, first.params)
        np.testing.assert_array_equal(second.fields, first.fields)
        np.testing.assert_array_equal(second.chart["x"], first.chart["x"])
        for before, after in zip(first.chart["series"], second.chart["series"]):
            np.testing.assert_array_equal(after, before)
        self.assertEqual(second.chart["labels"], first.chart["labels"])

    def test_source_fingerprint_change_uses_a_fresh_cache_entry(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(
                models, "source_fingerprint", side_effect=["a" * 64, "b" * 64]
            ):
                with patch.object(models, "_solve", side_effect=[self.fake_model(), self.fake_model()]) as solve:
                    first = prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)
                    second = prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)

        self.assertFalse(first.cache_hit)
        self.assertFalse(second.cache_hit)
        self.assertEqual(solve.call_count, 2)

    def test_corrupt_cache_is_recomputed(self):
        with tempfile.TemporaryDirectory() as directory:
            fingerprint = "c" * 64
            cache_path = Path(directory) / f"fixed_end_heat_transfer_1d-fast-{fingerprint[:16]}.npz"
            with patch.object(models, "source_fingerprint", return_value=fingerprint):
                with patch.object(models, "_solve", side_effect=[self.fake_model(), self.fake_model()]) as solve:
                    prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)
                    cache_path.write_bytes(b"not a valid numpy archive")
                    recomputed = prepare_model("fixed_end_heat_transfer_1d", cache_dir=directory)

        self.assertEqual(solve.call_count, 2)
        self.assertFalse(recomputed.cache_hit)
        self.assertAlmostEqual(recomputed.fields[0, 0, -1], 200.0)


if __name__ == "__main__":
    unittest.main()
