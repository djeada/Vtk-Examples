"""Stability, exactness, and conservation checks for the beginner CFD examples."""

import importlib
import unittest

import numpy as np

diffusion = importlib.import_module("src.09_cfd.transient_heat_diffusion_2d")
advection = importlib.import_module("src.09_cfd.rotating_blob_advection_2d")


class HeatDiffusionTests(unittest.TestCase):
    def test_stable_step_sits_at_the_fourier_limit(self):
        dt = diffusion.stable_time_step(2.0, 0.1, 0.2, safety=1.0)
        self.assertAlmostEqual(2.0 * dt / 0.1**2 + 2.0 * dt / 0.2**2, 0.5)
        with self.assertRaises(ValueError):
            diffusion.stable_time_step(0, 0.1, 0.1)

    def test_sine_mode_decays_by_the_exact_discrete_factor(self):
        n, length, alpha = 20, 1.0, 0.3
        x = np.linspace(0, length, n + 1)
        dx = length / n
        dt = diffusion.stable_time_step(alpha, dx, dx, safety=0.8)
        mode = np.outer(np.sin(np.pi * x), np.sin(np.pi * x))
        fo = alpha * dt / dx**2
        factor = 1 - 8 * fo * np.sin(np.pi * dx / 2) ** 2
        stepped = diffusion.ftcs_step(mode, alpha, dt, dx, dx)
        np.testing.assert_allclose(stepped, factor * mode, atol=1e-14)

    def test_fine_grid_follows_the_continuous_decay_rate(self):
        n, alpha, steps = 64, 1.0, 400
        x = np.linspace(0, 1, n + 1)
        dx = 1 / n
        dt = diffusion.stable_time_step(alpha, dx, dx)
        field = np.outer(np.sin(np.pi * x), np.sin(np.pi * x))
        for _ in range(steps):
            field = diffusion.ftcs_step(field, alpha, dt, dx, dx)
        expected = np.exp(-2 * np.pi**2 * alpha * steps * dt)
        self.assertAlmostEqual(field.max() / expected, 1, delta=2e-3)

    def test_exceeding_the_limit_blows_up(self):
        n = 16
        dx = 1 / n
        rng = np.random.default_rng(0)
        field = np.zeros((n + 1, n + 1))
        field[1:-1, 1:-1] = rng.uniform(-1, 1, (n - 1, n - 1))
        dt = diffusion.stable_time_step(1.0, dx, dx, safety=1.5)
        for _ in range(200):
            field = diffusion.ftcs_step(field, 1.0, dt, dx, dx)
        self.assertGreater(np.abs(field).max(), 1e3)

    def test_simulation_obeys_the_maximum_principle(self):
        x, y, times, steps, frames, dt = diffusion.simulate(n=32, duration=20, snapshots=11)
        self.assertEqual(frames.shape, (len(times), 33, 33))
        self.assertTrue(np.all(np.diff(times) > 0))
        np.testing.assert_allclose(times, steps * dt)
        peaks = frames.max(axis=(1, 2))
        self.assertTrue(np.all(np.diff(peaks) <= 1e-12))
        self.assertGreaterEqual(frames.min(), 20.0 - 1e-12)
        np.testing.assert_allclose(frames[:, 0, :], 20.0)
        np.testing.assert_allclose(frames[:, :, -1], 20.0)


class BlobAdvectionTests(unittest.TestCase):
    def test_cfl_step_hits_the_requested_courant_number(self):
        x = np.linspace(0, 1, 41)
        u, v = advection.rotation_velocity(x, x)
        dt = advection.cfl_time_step(u, v, 1 / 40, 1 / 40, cfl=0.7)
        self.assertAlmostEqual(np.max(np.abs(u) * 40 + np.abs(v) * 40) * dt, 0.7)

    def test_exact_solution_returns_after_one_turn(self):
        x = np.linspace(0, 1, 51)
        np.testing.assert_allclose(
            advection.exact_blob(x, x, 1.0), advection.gaussian_blob(x, x, (0.5, 0.7), 0.07), atol=1e-12
        )

    def test_upwind_is_monotone_and_conserves_dye(self):
        x, y, times, steps, frames, dt, _ = advection.simulate(n=60, revolutions=0.25, snapshots=5)
        self.assertGreaterEqual(frames.min(), 0)
        self.assertLessEqual(frames.max(), frames[0].max() + 1e-12)
        self.assertAlmostEqual(times[-1], 0.25)
        # Flux form: interior faces cancel, so dye is conserved until the
        # smeared tails reach an outflow edge.
        n = 60
        x = np.linspace(0, 1, n + 1)
        u, v = advection.rotation_velocity(x, x)
        dt = advection.cfl_time_step(u, v, 1 / n, 1 / n)
        phi = advection.gaussian_blob(x, x)
        for _ in range(10):
            phi = advection.upwind_step(phi, u, v, dt, 1 / n, 1 / n)
        self.assertAlmostEqual(phi.sum() / frames[0].sum(), 1, delta=1e-8)

    def test_numerical_diffusion_shrinks_with_refinement(self):
        errors = []
        for n in (40, 80):
            x, y, times, _, frames, _, _ = advection.simulate(n=n, revolutions=0.25, snapshots=2)
            errors.append(np.abs(frames[-1] - advection.exact_blob(x, y, times[-1])).max())
        self.assertLess(errors[1], 0.75 * errors[0])

    def test_peak_drops_even_though_the_exact_peak_stays_one(self):
        x, y, times, _, frames, _, _ = advection.simulate(n=60, snapshots=3)
        self.assertLess(frames[-1].max(), 0.5)
        self.assertAlmostEqual(advection.exact_blob(x, y, times[-1]).max(), 1, delta=1e-3)


if __name__ == "__main__":
    unittest.main()
