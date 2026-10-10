"""Numerical regression checks for the two CFD examples (no GUI required)."""

import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def _load_example(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "src/09_cfd" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cavity = _load_example("cavity_flow_example", "fluid_flow_simulator.py")
lbm = _load_example("obstacle_lbm_example", "obstacle_flow_simulation.py")


class CavityFlowTests(unittest.TestCase):
    def test_manufactured_pressure_gradient_and_checkerboard_mode(self):
        nx, ny = 8, 6
        dx, dy = 0.2, 0.3
        x = np.arange(nx) * dx
        y = np.arange(ny) * dy
        X, Y = np.meshgrid(x, y)
        u = np.zeros_like(X)
        v = np.zeros_like(X)
        rhs_u, rhs_v = cavity.momentum_rhs(u, v, 3.0 * X - 2.0 * Y, dx, dy, 0.1)
        np.testing.assert_allclose(rhs_u, -3.0)
        np.testing.assert_allclose(rhs_v, 2.0)

        # Centered collocated gradients annihilate this pressure mode; this is
        # the documented reason the pressure field is not benchmark grade.
        checkerboard = (-1.0) ** np.indices((ny, nx)).sum(axis=0)
        rhs_u, rhs_v = cavity.momentum_rhs(u, v, checkerboard, dx, dy, 0.1)
        np.testing.assert_allclose(rhs_u, 0.0)
        np.testing.assert_allclose(rhs_v, 0.0)

    def test_boundary_conditions_and_timestep_use_field_speeds(self):
        u = np.ones((4, 5))
        v = np.ones((4, 5))
        cavity.apply_boundary_conditions(u, v, -2.0)
        self.assertTrue(np.all(u[:, 0] == 0) and np.all(u[:, -1] == 0))
        self.assertTrue(np.all(u[-1, 1:-1] == -2.0))
        self.assertTrue(np.all(v[-1, :] == 0))
        zeros = np.zeros((3, 4))
        dt_slow = cavity.compute_stable_timestep(0.1, 0.2, 0.01, 1.0, zeros, zeros)
        fast_u = np.full((3, 4), 4.0)
        dt_fast = cavity.compute_stable_timestep(0.1, 0.2, 0.01, 1.0, fast_u, zeros)
        self.assertLess(dt_fast, dt_slow)

    def test_rectangular_vtk_mapping_and_centerline_profile_coordinates(self):
        x = np.array([0.0, 2.0, 4.0])
        y = np.array([0.0, 2.0, 4.0, 6.0])
        field = np.arange(12.0).reshape(4, 3)
        points, scalars = cavity.create_points(x, y, field)
        mesh = cavity.create_poly_data(3, 4, points, scalars)
        self.assertEqual(mesh.GetNumberOfPoints(), 12)
        self.assertEqual(mesh.GetNumberOfCells(), 6)
        self.assertEqual(mesh.GetPoint(5), (4.0, 2.0, 0.0))
        self.assertEqual(mesh.GetCell(0).GetPointIds().GetId(0), 0)
        self.assertEqual(mesh.GetCell(0).GetPointIds().GetId(1), 1)
        self.assertEqual(mesh.GetCell(0).GetPointIds().GetId(2), 4)
        self.assertEqual(mesh.GetCell(0).GetPointIds().GetId(3), 3)

        u = np.tile(np.arange(4.0)[:, None], (1, 3))
        v = np.tile(np.arange(3.0)[None, :], (4, 1))
        p = np.zeros_like(u)
        with patch.object(cavity.plt, "show"):
            cavity.plot_results(x, y, u, v, p)
        axes = cavity.plt.gcf().axes
        profile_lines = axes[3].lines
        np.testing.assert_array_equal(profile_lines[0].get_xdata(), u[:, 1])
        np.testing.assert_array_equal(profile_lines[0].get_ydata(), y)
        np.testing.assert_array_equal(profile_lines[1].get_xdata(), v[2, :])
        np.testing.assert_array_equal(profile_lines[1].get_ydata(), x)
        cavity.plt.close("all")

    def test_low_re_run_converges_and_exhaustion_is_explicit(self):
        x, y, u, v, p = cavity.run_simulation(
            nx=15, ny=15, reynolds=1.0, u_lid=0.1, max_iter=3000,
            tolerance=1e-3, print_interval=3000,
        )
        dx, dy = x[1] - x[0], y[1] - y[0]
        self.assertLessEqual(cavity.compute_divergence(u, v, dx, dy), 1e-3)
        self.assertLessEqual(cavity.compute_momentum_residual(u, v, p, dx, dy, 0.1), 1e-3)
        with self.assertRaises(cavity.ConvergenceError):
            cavity.run_simulation(
                nx=9, ny=9, reynolds=1.0, u_lid=0.1, max_iter=1, print_interval=1
            )

    def test_zero_lid_requires_viscosity_and_returns_rest_state(self):
        with self.assertRaisesRegex(ValueError, "nu must be supplied"):
            cavity.run_simulation(nx=5, ny=5, u_lid=0, max_iter=1)
        _, _, u, v, _ = cavity.run_simulation(
            nx=5, ny=5, u_lid=0, nu=0.1, max_iter=2, print_interval=1
        )
        self.assertEqual(np.count_nonzero(u), 0)
        self.assertEqual(np.count_nonzero(v), 0)


class LatticeBoltzmannTests(unittest.TestCase):
    def test_equilibrium_moments_for_d2q9_and_d1q3(self):
        rho = np.full((4, 5), 1.2)
        velocity = np.zeros((2, 4, 5))
        velocity[0] = 0.03
        velocity[1] = -0.02
        populations = lbm.equilibrium(rho, velocity)
        np.testing.assert_allclose(populations.sum(axis=0), rho, atol=1e-14)
        np.testing.assert_allclose(
            np.einsum("qd,qxy->dxy", lbm.C, populations), rho * velocity, atol=1e-14
        )

        c = np.array([[-1], [0], [1]])
        weights = np.array([1 / 6, 2 / 3, 1 / 6])
        rho_1d = np.array([0.8, 1.1, 1.3])
        u_1d = np.array([[0.1, -0.04, 0.0]])
        populations_1d = lbm.equilibrium(rho_1d, u_1d, c, weights)
        np.testing.assert_allclose(populations_1d.sum(axis=0), rho_1d)
        np.testing.assert_allclose(
            np.einsum("qd,qx->dx", c, populations_1d), rho_1d * u_1d
        )

    def test_zou_he_inlet_preserves_velocity_with_neq_north_south(self):
        rho0 = np.ones((8, 6))
        u0 = np.zeros((2, 8, 6))
        u0[0] = 0.03
        u0[1] = -0.01
        fin = lbm.equilibrium(rho0, u0)
        fin[2, 0, :] += 0.02  # north
        fin[4, 0, :] -= 0.005  # south: deliberate non-equilibrium imbalance
        prescribed = np.tile(np.array([[0.04], [0.02]]), (1, 6))
        lbm.apply_zou_he_velocity_inlet(fin, prescribed)
        rho, velocity = lbm.compute_macroscopic(fin[:, 0, :])
        self.assertTrue(np.all(rho > 0))
        np.testing.assert_allclose(velocity[0], 0.04, atol=1e-14)
        np.testing.assert_allclose(velocity[1], 0.02, atol=1e-14)

    def test_streaming_has_periodic_y_but_no_x_wrap_and_link_bounceback(self):
        fout = np.zeros((9, 5, 4))
        fout[1, 0, 1] = 0.25  # eastward at left edge
        fout[3, 0, 2] = 0.5  # westward at left edge must leave, not wrap
        fout[2, 2, 3] = 0.75  # north wraps periodically in y
        streamed = lbm.stream_distributions(fout)
        self.assertEqual(streamed[1, 1, 1], 0.25)
        self.assertEqual(streamed[3, -1, 2], 0.0)
        self.assertEqual(streamed[2, 2, 0], 0.75)

        obstacle = np.zeros((5, 4), dtype=bool)
        obstacle[2, 1] = True
        post = np.zeros_like(fout)
        post[1, 1, 1] = 0.6  # eastward fluid link hits the solid at (2,1)
        bounced = lbm.apply_link_bounce_back(post, lbm.stream_distributions(post), obstacle)
        self.assertEqual(bounced[3, 1, 1], 0.6)

    def test_small_runtime_snapshots_and_vtk_point_mapping(self):
        snapshots, obstacle_mask, fin, velocity = lbm.simulate_flow(
            nx=36, ny=20, max_iter=8, reynolds=10, inlet_speed=0.04,
            obstacle_length=6, obstacle_height=4, snapshot_interval=3,
        )
        self.assertEqual(len(snapshots), 3)  # steps 3, 6, and final step 8
        self.assertTrue(all(snapshot.shape == (2, 20, 36) for snapshot in snapshots))
        self.assertTrue(np.isfinite(fin).all() and np.isfinite(velocity).all())
        self.assertTrue(np.all(velocity[:, obstacle_mask] == 0.0))

        mask = np.zeros((4, 5), dtype=bool)
        mask[1, 2] = True
        field = np.ones((2, 4, 5))
        grid = lbm.create_structured_grid(5, 4, field, mask)
        self.assertEqual(grid.GetNumberOfPoints(), 20)
        self.assertEqual(grid.GetNumberOfCells(), 12)
        self.assertEqual(grid.GetPoint(7), (2.0, 1.0, 0.0))
        self.assertTrue(np.isnan(grid.GetPointData().GetScalars().GetTuple1(7)))
        self.assertEqual(grid.GetPointData().GetVectors().GetTuple3(7), (0.0, 0.0, 0.0))
        self.assertEqual(grid.GetPointData().GetScalars().GetName(), "Velocity Magnitude (lattice units)")
        self.assertEqual(grid.GetPointData().GetVectors().GetName(), "Velocity (lattice units)")

    def test_left_inlet_rejects_nonpositive_flow_direction(self):
        with self.assertRaisesRegex(ValueError, "inlet_speed must be finite, positive"):
            lbm.simulate_flow(
                nx=20, ny=12, max_iter=2, inlet_speed=-0.04,
                obstacle_length=4, obstacle_height=3,
            )


if __name__ == "__main__":
    unittest.main()
