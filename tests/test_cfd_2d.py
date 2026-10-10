"""Analytical, conservation, and field-to-geometry checks for the 2D examples."""

import importlib
import unittest

import numpy as np

from src.common.heat_transfer import (
    channel_velocity,
    channel_weights,
    solve_channel,
    solve_quarter_conduction,
    solve_transport,
)
from src.common.numerics import ConvergenceError

channel = importlib.import_module("src.09_cfd.heat_conduction_solver_2d")
conduction = importlib.import_module("src.09_cfd.enhanced_heat_transfer_solver_2d")
transport = importlib.import_module("src.09_cfd.heat_convection_solver_2d")


class TransportTests(unittest.TestCase):
    def test_harmonic_quadratic_on_rectangular_mesh(self):
        # x²-y² is harmonic, and central differences represent it exactly.
        for nx, ny in ((12, 6), (6, 12), (20, 10)):
            x, y, field, *_ = solve_transport(
                nx,
                ny,
                2,
                1,
                3,
                1,
                0,
                0,
                lambda y: -(y**2),
                lambda y: 4 - y**2,
                lambda x: x**2,
                lambda x: x**2 - 1,
                tolerance=1e-11,
            )
            np.testing.assert_allclose(
                field, x[None, :] ** 2 - y[:, None] ** 2, atol=2e-9
            )

    def test_affine_solution_with_opposing_velocity_components(self):
        # u*phi_x+v*phi_y = 2*1 - 1*2 = 0, as does the Laplacian.
        x, y, field, *_ = solve_transport(
            18,
            9,
            2,
            1,
            1,
            1,
            2,
            -1,
            lambda y: 2 * y,
            lambda y: 2 + 2 * y,
            lambda x: x,
            lambda x: x + 2,
            tolerance=1e-11,
        )
        np.testing.assert_allclose(field, x[None, :] + 2 * y[:, None], atol=2e-9)

    def test_advection_diffusion_refines_to_exponential_solution(self):
        errors = []

        def exact(x):
            return np.expm1(3 * x) / np.expm1(3)

        for nx in (12, 24, 48):
            x, _, field, *_ = solve_transport(
                nx,
                6,
                1,
                1,
                1,
                1,
                3,
                0,
                0,
                1,
                exact,
                exact,
                tolerance=1e-10,
            )
            self.assertGreaterEqual(field.min(), 0)
            self.assertLessEqual(field.max(), 1)
            errors.append(np.max(np.abs(field - exact(x)[None, :])))
        self.assertLess(errors[1], errors[0] * 0.6)
        self.assertLess(errors[2], errors[1] * 0.6)

    def test_symmetric_diffusion_center_is_independent_of_cell_aspect_ratio(self):
        for nx, ny in ((20, 20), (20, 10), (40, 20)):
            _, _, field, *_ = solve_transport(
                nx,
                ny,
                1,
                1,
                1,
                1,
                0,
                0,
                0,
                1,
                0,
                0,
                tolerance=1e-10,
            )
            self.assertAlmostEqual(field[ny // 2, nx // 2], 0.25, delta=8e-4)

    def test_invalid_inputs_and_iteration_limit_are_explicit(self):
        args = [12, 6, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0]
        with self.assertRaises(ConvergenceError) as caught:
            solve_transport(*args, tolerance=1e-12, max_iter=1)
        self.assertEqual(caught.exception.iterations, 1)
        for index, value in (
            (0, 1),
            (1, 2.5),
            (2, -1),
            (4, 0),
            (6, np.nan),
            (8, [1, 2]),
        ):
            invalid = args.copy()
            invalid[index] = value
            with self.assertRaises(ValueError):
                solve_transport(*invalid)


class QuarterConductionTests(unittest.TestCase):
    def test_zero_source_preserves_constant_including_insulated_corners(self):
        _, _, field, iterations = solve_quarter_conduction(8, 1, 3, 0, 100, 1e-10)
        np.testing.assert_array_equal(field, 100)
        self.assertEqual(iterations, 0)

    def test_second_order_refinement_against_analytical_cosine_series(self):
        # Separation of variables: cos((m+1/2)πx) obeys insulation at x=0
        # and fixed T=0 at x=1. Expand the uniform source in this basis.
        m = np.arange(180)
        eigen = (m + 0.5) * np.pi
        signs = (-1.0) ** m
        reference = np.sum(
            4
            * signs[:, None]
            * signs[None, :]
            / (
                eigen[:, None]
                * eigen[None, :]
                * (eigen[:, None] ** 2 + eigen[None, :] ** 2)
            )
        )
        errors = []
        for n in (8, 16, 32):
            _, _, field, _ = solve_quarter_conduction(n, 1, 1, 1, 0, 1e-11)
            self.assertEqual(field.argmax(), 0)
            self.assertGreaterEqual(field.min(), 0)
            np.testing.assert_allclose(field, field.T, atol=1e-12)
            errors.append(abs(field[0, 0] - reference))
        self.assertLess(errors[1], errors[0] * 0.27)
        self.assertLess(errors[2], errors[1] * 0.27)

    def test_generation_balances_flux_into_both_fixed_edges(self):
        n, L, k, source = 16, 1.3, 2.7, 4.1
        _, _, field, _ = solve_quarter_conduction(n, L, k, source, 20, 1e-11)
        dx = L / n
        widths = np.full(n, dx)
        widths[0] /= 2
        outward = k / dx * np.dot(widths, field[:-1, -2] - field[:-1, -1])
        outward += k / dx * np.dot(widths, field[-2, :-1] - field[-1, :-1])
        self.assertAlmostEqual(outward, source * (L - dx / 2) ** 2, delta=1e-7)


class ChannelTests(unittest.TestCase):
    def test_grid_roundoff_and_nondivisible_spacing(self):
        nx, ny, x, y = channel.initialize_grid(0.3, 0.3, 0.1, 0.1)
        self.assertEqual((nx, ny, len(x), len(y)), (3, 3, 4, 4))
        with self.assertRaises(ValueError):
            channel.initialize_grid(1, 1, 0.3, 0.1)

    def test_velocity_is_forward_no_slip_with_unit_mean(self):
        y = np.linspace(0, 2, 21)
        u = channel_velocity(y, 2)
        np.testing.assert_allclose(u[[0, 10, -1]], [0, 1.5, 0])
        self.assertGreaterEqual(u.min(), 0)
        weights, heights = channel_weights(y, 2)
        self.assertAlmostEqual(weights.sum(), 2)
        self.assertAlmostEqual(heights.sum(), 2)

    def test_insulated_channel_and_legacy_wall_relaxation_preserve_constant(self):
        initial = np.full((7, 13), 100.0)
        field, iterations, _ = solve_channel(
            initial,
            0.2,
            1 / 6,
            np.linspace(0, 1, 7),
            100,
            0,
            1,
            10,
            1,
            1,
            1e-10,
        )
        np.testing.assert_array_equal(field, initial)
        self.assertEqual(iterations, 0)
        channel.apply_boundary_conditions(initial, 100, 12, 6, 1, 0.2, 1 / 6, 0)
        np.testing.assert_allclose(initial, 100, atol=1e-13)

    def test_wall_heating_conserves_energy_and_matches_developed_nusselt(self):
        # Dimensionless h=k=q=1, rho*c=10, L=30: enough thermal development.
        y = np.linspace(0, 1, 21)
        initial = np.full((21, 101), 20.0)
        field, _, _ = solve_channel(initial, 0.3, 0.05, y, 20, 1, 1, 10, 1, 1, 1e-10)
        np.testing.assert_array_equal(initial, 20)  # kernel must not mutate callers
        np.testing.assert_allclose(field, field[::-1, :], atol=1e-10)
        self.assertGreaterEqual(field.min(), 20)
        weights, heights = channel_weights(y, 1)
        outlet_advected = 10 * np.dot(weights, field[:, -1] - 20)
        inlet_conduction = np.dot(heights, field[:, 1] - 20) / 0.3
        self.assertAlmostEqual(
            outlet_advected + inlet_conduction, 2 * (30 - 0.3 / 2), delta=1e-6
        )
        bulk = channel.calculate_bulk_temperature(field, y, 0.05, 1, 100)
        _, nu = channel.calculate_ht_and_nusselt(bulk, np.ones(101), 1, 1, field)
        self.assertTrue(np.isnan(nu[0]))
        # Analytical fully developed equal-flux parallel plates: Nu=140/17.
        self.assertAlmostEqual(nu[-10], 140 / 17, delta=0.02)
        self.assertAlmostEqual(
            channel.calculate_heat_flux(field, 1, 0.05)[-10], 1, delta=0.015
        )

    def test_cooling_sign_and_explicit_nonconvergence(self):
        y = np.linspace(0, 1, 7)
        initial = np.full((7, 13), 20.0)
        args = (initial, 0.2, 1 / 6, y, 20, -1, 1, 10, 1, 1)
        field, _, _ = solve_channel(*args, 1e-10)
        self.assertLessEqual(field.max(), 20)
        self.assertLess(field.min(), 20)
        with self.assertRaises(ConvergenceError):
            solve_channel(*args, 1e-12, max_iter=1)


class FieldVisualizationTests(unittest.TestCase):
    def test_rectangular_scalar_coordinates_heights_and_cell_orientation(self):
        x, y = np.array([2.0, 3.0, 5.0, 8.0]), np.array([-1.0, 0.0, 4.0])
        field = np.arange(12.0).reshape(3, 4)
        for example in (channel, conduction):
            points, scalars = example.create_points(x, y, field)
            mesh = example.create_poly_data(4, 3, points, scalars)
            self.check_mesh(mesh, x, y, field, height=0)
        mesh = transport.create_2d_surface(x, y, field, height_scale=2)
        self.check_mesh(mesh, x, y, field, height=2)

    def check_mesh(self, mesh, x, y, field, height):
        self.assertEqual(mesh.GetNumberOfPoints(), field.size)
        self.assertEqual(mesh.GetNumberOfCells(), (len(x) - 1) * (len(y) - 1))
        for j, py in enumerate(y):
            for i, px in enumerate(x):
                index = j * len(x) + i
                actual = mesh.GetPoint(index)
                np.testing.assert_allclose(actual, [px, py, field[j, i] / 11 * height])
                self.assertEqual(
                    mesh.GetPointData().GetScalars().GetValue(index), field[j, i]
                )
        cell = mesh.GetCell(0)
        a, b, c = (np.array(mesh.GetPoint(cell.GetPointId(i))) for i in range(3))
        self.assertGreater(np.cross(b - a, c - a)[2], 0)

    def test_uniform_surface_is_flat_and_malformed_fields_rejected(self):
        mesh = transport.create_2d_surface([0, 1], [0, 1], np.full((2, 2), 100.0))
        self.assertEqual(mesh.GetBounds()[4:], (0, 0))
        for example in (channel, conduction):
            with self.assertRaises(ValueError):
                example.create_points([0, 1, 2], [0, 1], np.ones((3, 2)))
        with self.assertRaises(ValueError):
            transport.create_2d_surface([0, 1], [0, 1], np.array([[0, np.nan], [1, 2]]))

    def test_color_legend_preserves_actual_field_range(self):
        field = np.array([[30.0, 60.0], [90.0, 120.0]])
        for example in (channel, conduction, transport):
            lookup = example.create_color_map(field)
            self.assertEqual(lookup.GetRange(), (30, 120))
            self.assertGreater(lookup.GetTableValue(0)[2], lookup.GetTableValue(0)[0])
            self.assertGreater(
                lookup.GetTableValue(255)[0], lookup.GetTableValue(255)[2]
            )


if __name__ == "__main__":
    unittest.main()
