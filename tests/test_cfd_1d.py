"""Numerical checks for the shared 1D heat and transport examples."""

import importlib
import unittest

import numpy as np

from src.common.numerics import (
    ConvergenceError,
    tdma_solve,
    validate_count,
    validate_iteration_controls,
    validate_positive,
)

fixed = importlib.import_module("src.09_cfd.fixed_end_heat_transfer_1d")
robin = importlib.import_module("src.09_cfd.convective_end_heat_transfer_1d")
nonlinear = importlib.import_module("src.09_cfd.enhanced_heat_transfer_solver_1d")
convection = importlib.import_module("src.09_cfd.heat_convection_solver_1d")


class SharedNumericsTests(unittest.TestCase):
    def test_tdma_matches_independent_dense_solve_without_mutation(self):
        main = np.array([4.0, 5.0, 6.0, 7.0])
        lower = np.array([0.0, -1.0, -2.0, -1.5])
        upper = np.array([-0.5, -1.0, -1.25, 0.0])
        rhs = np.array([1.0, 2.0, -3.0, 4.0])
        originals = [array.copy() for array in (main, lower, upper, rhs)]
        dense = np.diag(main) + np.diag(lower[1:], -1) + np.diag(upper[:-1], 1)

        result = tdma_solve(main, lower, upper, rhs)

        np.testing.assert_allclose(result, np.linalg.solve(dense, rhs))
        for actual, expected in zip((main, lower, upper, rhs), originals):
            np.testing.assert_array_equal(actual, expected)

    def test_tdma_validates_vectors_and_pivots(self):
        with self.assertRaisesRegex(ValueError, "same length"):
            tdma_solve([1, 1], [0], [0, 0], [1, 1])
        with self.assertRaisesRegex(ValueError, "one-dimensional"):
            tdma_solve([[1]], [0], [0], [1])
        with self.assertRaisesRegex(ValueError, "finite"):
            tdma_solve([np.nan], [0], [0], [1])
        with self.assertRaisesRegex(np.linalg.LinAlgError, "pivot"):
            tdma_solve([0], [0], [0], [1])

    def test_validators_reject_invalid_controls(self):
        for value in (0, -1, np.inf, np.nan):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_positive(value=value)
        for value in (True, 0, 1.5):
            with self.subTest(value=value), self.assertRaises(ValueError):
                validate_count("nodes", value)
        self.assertEqual(validate_count("nodes", np.int64(3)), 3)
        self.assertEqual(validate_iteration_controls(1e-8, 50), (1e-8, 50))
        with self.assertRaises(ValueError):
            validate_iteration_controls(0, 50)

    def test_convergence_error_exposes_context(self):
        error = ConvergenceError("solver stopped", 7, 1.25e-3)
        self.assertEqual(error.iterations, 7)
        self.assertEqual(error.residual, 1.25e-3)
        self.assertIn("7 iterations", str(error))
        self.assertIn("0.00125", str(error))


class LinearHeatTests(unittest.TestCase):
    def test_fixed_end_solution_matches_parabolic_analytic_solution(self):
        n, length, conductivity, source = 24, 1.7, 2.4, 3.2
        diagonals, rhs = fixed.setup_problem(
            n, length, conductivity, source, 15.0, 70.0
        )
        numerical = tdma_solve(*diagonals, rhs)
        x = np.linspace(0, length, n + 1)
        analytic = fixed.analytical_solution(
            x, length, conductivity, source, 15.0, 70.0
        )
        np.testing.assert_allclose(numerical, analytic, rtol=1e-12, atol=1e-12)

    def test_fixed_end_energy_balance_uses_second_order_end_derivatives(self):
        n, length, conductivity, source = 16, 1.3, 4.0, 8.5
        diagonals, rhs = fixed.setup_problem(
            n, length, conductivity, source, 10.0, 32.0
        )
        temperature = tdma_solve(*diagonals, rhs)
        dx = length / n
        _, _, total_out = fixed.calculate_energy_balance(
            temperature, dx, conductivity
        )
        self.assertAlmostEqual(total_out, source * length, places=10)

    def test_robin_solution_matches_linear_analytic_solution(self):
        n, length, conductivity = 32, 1.2, 3.5
        h_left, h_right, ambient_left, ambient_right = 70.0, 45.0, 120.0, 25.0
        dx = length / n
        system = robin.initialize_matrices(
            n, dx, conductivity, h_left, h_right, ambient_left, ambient_right
        )
        numerical = tdma_solve(*system[:3], system[3])
        x = np.linspace(0, length, n + 1)
        analytic = robin.analytical_solution(
            x,
            length,
            conductivity,
            h_left,
            h_right,
            ambient_left,
            ambient_right,
        )
        np.testing.assert_allclose(numerical, analytic, rtol=1e-12, atol=1e-12)

    def test_vtk_rod_scalar_values_match_their_axial_rings(self):
        coordinates = np.array([0.0, 0.5, 1.0])
        scalar_profiles = (
            (fixed, np.array([10.0, 20.0, 30.0]), "Temperature"),
            (robin, np.array([11.0, 21.0, 31.0]), "Temperature"),
            (nonlinear, np.array([12.0, 22.0, 32.0]), "Temperature"),
            (convection, np.array([1.0, 2.0, 3.0]), "Phi"),
        )
        resolution = 8
        for module, scalar_profile, field_name in scalar_profiles:
            with self.subTest(module=module.__name__):
                polydata, _ = module.create_cylinder_rod(
                    coordinates, scalar_profile, resolution=resolution
                )
                scalars = polydata.GetPointData().GetScalars()
                self.assertEqual(scalars.GetName(), field_name)
                actual = np.array(
                    [
                        scalars.GetTuple1(index)
                        for index in range(scalars.GetNumberOfTuples())
                    ]
                )
                expected = np.repeat(scalar_profile, resolution)
                np.testing.assert_array_equal(actual, expected)

        with self.assertRaisesRegex(ValueError, "equal length"):
            convection.create_cylinder_rod(coordinates, np.array([0.0, 1.0]))

    def test_convection_analytic_solution_is_stable_for_extreme_and_tiny_peclet(self):
        x = np.array([0.0, 0.2, 0.5, 0.8, 1.0])
        for peclet in (-1000.0, 1000.0):
            result = convection.analytical_solution(x, 1.0, peclet, 2.0, 5.0)
            self.assertTrue(np.all(np.isfinite(result)))
            self.assertAlmostEqual(result[0], 2.0)
            self.assertAlmostEqual(result[-1], 5.0)
            self.assertTrue(np.all(np.diff(result) >= 0))

        tiny = convection.analytical_solution(x, 1.0, 1e-14, 2.0, 5.0)
        np.testing.assert_allclose(tiny, 2.0 + 3.0 * x, rtol=1e-14, atol=1e-14)
        zero = convection.analytical_solution(x, 1.0, 0.0, 2.0, 5.0)
        np.testing.assert_array_equal(zero, 2.0 + 3.0 * x)


class NonlinearHeatTests(unittest.TestCase):
    @staticmethod
    def solve(n, beta=0.0015, tolerance=1e-11):
        return nonlinear.solve_heat_transfer_equation(
            n=n,
            L=1.5,
            k_c=4.0,
            beta=beta,
            h_l=90.0,
            h_r=60.0,
            T_l=180.0,
            T_r=20.0,
            Rtol=tolerance,
        )

    def test_zero_temperature_coefficient_matches_robin_solution(self):
        x, temperature, _, profile = self.solve(30, beta=0.0)
        analytic = robin.analytical_solution(
            x, 1.5, 4.0, 90.0, 60.0, 180.0, 20.0
        )
        np.testing.assert_allclose(temperature, analytic, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(profile, 4.0)

    def test_nonlinear_face_and_boundary_fluxes_are_conserved(self):
        x, temperature, _, _ = self.solve(50)
        profile = nonlinear.conductivity(temperature, 4.0, 0.0015)
        lower = np.minimum(profile[:-1], profile[1:])
        upper = np.maximum(profile[:-1], profile[1:])
        faces = lower * (2.0 / (1.0 + lower / upper))
        face_flux = faces * (temperature[:-1] - temperature[1:]) / np.diff(x)
        left_flux = 90.0 * (180.0 - temperature[0])
        right_flux = 60.0 * (temperature[-1] - 20.0)
        np.testing.assert_allclose(face_flux, left_flux, rtol=1e-8, atol=1e-8)
        self.assertAlmostEqual(left_flux, right_flux, places=7)

    def test_nonlinear_solution_refines_toward_integrated_conductivity_solution(self):
        length, k_c, beta = 1.5, 4.0, 0.0015
        h_left, h_right, ambient_left, ambient_right = 90.0, 60.0, 180.0, 20.0

        def primitive(value):
            return k_c * (value + 0.5 * beta * value**2)

        def flux_residual(flux):
            left_surface = ambient_left - flux / h_left
            right_surface = ambient_right + flux / h_right
            return primitive(left_surface) - primitive(right_surface) - flux * length

        low, high = 0.0, 1000.0
        for _ in range(100):
            midpoint = (low + high) / 2
            if flux_residual(midpoint) > 0:
                low = midpoint
            else:
                high = midpoint
        flux = (low + high) / 2
        left_surface = ambient_left - flux / h_left
        def exact_temperature(x):
            target_primitive = primitive(left_surface) - flux * x
            return (-1.0 + np.sqrt(1.0 + 2.0 * beta * target_primitive / k_c)) / beta

        coarse_x, coarse_temperature, _, _ = self.solve(20)
        fine_x, fine_temperature, _, _ = self.solve(80)
        coarse_exact = exact_temperature(coarse_x)
        fine_exact = exact_temperature(fine_x)
        coarse_error = np.max(np.abs(coarse_temperature - coarse_exact))
        fine_error = np.max(np.abs(fine_temperature - fine_exact))
        self.assertLess(fine_error, coarse_error)

    def test_rejects_nonpositive_conductivity_and_reports_iteration_limit(self):
        with self.assertRaisesRegex(ValueError, "conductivity"):
            nonlinear.conductivity(np.array([0.0, 100.0]), 2.0, -0.02)
        with self.assertRaises(ConvergenceError) as caught:
            nonlinear.solve_heat_transfer_equation(
                n=30,
                L=1.5,
                k_c=4.0,
                beta=0.01,
                h_l=90.0,
                h_r=60.0,
                T_l=180.0,
                T_r=20.0,
                Rtol=1e-14,
                max_iterations=1,
            )
        self.assertEqual(caught.exception.iterations, 1)


if __name__ == "__main__":
    unittest.main()
