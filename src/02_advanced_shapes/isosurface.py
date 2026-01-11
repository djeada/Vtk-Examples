"""
3D Temperature Isosurface Visualization with VTK

This module demonstrates isosurface extraction and visualization from a physically
meaningful CFD problem: steady-state 3D heat conduction with internal heat generation.
Isosurfaces (isotherms) are fundamental visualization tools in CFD for identifying
regions of equal temperature, pressure, or concentration.

Physical Problem:
-----------------
We solve the steady-state heat equation with uniform volumetric heat generation
in a cubic domain:

    ∇²T + Q/k = 0

    where:
    - T = temperature (°C)
    - Q = volumetric heat generation rate (W/m³)
    - k = thermal conductivity (W/(m·K))

For this example, we solve with:
- Dirichlet boundary conditions: T = 0°C on all six faces
- Uniform internal heat generation: Q = constant

The analytical solution for a cube of side L with zero boundary temperature
and uniform heat generation is approximately:

    T(x,y,z) ≈ (Q/k) * f(x,y,z)

where f(x,y,z) captures the spatial variation peaking at the center.

Numerical Method:
-----------------
The Poisson equation is discretized using second-order central finite differences:

    (T[i+1,j,k] - 2T[i,j,k] + T[i-1,j,k])/dx² +
    (T[i,j+1,k] - 2T[i,j,k] + T[i,j-1,k])/dy² +
    (T[i,j,k+1] - 2T[i,j,k] + T[i,j,k-1])/dz² = -Q/k

Rearranging for the Gauss-Seidel iteration:

    T[i,j,k] = (1/6)(T[i+1,j,k] + T[i-1,j,k] + T[i,j+1,k] +
                     T[i,j-1,k] + T[i,j,k+1] + T[i,j,k-1] + h²Q/k)

Isosurface Concepts:
--------------------
An isosurface connects all points in a 3D scalar field that share the same value.
For temperature fields, isosurfaces are called isotherms. Key properties:

1. Marching Cubes Algorithm: VTK's vtkContourFilter uses this algorithm to
   extract polygonal meshes representing isosurfaces.

2. Multiple Isovalues: Extracting several isosurfaces at different temperatures
   reveals the thermal structure of the domain.

3. Heat Flux Direction: Heat flows perpendicular to isotherms, from hot to cold.
   Closely spaced isotherms indicate high temperature gradients (high heat flux).

CFD Educational Value:
----------------------
This example demonstrates:
- Poisson equation with source term (heat generation)
- Dirichlet boundary conditions on all faces
- Finite difference discretization of elliptic PDEs
- Isosurface extraction using marching cubes
- Color mapping for scalar field visualization
- Relationship between isotherms and heat flux

References:
-----------
1. Incropera, F.P., et al. "Fundamentals of Heat and Mass Transfer"
2. Patankar, S.V. "Numerical Heat Transfer and Fluid Flow"
3. Lorensen, W.E., Cline, H.E. (1987) "Marching Cubes" SIGGRAPH
"""

from typing import List, Tuple

import numpy as np
import vtk


class HeatGenerationSolver:
    """
    3D steady-state heat conduction solver with internal heat generation.

    Solves the Poisson equation ∇²T = -Q/k in a cubic domain with
    Dirichlet boundary conditions (T=0) on all six faces.

    Attributes:
        nx, ny, nz: Number of grid points in each direction.
        length: Physical size of the cubic domain (m).
        q_over_k: Heat generation rate divided by thermal conductivity (K/m²).
        omega: SOR relaxation parameter.
        temperature: 3D numpy array of temperature values.
    """

    def __init__(
        self,
        nx: int = 50,
        ny: int = 50,
        nz: int = 50,
        length: float = 1.0,
        q_over_k: float = 100.0,
        omega: float = 1.6,
    ):
        """
        Initialize the heat generation solver.

        Args:
            nx: Number of grid points in x-direction.
            ny: Number of grid points in y-direction.
            nz: Number of grid points in z-direction.
            length: Physical size of cubic domain (m).
            q_over_k: Ratio of heat generation to conductivity Q/k (K/m²).
            omega: SOR relaxation parameter (optimal ~1.5-1.9).
        """
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.length = length
        self.q_over_k = q_over_k
        self.omega = omega

        # Grid spacing (uniform)
        self.h = length / (nx - 1)

        # Source term contribution: h²·Q/k
        self.source_term = self.h**2 * q_over_k

        # Initialize temperature field to zero
        self.temperature = np.zeros((nx, ny, nz))

    def _apply_boundary_conditions(self) -> None:
        """
        Apply Dirichlet boundary conditions (T=0) on all six faces.

        For heat generation problems, the boundaries act as heat sinks,
        removing the internally generated heat.
        """
        self.temperature[0, :, :] = 0.0  # x = 0
        self.temperature[-1, :, :] = 0.0  # x = L
        self.temperature[:, 0, :] = 0.0  # y = 0
        self.temperature[:, -1, :] = 0.0  # y = L
        self.temperature[:, :, 0] = 0.0  # z = 0
        self.temperature[:, :, -1] = 0.0  # z = L

    def solve(self, max_iterations: int = 5000, tolerance: float = 1e-5) -> int:
        """
        Solve the Poisson equation using Gauss-Seidel with SOR.

        The iteration formula includes the source term from heat generation:
            T_new = (1/6)(sum_of_neighbors + h²·Q/k)

        Args:
            max_iterations: Maximum number of iterations.
            tolerance: Convergence criterion (max equation residual).

        Returns:
            Number of iterations to convergence.
        """
        print("Solving 3D heat equation with internal heat generation...")
        print(f"Grid: {self.nx} x {self.ny} x {self.nz}")
        print(f"Q/k = {self.q_over_k} K/m²")
        print(f"SOR relaxation parameter: ω = {self.omega}")
        print("-" * 50)

        for iteration in range(1, max_iterations + 1):
            max_residual = 0.0

            for i in range(1, self.nx - 1):
                for j in range(1, self.ny - 1):
                    for k in range(1, self.nz - 1):
                        # Sum of six neighbors plus source term
                        neighbor_sum = (
                            self.temperature[i + 1, j, k]
                            + self.temperature[i - 1, j, k]
                            + self.temperature[i, j + 1, k]
                            + self.temperature[i, j - 1, k]
                            + self.temperature[i, j, k + 1]
                            + self.temperature[i, j, k - 1]
                        )

                        # New value from Poisson equation
                        t_new = (neighbor_sum + self.source_term) / 6.0

                        # Compute residual before update
                        residual = abs(t_new - self.temperature[i, j, k])
                        max_residual = max(max_residual, residual)

                        # SOR update
                        self.temperature[i, j, k] = (1 - self.omega) * self.temperature[
                            i, j, k
                        ] + self.omega * t_new

            # Apply boundary conditions
            self._apply_boundary_conditions()

            if iteration % 100 == 0:
                print(f"Iteration {iteration:4d}: max residual = {max_residual:.2e}")

            if max_residual < tolerance:
                print(f"\nConverged after {iteration} iterations")
                print(f"Final residual: {max_residual:.2e}")
                print(f"Max temperature: {np.max(self.temperature):.2f}°C")
                return iteration

        print(f"\nMax iterations ({max_iterations}) reached")
        return max_iterations


def create_vtk_image_data(
    temperature: np.ndarray, spacing: Tuple[float, float, float] = (1.0, 1.0, 1.0)
) -> vtk.vtkImageData:
    """
    Convert numpy temperature array to VTK image data for isosurface extraction.

    Args:
        temperature: 3D numpy array of temperature values.
        spacing: Physical spacing between grid points (dx, dy, dz).

    Returns:
        vtk.vtkImageData containing the temperature field.
    """
    nx, ny, nz = temperature.shape
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(nx, ny, nz)
    image_data.SetSpacing(*spacing)
    image_data.SetOrigin(0, 0, 0)

    # Flatten temperature array in Fortran order for VTK
    flat_temp = temperature.flatten(order="F").astype(np.float32)

    scalars = vtk.vtkFloatArray()
    scalars.SetName("Temperature")
    scalars.SetNumberOfValues(len(flat_temp))
    for i, val in enumerate(flat_temp):
        scalars.SetValue(i, val)

    image_data.GetPointData().SetScalars(scalars)
    return image_data


def create_isosurfaces(data: vtk.vtkImageData, isovalues: List[float]) -> vtk.vtkActor:
    """
    Extract isosurfaces (isotherms) at specified temperature values.

    Uses the Marching Cubes algorithm via vtkContourFilter to generate
    polygonal surfaces where the scalar field equals each isovalue.

    Args:
        data: VTK image data containing the temperature field.
        isovalues: List of temperature values for isosurface extraction.

    Returns:
        vtkActor representing the isosurfaces with color mapping.
    """
    contour = vtk.vtkContourFilter()
    contour.SetInputData(data)
    contour.GenerateValues(len(isovalues), min(isovalues), max(isovalues))

    # Compute normals for better shading
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(contour.GetOutputPort())
    normals.ComputePointNormalsOn()
    normals.SplittingOff()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.SetScalarRange(min(isovalues), max(isovalues))

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetOpacity(0.7)

    return actor


def create_color_map(t_min: float, t_max: float) -> vtk.vtkColorTransferFunction:
    """
    Create a cool-to-warm color transfer function for temperature.

    Maps temperature values to colors:
    - Blue (cold) → Cyan → Green → Yellow → Red (hot)

    Args:
        t_min: Minimum temperature value.
        t_max: Maximum temperature value.

    Returns:
        vtkColorTransferFunction for temperature-to-color mapping.
    """
    color_tf = vtk.vtkColorTransferFunction()
    t_range = t_max - t_min

    color_tf.AddRGBPoint(t_min, 0.0, 0.0, 1.0)  # Blue
    color_tf.AddRGBPoint(t_min + 0.25 * t_range, 0.0, 1.0, 1.0)  # Cyan
    color_tf.AddRGBPoint(t_min + 0.5 * t_range, 0.0, 1.0, 0.0)  # Green
    color_tf.AddRGBPoint(t_min + 0.75 * t_range, 1.0, 1.0, 0.0)  # Yellow
    color_tf.AddRGBPoint(t_max, 1.0, 0.0, 0.0)  # Red

    return color_tf


def visualize_isosurfaces(
    actor: vtk.vtkActor,
    t_min: float,
    t_max: float,
) -> None:
    """
    Visualize isosurfaces with color legend and interactive controls.

    Sets up a complete VTK visualization pipeline with:
    - Renderer with dark background for contrast
    - Scalar bar showing temperature scale
    - Orientation axes widget
    - Trackball camera interaction

    Args:
        actor: The isosurface actor to visualize.
        t_min: Minimum temperature for color scale.
        t_max: Maximum temperature for color scale.
    """
    # Apply color map
    color_tf = create_color_map(t_min, t_max)
    actor.GetMapper().SetLookupTable(color_tf)

    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.1, 0.1, 0.15)

    # Create scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_tf)
    scalar_bar.SetTitle("Temperature (°C)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.08)
    scalar_bar.SetHeight(0.4)
    scalar_bar.SetPosition(0.88, 0.3)
    renderer.AddActor(scalar_bar)

    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1024, 768)
    render_window.SetWindowName("3D Temperature Isosurfaces (Isotherms)")

    # Create interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

    # Add orientation axes
    axes = vtk.vtkAxesActor()
    widget = vtk.vtkOrientationMarkerWidget()
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(interactor)
    widget.SetViewport(0, 0, 0.2, 0.2)
    widget.SetEnabled(1)
    widget.InteractiveOff()

    # Set up camera
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    camera.Azimuth(30)
    camera.Elevation(20)

    # Start visualization
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


def main():
    """
    Main function to run the heat generation simulation and isosurface visualization.

    Workflow:
    1. Solve the 3D Poisson equation for temperature with heat generation
    2. Convert solution to VTK image data
    3. Extract isosurfaces at multiple temperature levels
    4. Display with color mapping and interactive controls
    """
    print("=" * 60)
    print("3D Heat Generation - Isosurface Visualization")
    print("=" * 60)
    print()

    # Create and solve the heat generation problem
    solver = HeatGenerationSolver(
        nx=50,
        ny=50,
        nz=50,
        length=1.0,
        q_over_k=500.0,  # Moderate heat generation
        omega=1.6,
    )
    solver.solve(max_iterations=5000, tolerance=1e-5)

    # Convert to VTK format
    h = solver.length / (solver.nx - 1)
    image_data = create_vtk_image_data(solver.temperature, spacing=(h, h, h))

    # Create isosurfaces at several temperature levels
    t_max = np.max(solver.temperature)
    t_min = np.min(solver.temperature)
    isovalues = np.linspace(t_min + 0.1 * t_max, 0.9 * t_max, 8).tolist()

    print(f"\nExtracting {len(isovalues)} isosurfaces...")
    print(f"Temperature range: {t_min:.2f}°C to {t_max:.2f}°C")
    print("Use mouse to rotate, zoom, and pan.")
    print()

    iso_actor = create_isosurfaces(image_data, isovalues)
    visualize_isosurfaces(iso_actor, t_min, t_max)


if __name__ == "__main__":
    main()
