"""
3D Heat Conduction Volume Rendering with VTK

This module demonstrates volume rendering of a physically meaningful CFD problem:
steady-state 3D heat conduction in a cubic domain. Volume rendering is a powerful
visualization technique in CFD for displaying scalar fields (temperature, density,
concentration) throughout a 3D volume without slicing or isosurfaces.

Physical Problem:
-----------------
We solve the steady-state heat equation (Laplace equation) in a cubic domain:

    ∇²T = ∂²T/∂x² + ∂²T/∂y² + ∂²T/∂z² = 0

with Dirichlet boundary conditions:
- T = T_hot (100°C) on the bottom face (z = 0)
- T = T_cold (0°C) on the top face (z = L)
- Adiabatic (zero heat flux) on all side walls

This represents heat conduction through a solid material (e.g., metal block)
with fixed temperatures on top and bottom surfaces.

Numerical Method:
-----------------
The heat equation is discretized using a second-order central finite difference
scheme on a uniform Cartesian grid:

    (T[i+1,j,k] - 2T[i,j,k] + T[i-1,j,k])/dx² +
    (T[i,j+1,k] - 2T[i,j,k] + T[i,j-1,k])/dy² +
    (T[i,j,k+1] - 2T[i,j,k] + T[i,j,k-1])/dz² = 0

For a uniform grid (dx = dy = dz = h), this simplifies to:

    T[i,j,k] = (1/6)(T[i+1,j,k] + T[i-1,j,k] + T[i,j+1,k] +
                     T[i,j-1,k] + T[i,j,k+1] + T[i,j,k-1])

The system is solved iteratively using the Gauss-Seidel method with
Successive Over-Relaxation (SOR) for faster convergence.

Volume Rendering Concepts:
--------------------------
Volume rendering visualizes 3D scalar data by casting rays through the volume
and accumulating color and opacity based on transfer functions:

1. Color Transfer Function: Maps scalar values to RGB colors
   - Blue (0°C) → Cyan → Green → Yellow → Red (100°C)
   - Uses a physically intuitive "cold-to-hot" colormap

2. Opacity Transfer Function: Controls visibility of different temperature ranges
   - Higher opacity for extreme temperatures (boundaries)
   - Lower opacity for intermediate values (reveals internal structure)

3. Shading: Phong shading based on temperature gradients for 3D depth perception

CFD Educational Value:
----------------------
This example demonstrates:
- Finite difference discretization of the Laplace equation
- Iterative solution methods (Gauss-Seidel with SOR)
- Boundary condition implementation in 3D
- Convergence monitoring using residual norms
- Volume rendering for 3D scalar field visualization
- Transfer function design for thermal visualization

References:
-----------
1. Incropera, F.P., et al. "Fundamentals of Heat and Mass Transfer"
2. Patankar, S.V. "Numerical Heat Transfer and Fluid Flow"
3. Drebin, R.A., et al. (1988) "Volume Rendering" Computer Graphics, 22(4)
"""

import numpy as np
import vtk


class HeatConductionSolver:
    """
    3D steady-state heat conduction solver using finite differences.

    Solves ∇²T = 0 in a cubic domain with Dirichlet boundary conditions
    on top and bottom faces and adiabatic (Neumann) conditions on side walls.

    Attributes:
        nx, ny, nz: Number of grid points in each direction.
        length: Physical size of the cubic domain (m).
        t_hot: Temperature at bottom boundary (°C).
        t_cold: Temperature at top boundary (°C).
        omega: SOR relaxation parameter (1 < omega < 2 for over-relaxation).
        temperature: 3D numpy array of temperature values.
    """

    def __init__(
        self,
        nx: int = 50,
        ny: int = 50,
        nz: int = 50,
        length: float = 1.0,
        t_hot: float = 100.0,
        t_cold: float = 0.0,
        omega: float = 1.5,
    ):
        """
        Initialize the heat conduction solver.

        Args:
            nx: Number of grid points in x-direction.
            ny: Number of grid points in y-direction.
            nz: Number of grid points in z-direction.
            length: Physical size of cubic domain (m).
            t_hot: Hot boundary temperature at z=0 (°C).
            t_cold: Cold boundary temperature at z=L (°C).
            omega: SOR relaxation parameter (optimal is ~1.5-1.9).
        """
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.length = length
        self.t_hot = t_hot
        self.t_cold = t_cold
        self.omega = omega

        # Grid spacing
        self.dx = length / (nx - 1)
        self.dy = length / (ny - 1)
        self.dz = length / (nz - 1)

        # Initialize temperature field with linear interpolation (good initial guess)
        self.temperature = self._initialize_temperature()

    def _initialize_temperature(self) -> np.ndarray:
        """
        Initialize temperature field with linear profile between boundaries.

        A linear initial guess reduces the number of iterations needed
        for convergence compared to a uniform initial condition.

        Returns:
            3D numpy array with initial temperature distribution.
        """
        temperature = np.zeros((self.nx, self.ny, self.nz))

        # Linear interpolation from t_hot (z=0) to t_cold (z=L)
        for k in range(self.nz):
            z_ratio = k / (self.nz - 1)
            temperature[:, :, k] = self.t_hot + (self.t_cold - self.t_hot) * z_ratio

        return temperature

    def _apply_boundary_conditions(self) -> None:
        """
        Apply boundary conditions to the temperature field.

        Dirichlet conditions:
        - Bottom face (k=0): T = t_hot
        - Top face (k=nz-1): T = t_cold

        Neumann conditions (adiabatic, ∂T/∂n = 0):
        - Side walls: copy from adjacent interior cells
        """
        # Dirichlet: fixed temperature on top and bottom
        self.temperature[:, :, 0] = self.t_hot  # Bottom (hot)
        self.temperature[:, :, -1] = self.t_cold  # Top (cold)

        # Neumann: zero gradient on side walls (adiabatic)
        self.temperature[0, :, :] = self.temperature[1, :, :]  # x = 0
        self.temperature[-1, :, :] = self.temperature[-2, :, :]  # x = L
        self.temperature[:, 0, :] = self.temperature[:, 1, :]  # y = 0
        self.temperature[:, -1, :] = self.temperature[:, -2, :]  # y = L

    def solve(self, max_iterations: int = 5000, tolerance: float = 1e-5) -> int:
        """
        Solve the steady-state heat equation using Gauss-Seidel with SOR.

        The Gauss-Seidel method updates each point using the latest available
        values, providing faster convergence than Jacobi iteration.
        SOR accelerates convergence by over-relaxing the update.

        Update formula:
            T_new = (1-ω)T_old + (ω/6)(T_neighbors_sum)

        Args:
            max_iterations: Maximum number of iterations.
            tolerance: Convergence criterion (max equation residual).

        Returns:
            Number of iterations to convergence.
        """
        print("Solving 3D heat conduction equation...")
        print(f"Grid: {self.nx} x {self.ny} x {self.nz}")
        print(f"T_hot = {self.t_hot}°C, T_cold = {self.t_cold}°C")
        print(f"SOR relaxation parameter: ω = {self.omega}")
        print("-" * 50)

        for iteration in range(1, max_iterations + 1):
            max_residual = 0.0

            # Update interior points using Gauss-Seidel with SOR
            for i in range(1, self.nx - 1):
                for j in range(1, self.ny - 1):
                    for k in range(1, self.nz - 1):
                        # Sum of six neighbors
                        neighbor_sum = (
                            self.temperature[i + 1, j, k]
                            + self.temperature[i - 1, j, k]
                            + self.temperature[i, j + 1, k]
                            + self.temperature[i, j - 1, k]
                            + self.temperature[i, j, k + 1]
                            + self.temperature[i, j, k - 1]
                        )

                        # New value from Laplace equation
                        t_new = neighbor_sum / 6.0

                        # Compute equation residual before update
                        # Residual measures how well current value satisfies Laplace equation
                        residual = abs(t_new - self.temperature[i, j, k])
                        max_residual = max(max_residual, residual)

                        # SOR update
                        self.temperature[i, j, k] = (
                            1 - self.omega
                        ) * self.temperature[i, j, k] + self.omega * t_new

            # Apply boundary conditions after each iteration
            self._apply_boundary_conditions()

            # Check convergence
            if iteration % 100 == 0:
                print(f"Iteration {iteration:4d}: max residual = {max_residual:.2e}")

            if max_residual < tolerance:
                print(f"\nConverged after {iteration} iterations")
                print(f"Final residual: {max_residual:.2e}")
                return iteration

        print(f"\nMax iterations ({max_iterations}) reached")
        print(f"Final residual: {max_residual:.2e}")
        return max_iterations


class VolumeRenderer:
    """
    VTK-based volume renderer for 3D scalar field visualization.

    Creates a volume rendering of a 3D temperature field using ray casting
    with configurable transfer functions for color and opacity mapping.

    Attributes:
        temperature: 3D numpy array of temperature values.
        t_min: Minimum temperature for transfer function scaling.
        t_max: Maximum temperature for transfer function scaling.
    """

    def __init__(self, temperature: np.ndarray, t_min: float = 0.0, t_max: float = 100.0):
        """
        Initialize the volume renderer.

        Args:
            temperature: 3D numpy array of temperature values.
            t_min: Minimum temperature for color/opacity mapping.
            t_max: Maximum temperature for color/opacity mapping.
        """
        self.temperature = temperature
        self.t_min = t_min
        self.t_max = t_max

        # Create VTK pipeline components
        self.volume_data = self._create_volume_data()
        self.opacity_tf = self._create_opacity_transfer_function()
        self.color_tf = self._create_color_transfer_function()
        self.volume_property = self._create_volume_property()
        self.volume_mapper = self._create_volume_mapper()
        self.volume_actor = self._create_volume_actor()

    def _create_volume_data(self) -> vtk.vtkImageData:
        """
        Convert numpy temperature array to VTK image data.

        Uses VTK_FLOAT for accurate representation of temperature values.
        Sets spacing to maintain proper aspect ratio.

        Returns:
            vtkImageData object containing the temperature field.
        """
        nx, ny, nz = self.temperature.shape
        volume_data = vtk.vtkImageData()
        volume_data.SetDimensions(nx, ny, nz)
        volume_data.SetSpacing(1.0 / (nx - 1), 1.0 / (ny - 1), 1.0 / (nz - 1))
        volume_data.SetOrigin(0, 0, 0)
        volume_data.AllocateScalars(vtk.VTK_FLOAT, 1)

        # Copy temperature data to VTK format
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    volume_data.SetScalarComponentFromFloat(
                        i, j, k, 0, self.temperature[i, j, k]
                    )

        return volume_data

    def _create_opacity_transfer_function(self) -> vtk.vtkPiecewiseFunction:
        """
        Create opacity transfer function for temperature visualization.

        Design philosophy for heat conduction visualization:
        - Low opacity for intermediate temperatures (reveals internal structure)
        - Higher opacity near hot and cold boundaries (emphasizes sources/sinks)
        - Smooth transitions prevent visual artifacts

        Returns:
            vtkPiecewiseFunction mapping temperature to opacity [0, 1].
        """
        opacity_tf = vtk.vtkPiecewiseFunction()

        # Temperature range parameters
        t_range = self.t_max - self.t_min
        t_25 = self.t_min + 0.25 * t_range
        t_50 = self.t_min + 0.50 * t_range
        t_75 = self.t_min + 0.75 * t_range

        # Opacity profile: higher at boundaries, lower in middle
        opacity_tf.AddPoint(self.t_min, 0.4)  # Cold boundary
        opacity_tf.AddPoint(t_25, 0.15)  # Transition
        opacity_tf.AddPoint(t_50, 0.08)  # Middle (most transparent)
        opacity_tf.AddPoint(t_75, 0.15)  # Transition
        opacity_tf.AddPoint(self.t_max, 0.4)  # Hot boundary

        return opacity_tf

    def _create_color_transfer_function(self) -> vtk.vtkColorTransferFunction:
        """
        Create color transfer function for thermal visualization.

        Uses a physically intuitive "cool-to-warm" colormap:
        - Blue (cold) → Cyan → Green → Yellow → Red (hot)
        - This rainbow-like progression is widely recognized in thermal analysis

        Returns:
            vtkColorTransferFunction mapping temperature to RGB colors.
        """
        color_tf = vtk.vtkColorTransferFunction()

        # Temperature range parameters
        t_range = self.t_max - self.t_min
        t_25 = self.t_min + 0.25 * t_range
        t_50 = self.t_min + 0.50 * t_range
        t_75 = self.t_min + 0.75 * t_range

        # Cool-to-warm colormap
        color_tf.AddRGBPoint(self.t_min, 0.0, 0.0, 1.0)  # Blue (cold)
        color_tf.AddRGBPoint(t_25, 0.0, 1.0, 1.0)  # Cyan
        color_tf.AddRGBPoint(t_50, 0.0, 1.0, 0.0)  # Green
        color_tf.AddRGBPoint(t_75, 1.0, 1.0, 0.0)  # Yellow
        color_tf.AddRGBPoint(self.t_max, 1.0, 0.0, 0.0)  # Red (hot)

        return color_tf

    def _create_volume_property(self) -> vtk.vtkVolumeProperty:
        """
        Configure volume rendering properties.

        Enables Phong shading based on scalar gradients to enhance
        3D depth perception and surface details.

        Returns:
            vtkVolumeProperty with shading and transfer functions.
        """
        volume_property = vtk.vtkVolumeProperty()
        volume_property.SetColor(self.color_tf)
        volume_property.SetScalarOpacity(self.opacity_tf)
        volume_property.SetInterpolationTypeToLinear()

        # Enable gradient-based shading for 3D appearance
        volume_property.ShadeOn()
        volume_property.SetAmbient(0.4)
        volume_property.SetDiffuse(0.6)
        volume_property.SetSpecular(0.2)
        volume_property.SetSpecularPower(10.0)

        return volume_property

    def _create_volume_mapper(self) -> vtk.vtkSmartVolumeMapper:
        """
        Create volume mapper for ray casting.

        vtkSmartVolumeMapper automatically selects the best available
        rendering method (GPU ray casting, CPU ray casting, etc.).

        Returns:
            vtkSmartVolumeMapper connected to volume data.
        """
        volume_mapper = vtk.vtkSmartVolumeMapper()
        volume_mapper.SetInputData(self.volume_data)
        volume_mapper.SetBlendModeToComposite()

        return volume_mapper

    def _create_volume_actor(self) -> vtk.vtkVolume:
        """
        Create the volume actor for rendering.

        Returns:
            vtkVolume with mapper and properties configured.
        """
        volume_actor = vtk.vtkVolume()
        volume_actor.SetMapper(self.volume_mapper)
        volume_actor.SetProperty(self.volume_property)

        return volume_actor

    def _create_scalar_bar(self) -> vtk.vtkScalarBarActor:
        """
        Create a color legend (scalar bar) for the visualization.

        Returns:
            vtkScalarBarActor showing temperature-to-color mapping.
        """
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(self.color_tf)
        scalar_bar.SetTitle("Temperature (°C)")
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetWidth(0.08)
        scalar_bar.SetHeight(0.4)
        scalar_bar.SetPosition(0.88, 0.3)

        # Customize text properties
        title_prop = scalar_bar.GetTitleTextProperty()
        title_prop.SetFontSize(12)
        title_prop.SetColor(1, 1, 1)

        label_prop = scalar_bar.GetLabelTextProperty()
        label_prop.SetFontSize(10)
        label_prop.SetColor(1, 1, 1)

        return scalar_bar

    def _create_axes_widget(
        self, interactor: vtk.vtkRenderWindowInteractor
    ) -> vtk.vtkOrientationMarkerWidget:
        """
        Create an orientation axes widget for spatial reference.

        Args:
            interactor: The render window interactor.

        Returns:
            vtkOrientationMarkerWidget with axes display.
        """
        axes = vtk.vtkAxesActor()
        axes.SetTotalLength(0.3, 0.3, 0.3)

        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        widget.SetInteractor(interactor)
        widget.SetViewport(0.0, 0.0, 0.2, 0.2)
        widget.SetEnabled(1)
        widget.InteractiveOff()

        return widget

    def render(self) -> None:
        """
        Display the volume rendering with interactive controls.

        Sets up:
        - Renderer with dark background for contrast
        - Scalar bar legend
        - Orientation axes widget
        - Interactive camera controls (trackball style)
        """
        # Create renderer
        renderer = vtk.vtkRenderer()
        renderer.AddVolume(self.volume_actor)
        renderer.AddActor(self._create_scalar_bar())
        renderer.SetBackground(0.1, 0.1, 0.15)

        # Create render window
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(1024, 768)
        render_window.SetWindowName("3D Heat Conduction - Volume Rendering")

        # Create interactor
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)

        # Set interaction style
        style = vtk.vtkInteractorStyleTrackballCamera()
        interactor.SetInteractorStyle(style)

        # Add axes widget
        axes_widget = self._create_axes_widget(interactor)

        # Position camera for good initial view
        camera = renderer.GetActiveCamera()
        camera.SetPosition(2.0, 2.0, 1.5)
        camera.SetFocalPoint(0.5, 0.5, 0.5)
        camera.SetViewUp(0, 0, 1)
        renderer.ResetCamera()

        # Start rendering
        interactor.Initialize()
        render_window.Render()
        interactor.Start()

        # Cleanup
        render_window.Finalize()
        interactor.TerminateApp()


def main():
    """
    Main function to run the 3D heat conduction simulation and visualization.

    Steps:
    1. Create and solve the heat conduction problem
    2. Initialize volume renderer with computed temperature field
    3. Display interactive volume rendering
    """
    print("=" * 60)
    print("3D Heat Conduction Volume Rendering")
    print("=" * 60)
    print()

    # Create solver with physically meaningful parameters
    solver = HeatConductionSolver(
        nx=50,  # Grid resolution
        ny=50,
        nz=50,
        length=1.0,  # 1m cube
        t_hot=100.0,  # Bottom boundary (°C)
        t_cold=0.0,  # Top boundary (°C)
        omega=1.6,  # SOR parameter (optimal ~ 1.5-1.9)
    )

    # Solve the heat equation
    solver.solve(max_iterations=5000, tolerance=1e-5)

    print()
    print("Starting volume rendering visualization...")
    print("Use mouse to rotate, zoom, and pan the volume.")
    print()

    # Create and display volume rendering
    renderer = VolumeRenderer(
        temperature=solver.temperature, t_min=solver.t_cold, t_max=solver.t_hot
    )
    renderer.render()


if __name__ == "__main__":
    main()
