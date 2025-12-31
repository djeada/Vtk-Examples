"""
Lid-Driven Cavity Flow Streamline Visualization with VTK

This module demonstrates streamline visualization of a physically meaningful CFD
problem: the 2D lid-driven cavity flow extended to 3D. Streamlines are fundamental
tools in fluid dynamics for visualizing flow patterns, identifying vortices, and
understanding mass transport.

Physical Problem:
-----------------
The lid-driven cavity is a classical benchmark problem in computational fluid dynamics.
A rectangular/cubic cavity has:
- Three stationary walls (left, right, bottom) with no-slip conditions: u = v = w = 0
- A top wall (lid) moving at constant velocity U_lid in the x-direction

The governing equations are the incompressible Navier-Stokes equations:

    ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u   (momentum)
    ∇·u = 0                          (continuity/incompressibility)

where:
- u = (u, v, w) is the velocity vector field
- p = pressure
- ρ = fluid density
- ν = kinematic viscosity

Reynolds Number:
    Re = U_lid * L / ν

The flow pattern depends strongly on Re:
- Re < 100: Single primary vortex, laminar
- Re ~ 100-1000: Primary vortex with corner eddies
- Re > ~8000: Transition to turbulence

Numerical Method:
-----------------
For this educational example, we use a simplified steady-state formulation.
The velocity field is computed using the Artificial Compressibility Method (ACM):

    (1/β) ∂p/∂t + ∇·u = 0   (modified continuity)

where β is the artificial compressibility parameter. As the solution converges
to steady state, ∂p/∂t → 0, recovering true incompressibility.

The momentum equations are solved using:
- Second-order central differences for diffusion and pressure gradient
- First-order upwind for convection (stable for all Re)
- Explicit time stepping to steady state

Streamline Concepts:
--------------------
Streamlines are curves that are everywhere tangent to the instantaneous velocity
field. Key properties:

1. Definition: A streamline satisfies dr/ds = u(r) where r is position and s is
   arc length along the streamline.

2. Visualization: VTK's vtkStreamTracer integrates the velocity field to trace
   streamlines from seed points.

3. Physical Meaning:
   - In steady flow, streamlines = pathlines = streaklines
   - Fluid particles travel along streamlines
   - Streamlines never cross (in steady flow)

4. Seed Points: Starting positions for streamline integration. Strategic
   placement reveals flow structure.

CFD Educational Value:
----------------------
This example demonstrates:
- Navier-Stokes equations for incompressible flow
- Lid-driven cavity as a CFD benchmark
- Artificial Compressibility Method
- No-slip boundary conditions
- Reynolds number and its effect on flow
- Streamline visualization techniques
- Vortex identification through streamline patterns

References:
-----------
1. Ghia, U., Ghia, K.N., & Shin, C.T. (1982) "High-Re solutions for incompressible
   flow using the Navier-Stokes equations" J. Comp. Physics 48(3), 387-411
2. Chorin, A.J. (1967) "A numerical method for solving incompressible viscous
   flow problems" J. Comp. Physics 2(1), 12-26
3. Kundu, P.K., Cohen, I.M. "Fluid Mechanics"
"""

import numpy as np
import vtk


class CavityFlowSolver:
    """
    2D Lid-driven cavity flow solver using Artificial Compressibility Method.

    Solves the incompressible Navier-Stokes equations in a square cavity with
    a moving top lid. The solution is extended to 3D for visualization.

    Attributes:
        nx, ny: Number of grid points in each direction.
        length: Cavity size (m).
        u_lid: Lid velocity (m/s).
        reynolds: Reynolds number.
        u, v: Velocity component arrays.
        p: Pressure array.
    """

    def __init__(
        self,
        nx: int = 41,
        ny: int = 41,
        length: float = 1.0,
        reynolds: float = 100.0,
        u_lid: float = 1.0,
    ):
        """
        Initialize the cavity flow solver.

        Args:
            nx: Number of grid points in x-direction.
            ny: Number of grid points in y-direction.
            length: Cavity size (m).
            reynolds: Reynolds number Re = U*L/ν.
            u_lid: Lid velocity (m/s).
        """
        self.nx = nx
        self.ny = ny
        self.length = length
        self.reynolds = reynolds
        self.u_lid = u_lid

        # Compute kinematic viscosity from Reynolds number
        self.nu = u_lid * length / reynolds

        # Grid spacing
        self.dx = length / (nx - 1)
        self.dy = length / (ny - 1)

        # Artificial compressibility parameter
        self.beta = max(1.0, u_lid**2)

        # Initialize velocity and pressure fields
        self.u = np.zeros((ny, nx))
        self.v = np.zeros((ny, nx))
        self.p = np.zeros((ny, nx))

        # Apply initial boundary conditions
        self._apply_boundary_conditions()

    def _apply_boundary_conditions(self) -> None:
        """
        Apply no-slip boundary conditions on walls and moving lid.

        - Left, right, bottom walls: u = v = 0
        - Top wall (lid): u = u_lid, v = 0
        """
        # Left and right walls
        self.u[:, 0] = 0.0
        self.u[:, -1] = 0.0
        self.v[:, 0] = 0.0
        self.v[:, -1] = 0.0

        # Bottom wall
        self.u[0, :] = 0.0
        self.v[0, :] = 0.0

        # Top wall (moving lid)
        self.u[-1, :] = self.u_lid
        self.v[-1, :] = 0.0

    def _apply_pressure_bc(self) -> None:
        """Apply Neumann (zero-gradient) pressure boundary conditions."""
        self.p[:, 0] = self.p[:, 1]
        self.p[:, -1] = self.p[:, -2]
        self.p[0, :] = self.p[1, :]
        self.p[-1, :] = self.p[-2, :]

    def _compute_timestep(self) -> float:
        """
        Compute stable time step based on CFL and diffusion constraints.

        Returns:
            Stable time step (s).
        """
        h = min(self.dx, self.dy)
        c = np.sqrt(self.beta)

        # CFL constraint
        dt_cfl = 0.2 * h / (self.u_lid + c)

        # Diffusion constraint
        dt_diff = 0.2 * h**2 / (4 * self.nu)

        return min(dt_cfl, dt_diff)

    def solve(self, max_iterations: int = 5000, tolerance: float = 1e-6) -> int:
        """
        Solve the Navier-Stokes equations using Artificial Compressibility.

        Uses explicit time stepping to reach steady state. The ACM introduces
        a pseudo-pressure time derivative that vanishes at convergence.

        Args:
            max_iterations: Maximum number of iterations.
            tolerance: Convergence criterion (max velocity change).

        Returns:
            Number of iterations to convergence.
        """
        print("Solving lid-driven cavity flow...")
        print(f"Grid: {self.nx} x {self.ny}")
        print(f"Reynolds number: {self.reynolds}")
        print(f"Lid velocity: {self.u_lid} m/s")
        print(f"Kinematic viscosity: {self.nu:.6f} m²/s")
        print("-" * 50)

        dt = self._compute_timestep()

        for iteration in range(1, max_iterations + 1):
            u_old = self.u.copy()
            v_old = self.v.copy()

            # Update interior points
            for j in range(1, self.ny - 1):
                for i in range(1, self.nx - 1):
                    # Convection (upwind)
                    if self.u[j, i] >= 0:
                        dudx = (self.u[j, i] - self.u[j, i - 1]) / self.dx
                    else:
                        dudx = (self.u[j, i + 1] - self.u[j, i]) / self.dx

                    if self.v[j, i] >= 0:
                        dudy = (self.u[j, i] - self.u[j - 1, i]) / self.dy
                    else:
                        dudy = (self.u[j + 1, i] - self.u[j, i]) / self.dy

                    if self.u[j, i] >= 0:
                        dvdx = (self.v[j, i] - self.v[j, i - 1]) / self.dx
                    else:
                        dvdx = (self.v[j, i + 1] - self.v[j, i]) / self.dx

                    if self.v[j, i] >= 0:
                        dvdy = (self.v[j, i] - self.v[j - 1, i]) / self.dy
                    else:
                        dvdy = (self.v[j + 1, i] - self.v[j, i]) / self.dy

                    conv_u = self.u[j, i] * dudx + self.v[j, i] * dudy
                    conv_v = self.u[j, i] * dvdx + self.v[j, i] * dvdy

                    # Diffusion (central)
                    diff_u = self.nu * (
                        (self.u[j, i + 1] - 2 * self.u[j, i] + self.u[j, i - 1])
                        / self.dx**2
                        + (self.u[j + 1, i] - 2 * self.u[j, i] + self.u[j - 1, i])
                        / self.dy**2
                    )
                    diff_v = self.nu * (
                        (self.v[j, i + 1] - 2 * self.v[j, i] + self.v[j, i - 1])
                        / self.dx**2
                        + (self.v[j + 1, i] - 2 * self.v[j, i] + self.v[j - 1, i])
                        / self.dy**2
                    )

                    # Pressure gradient
                    dpdx = (self.p[j, i + 1] - self.p[j, i - 1]) / (2 * self.dx)
                    dpdy = (self.p[j + 1, i] - self.p[j - 1, i]) / (2 * self.dy)

                    # Divergence
                    div = (self.u[j, i + 1] - self.u[j, i - 1]) / (
                        2 * self.dx
                    ) + (self.v[j + 1, i] - self.v[j - 1, i]) / (2 * self.dy)

                    # Update
                    self.u[j, i] = u_old[j, i] + dt * (-conv_u - dpdx + diff_u)
                    self.v[j, i] = v_old[j, i] + dt * (-conv_v - dpdy + diff_v)
                    self.p[j, i] = self.p[j, i] - dt * self.beta * div

            self._apply_boundary_conditions()
            self._apply_pressure_bc()

            # Check convergence
            max_change = max(
                np.max(np.abs(self.u - u_old)), np.max(np.abs(self.v - v_old))
            )

            if iteration % 500 == 0:
                print(f"Iteration {iteration:4d}: max change = {max_change:.2e}")

            if max_change < tolerance and iteration > 100:
                print(f"\nConverged after {iteration} iterations")
                return iteration

        print(f"\nMax iterations ({max_iterations}) reached")
        return max_iterations

    def get_velocity_magnitude(self) -> np.ndarray:
        """Compute velocity magnitude field."""
        return np.sqrt(self.u**2 + self.v**2)


class StreamlineVisualizer:
    """
    VTK-based streamline visualizer for 2D/3D velocity fields.

    Creates streamlines from a velocity field using vtkStreamTracer
    with configurable seed points and integration parameters.

    Attributes:
        u, v: 2D velocity component arrays.
        nx, ny: Grid dimensions.
        length: Physical domain size.
    """

    def __init__(
        self, u: np.ndarray, v: np.ndarray, length: float = 1.0, nz: int = 5
    ):
        """
        Initialize the streamline visualizer.

        Args:
            u: x-velocity component array (ny, nx).
            v: y-velocity component array (ny, nx).
            length: Physical domain size (m).
            nz: Number of z-planes for 3D extension.
        """
        self.u = u
        self.v = v
        self.ny, self.nx = u.shape
        self.length = length
        self.nz = nz

        # Create VTK objects
        self.grid = self._create_structured_grid()
        self.seeds = self._create_seed_points()
        self.streamlines = self._create_streamlines()

    def _create_structured_grid(self) -> vtk.vtkStructuredGrid:
        """
        Create a VTK structured grid with the velocity field.

        Extends the 2D field to 3D by replicating in the z-direction.

        Returns:
            vtkStructuredGrid with velocity vectors.
        """
        points = vtk.vtkPoints()
        vectors = vtk.vtkDoubleArray()
        vectors.SetNumberOfComponents(3)
        vectors.SetName("Velocity")

        dx = self.length / (self.nx - 1)
        dy = self.length / (self.ny - 1)
        dz = self.length / (self.nz - 1) if self.nz > 1 else 1.0

        # Create 3D grid by extending 2D field
        for k in range(self.nz):
            z = k * dz
            for j in range(self.ny):
                y = j * dy
                for i in range(self.nx):
                    x = i * dx
                    points.InsertNextPoint(x, y, z)
                    vectors.InsertNextTuple3(self.u[j, i], self.v[j, i], 0.0)

        grid = vtk.vtkStructuredGrid()
        grid.SetDimensions(self.nx, self.ny, self.nz)
        grid.SetPoints(points)
        grid.GetPointData().SetVectors(vectors)

        return grid

    def _create_seed_points(self) -> vtk.vtkPoints:
        """
        Create seed points for streamline integration.

        Places seeds strategically to capture flow structure:
        - Near the lid to trace the primary vortex
        - Near corners to capture secondary eddies

        Returns:
            vtkPoints containing seed positions.
        """
        seeds = vtk.vtkPoints()
        z_mid = self.length / 2.0

        # Seeds near the top (lid) to trace primary circulation
        for i in range(5):
            x = 0.1 + 0.2 * i
            seeds.InsertNextPoint(x * self.length, 0.9 * self.length, z_mid)

        # Seeds in the middle
        for i in range(5):
            x = 0.1 + 0.2 * i
            seeds.InsertNextPoint(x * self.length, 0.5 * self.length, z_mid)

        # Seeds near the bottom
        for i in range(5):
            x = 0.1 + 0.2 * i
            seeds.InsertNextPoint(x * self.length, 0.1 * self.length, z_mid)

        # Vertical line of seeds
        for j in range(5):
            y = 0.1 + 0.2 * j
            seeds.InsertNextPoint(0.5 * self.length, y * self.length, z_mid)

        return seeds

    def _create_streamlines(self) -> vtk.vtkStreamTracer:
        """
        Create streamline tracer with integration parameters.

        Uses RK45 integration for accurate streamline computation.

        Returns:
            vtkStreamTracer configured for the velocity field.
        """
        # Convert points to polydata for seeding
        seed_polydata = vtk.vtkPolyData()
        seed_polydata.SetPoints(self.seeds)

        streamer = vtk.vtkStreamTracer()
        streamer.SetInputData(self.grid)
        streamer.SetSourceData(seed_polydata)
        streamer.SetMaximumPropagation(self.length * 10)
        streamer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
        streamer.SetInitialIntegrationStep(self.length / 100)
        streamer.SetIntegrationDirectionToBoth()
        streamer.SetIntegratorTypeToRungeKutta45()
        streamer.SetComputeVorticity(True)

        return streamer

    def visualize(self) -> None:
        """
        Display streamlines with velocity magnitude coloring.

        Sets up VTK visualization with:
        - Streamlines colored by velocity magnitude
        - Scalar bar legend
        - Cavity outline for context
        - Interactive camera controls
        """
        # Create streamline tubes for better visibility
        tubes = vtk.vtkTubeFilter()
        tubes.SetInputConnection(self.streamlines.GetOutputPort())
        tubes.SetRadius(self.length / 100)
        tubes.SetNumberOfSides(8)

        # Compute velocity magnitude range
        vel_mag = np.sqrt(self.u**2 + self.v**2)
        v_min, v_max = 0.0, np.max(vel_mag)

        # Create color map
        color_tf = vtk.vtkColorTransferFunction()
        color_tf.AddRGBPoint(v_min, 0.0, 0.0, 1.0)
        color_tf.AddRGBPoint(v_max * 0.5, 0.0, 1.0, 0.0)
        color_tf.AddRGBPoint(v_max, 1.0, 0.0, 0.0)

        # Mapper for streamlines
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tubes.GetOutputPort())
        mapper.SetLookupTable(color_tf)
        mapper.SetScalarRange(v_min, v_max)

        stream_actor = vtk.vtkActor()
        stream_actor.SetMapper(mapper)

        # Create cavity outline
        outline = vtk.vtkOutlineFilter()
        outline.SetInputData(self.grid)

        outline_mapper = vtk.vtkPolyDataMapper()
        outline_mapper.SetInputConnection(outline.GetOutputPort())

        outline_actor = vtk.vtkActor()
        outline_actor.SetMapper(outline_mapper)
        outline_actor.GetProperty().SetColor(0.5, 0.5, 0.5)

        # Create renderer
        renderer = vtk.vtkRenderer()
        renderer.AddActor(stream_actor)
        renderer.AddActor(outline_actor)
        renderer.SetBackground(0.1, 0.1, 0.15)

        # Add scalar bar
        scalar_bar = vtk.vtkScalarBarActor()
        scalar_bar.SetLookupTable(color_tf)
        scalar_bar.SetTitle("Velocity (m/s)")
        scalar_bar.SetNumberOfLabels(5)
        scalar_bar.SetWidth(0.08)
        scalar_bar.SetHeight(0.4)
        scalar_bar.SetPosition(0.88, 0.3)
        renderer.AddActor(scalar_bar)

        # Create render window
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(1024, 768)
        render_window.SetWindowName("Lid-Driven Cavity - Streamlines")

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

        # Position camera
        camera = renderer.GetActiveCamera()
        camera.SetPosition(0.5, 0.5, 3.0)
        camera.SetFocalPoint(0.5, 0.5, 0.25)
        camera.SetViewUp(0, 1, 0)
        renderer.ResetCamera()

        # Start visualization
        interactor.Initialize()
        render_window.Render()
        interactor.Start()


def main():
    """
    Main function to run the lid-driven cavity simulation and streamline visualization.

    Workflow:
    1. Solve the 2D Navier-Stokes equations for cavity flow
    2. Extend the velocity field to 3D for visualization
    3. Trace streamlines from strategic seed points
    4. Display with velocity magnitude coloring
    """
    print("=" * 60)
    print("Lid-Driven Cavity Flow - Streamline Visualization")
    print("=" * 60)
    print()

    # Solve the cavity flow problem
    solver = CavityFlowSolver(
        nx=41,
        ny=41,
        length=1.0,
        reynolds=100.0,
        u_lid=1.0,
    )
    solver.solve(max_iterations=10000, tolerance=1e-6)

    print()
    print("Creating streamline visualization...")
    print("Streamlines show flow pattern driven by the moving lid.")
    print("Use mouse to rotate, zoom, and pan.")
    print()

    # Visualize streamlines
    visualizer = StreamlineVisualizer(
        u=solver.u,
        v=solver.v,
        length=solver.length,
        nz=5,
    )
    visualizer.visualize()


if __name__ == "__main__":
    main()
