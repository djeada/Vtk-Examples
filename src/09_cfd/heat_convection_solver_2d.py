"""
2D Steady-State Convection-Diffusion Solver

Logical Extension of the 1D Convection-Diffusion Solver to Two Dimensions

This module provides a solver for the 2D steady-state convection-diffusion equation,
which models the transport of a scalar quantity (e.g., temperature, concentration)
under the combined effects of convection (advection by fluid flow) and diffusion
in a two-dimensional domain.

Workflow:
1. Initialization: Set up the 2D computational grid, material properties (diffusivity Gamma),
   and flow characteristics (velocity components u, v, density rho).
2. Matrix Assembly: Construct the coefficient matrix using Upwind Difference Scheme (UDS)
   for the convection term to ensure stability.
3. Boundary Conditions: Apply Dirichlet boundary conditions on all edges.
4. Solver: Use iterative Gauss-Seidel method to solve the 2D system.
5. Post-Processing: Visualize the scalar field with contours and flow vectors.

Physics:
The steady-state convection-diffusion equation in 2D is:
    d/dx(rho * u * phi) + d/dy(rho * v * phi) = d/dx(Gamma * dphi/dx) + d/dy(Gamma * dphi/dy)

where:
- phi: transported scalar (temperature, concentration, etc.)
- rho: fluid density (kg/m³)
- u, v: fluid velocity components (m/s)
- Gamma: diffusion coefficient (kg/m·s or W/m·K for heat)

For constant properties with uniform flow (u, v = const), this simplifies to:
    rho * u * dphi/dx + rho * v * dphi/dy = Gamma * (d²phi/dx² + d²phi/dy²)

The Peclet number Pe = (rho * |U| * L) / Gamma characterizes the relative importance
of convection to diffusion. High Pe indicates convection-dominated flow.

Discretization Scheme:
Upwind Difference Scheme (UDS):
- Uses upwind approximation for convection: stable for all Pe
- First-order accurate, introduces numerical diffusion
- For 2D with uniform flow, the coefficients are:
  a_W = D_w + max(F_w, 0), a_E = D_e + max(-F_e, 0)
  a_S = D_s + max(F_s, 0), a_N = D_n + max(-F_n, 0)
  a_P = a_W + a_E + a_S + a_N

Relationship to 1D:
This 2D solver extends the 1D version by:
- Adding y-direction discretization with north/south coefficients
- Supporting 2D velocity field (u, v) instead of just u
- Using 2D iterative solver instead of TDMA
- Providing 3D surface visualization with contours and flow vectors
"""

import matplotlib.pyplot as plt
import numpy as np
import vtk

# ========== Animation and Visualization Constants ==========
ANIMATION_TIMER_MS = 50  # Animation update interval in milliseconds
PARTICLE_RADIUS_FACTOR = 0.8  # Particles confined to this fraction of domain
CONVECTION_SCALE_FACTOR = 0.3  # Scale factor for particle convection animation
DIFFUSION_SCALE_FACTOR = 0.2  # Scale factor for particle diffusion animation
N_PARTICLES_2D = 100  # Number of particles in 2D animation


def solve_convection_diffusion_2d(
    nx: int,
    ny: int,
    Lx: float,
    Ly: float,
    Gamma: float,
    rho: float,
    u: float,
    v: float,
    phi_left: float,
    phi_right: float,
    phi_bottom: float,
    phi_top: float,
    Rtol: float = 1e-6,
    max_iter: int = 50000,
) -> tuple:
    """
    Solves the 2D steady-state convection-diffusion equation using Gauss-Seidel iteration.

    The discretized equation for interior nodes using UDS is:
        a_P * phi[i,j] = a_E * phi[i,j+1] + a_W * phi[i,j-1] + a_N * phi[i+1,j] + a_S * phi[i-1,j]

    :param nx: Number of divisions in x-direction (int)
    :param ny: Number of divisions in y-direction (int)
    :param Lx: Domain length in x-direction (float, meters)
    :param Ly: Domain length in y-direction (float, meters)
    :param Gamma: Diffusion coefficient (float, kg/m·s)
    :param rho: Fluid density (float, kg/m³)
    :param u: x-velocity component (float, m/s)
    :param v: y-velocity component (float, m/s)
    :param phi_left: Left boundary value (float)
    :param phi_right: Right boundary value (float)
    :param phi_bottom: Bottom boundary value (float)
    :param phi_top: Top boundary value (float)
    :param Rtol: Residual tolerance for convergence (float)
    :param max_iter: Maximum number of iterations (int)
    :return: Tuple of (x, y, phi, iterations, Pe_x, Pe_y)
    """
    dx = Lx / nx
    dy = Ly / ny

    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, Ly, ny + 1)

    # Initialize phi with average of boundary values
    phi_avg = (phi_left + phi_right + phi_bottom + phi_top) / 4
    phi = phi_avg * np.ones((ny + 1, nx + 1))

    # Apply Dirichlet boundary conditions
    phi[:, 0] = phi_left      # Left (x = 0)
    phi[:, -1] = phi_right    # Right (x = Lx)
    phi[0, :] = phi_bottom    # Bottom (y = 0)
    phi[-1, :] = phi_top      # Top (y = Ly)

    # Convective fluxes (mass flux per unit area)
    Fx = rho * u  # x-direction
    Fy = rho * v  # y-direction

    # Diffusion conductances
    Dx = Gamma / dx
    Dy = Gamma / dy

    # UDS coefficients for interior nodes
    # For positive flow: upwind is west/south
    a_W = Dx + max(Fx, 0)
    a_E = Dx + max(-Fx, 0)
    a_S = Dy + max(Fy, 0)
    a_N = Dy + max(-Fy, 0)
    a_P = a_W + a_E + a_S + a_N

    # Peclet numbers
    Pe_x = abs(Fx) * Lx / Gamma
    Pe_y = abs(Fy) * Ly / Gamma

    # Gauss-Seidel iteration
    iterations = 0
    Rmax = 1.0

    while Rmax > Rtol and iterations < max_iter:
        Rmax = 0.0
        iterations += 1

        for j in range(1, ny):
            for i in range(1, nx):
                phi_new = (
                    a_E * phi[j, i + 1] +
                    a_W * phi[j, i - 1] +
                    a_N * phi[j + 1, i] +
                    a_S * phi[j - 1, i]
                ) / a_P

                R = abs(phi_new - phi[j, i])
                Rmax = max(R, Rmax)
                phi[j, i] = phi_new

        # Re-apply boundary conditions (they should not change, but for safety)
        phi[:, 0] = phi_left
        phi[:, -1] = phi_right
        phi[0, :] = phi_bottom
        phi[-1, :] = phi_top

        if iterations % 1000 == 0:
            print(f"Iteration {iterations}: Residual = {Rmax:.2e}")

    return x, y, phi, iterations, Pe_x, Pe_y


def plot_solution_2d(
    x: np.ndarray,
    y: np.ndarray,
    phi: np.ndarray,
    u: float,
    v: float,
    Pe_x: float,
    Pe_y: float,
):
    """
    Plots the 2D scalar field with contours and velocity vectors.

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param phi: Scalar field (numpy.ndarray)
    :param u: x-velocity component (float)
    :param v: y-velocity component (float)
    :param Pe_x: Peclet number in x-direction (float)
    :param Pe_y: Peclet number in y-direction (float)
    """
    X, Y = np.meshgrid(x, y)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Filled contour plot
    levels = 50
    cf = axes[0].contourf(X, Y, phi, levels, cmap="coolwarm")
    plt.colorbar(cf, ax=axes[0], label="Phi (Scalar)")
    axes[0].set_xlabel("x (m)", fontsize=12)
    axes[0].set_ylabel("y (m)", fontsize=12)
    axes[0].set_title(f"2D Convection-Diffusion: Pe_x={Pe_x:.1f}, Pe_y={Pe_y:.1f}", fontsize=14)
    axes[0].set_aspect("equal")

    # Add velocity vector
    center_x, center_y = x[len(x) // 2], y[len(y) // 2]
    axes[0].quiver(center_x, center_y, u, v, color='white', scale=10, width=0.02)
    axes[0].annotate(f"u={u}, v={v}", (center_x + 0.1, center_y + 0.1), color='white', fontsize=10)

    # Contour lines with streamlines
    cs = axes[1].contour(X, Y, phi, 20, colors='k', linewidths=0.5)
    axes[1].clabel(cs, inline=True, fontsize=8)
    cf2 = axes[1].contourf(X, Y, phi, levels, cmap="coolwarm", alpha=0.7)
    plt.colorbar(cf2, ax=axes[1], label="Phi (Scalar)")

    # Create velocity field for streamlines
    U = u * np.ones_like(X)
    V = v * np.ones_like(Y)
    axes[1].streamplot(X, Y, U, V, color='blue', linewidth=0.8, arrowsize=1)

    axes[1].set_xlabel("x (m)", fontsize=12)
    axes[1].set_ylabel("y (m)", fontsize=12)
    axes[1].set_title("Contours with Streamlines", fontsize=14)
    axes[1].set_aspect("equal")

    plt.tight_layout()
    plt.show()


def create_color_map(phi: np.ndarray) -> vtk.vtkLookupTable:
    """
    Creates a color map for VTK visualization.

    :param phi: Scalar field values (numpy.ndarray)
    :return: vtkLookupTable
    """
    color_map = vtk.vtkLookupTable()
    min_val, max_val = np.min(phi), np.max(phi)
    color_map.SetRange(min_val, max_val)
    color_map.SetHueRange(0.667, 0)  # Blue to red
    color_map.Build()
    return color_map


def create_2d_surface(
    x: np.ndarray, y: np.ndarray, phi: np.ndarray, height_scale: float = 1.0
) -> vtk.vtkPolyData:
    """
    Create a 3D surface representation of the 2D scalar field.

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param phi: Scalar field (numpy.ndarray)
    :param height_scale: Scale factor for z-height (float)
    :return: vtkPolyData
    """
    points = vtk.vtkPoints()
    phi_values = vtk.vtkFloatArray()
    phi_values.SetName("Phi")

    nx, ny = len(x), len(y)
    phi_min = np.min(phi)
    phi_range = np.max(phi) - phi_min
    if phi_range == 0:
        phi_range = 1.0

    for j in range(ny):
        for i in range(nx):
            z = (phi[j, i] - phi_min) / phi_range * height_scale
            points.InsertNextPoint(x[i], y[j], z)
            phi_values.InsertNextValue(phi[j, i])

    # Create cells (quads)
    cells = vtk.vtkCellArray()
    for j in range(ny - 1):
        for i in range(nx - 1):
            quad = vtk.vtkQuad()
            quad.GetPointIds().SetId(0, j * nx + i)
            quad.GetPointIds().SetId(1, j * nx + (i + 1))
            quad.GetPointIds().SetId(2, (j + 1) * nx + (i + 1))
            quad.GetPointIds().SetId(3, (j + 1) * nx + i)
            cells.InsertNextCell(quad)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(phi_values)

    return polyData


def create_velocity_glyph_field_2d(
    x: np.ndarray, y: np.ndarray, u: float, v: float, n_glyphs: int = 8
) -> vtk.vtkActor:
    """
    Create a field of velocity arrows showing flow direction in 2D.

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param u: x-velocity component (float)
    :param v: y-velocity component (float)
    :param n_glyphs: Number of glyphs per dimension (int)
    :return: vtkActor with velocity glyphs
    """
    Lx = x[-1] - x[0]
    Ly = y[-1] - y[0]
    speed = np.sqrt(u**2 + v**2)

    points = vtk.vtkPoints()
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetName("Velocity")

    for j in range(n_glyphs):
        for i in range(n_glyphs):
            px = Lx * (i + 0.5) / n_glyphs
            py = Ly * (j + 0.5) / n_glyphs
            points.InsertNextPoint(px, py, 0.02)  # Slightly above surface
            vectors.InsertNextTuple3(u, v, 0)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetVectors(vectors)

    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(12)
    arrow.SetShaftResolution(12)

    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputData(polyData)
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    base_scale = 0.1 * min(Lx, Ly)
    glyph.SetScaleFactor(base_scale / max(speed, 0.1))
    glyph.OrientOn()
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.2, 0.6, 1.0)
    actor.GetProperty().SetOpacity(0.9)

    return actor


def create_contour_lines(
    polyData: vtk.vtkPolyData, phi: np.ndarray, n_contours: int = 10
) -> vtk.vtkActor:
    """
    Create contour lines on the surface.

    :param polyData: Surface polydata (vtkPolyData)
    :param phi: Scalar field (numpy.ndarray)
    :param n_contours: Number of contour lines (int)
    :return: vtkActor with contour lines
    """
    contour = vtk.vtkContourFilter()
    contour.SetInputData(polyData)
    contour.GenerateValues(n_contours, np.min(phi), np.max(phi))
    contour.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(contour.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.1, 0.1, 0.1)
    actor.GetProperty().SetLineWidth(2.0)

    return actor


def create_boundary_annotation_2d(
    text: str, position: tuple, color: tuple = (1, 1, 1)
) -> vtk.vtkBillboardTextActor3D:
    """
    Create a 3D text annotation for boundary conditions.

    :param text: Text to display (str)
    :param position: Position (x, y, z)
    :param color: Text color (r, g, b)
    :return: vtkBillboardTextActor3D
    """
    text_actor = vtk.vtkBillboardTextActor3D()
    text_actor.SetInput(text)
    text_actor.SetPosition(position)
    text_actor.GetTextProperty().SetFontSize(12)
    text_actor.GetTextProperty().SetColor(color)
    text_actor.GetTextProperty().SetJustificationToCentered()
    text_actor.GetTextProperty().SetBold(True)

    return text_actor


def create_particle_system_2d(
    Lx: float, Ly: float, n_particles: int = 100
) -> tuple:
    """
    Create a particle system for 2D animation.

    :param Lx: Domain length in x (float)
    :param Ly: Domain length in y (float)
    :param n_particles: Number of particles (int)
    :return: Tuple of (points, glyph, actor, initial_positions)
    """
    points = vtk.vtkPoints()
    initial_positions = []

    for _ in range(n_particles):
        px = np.random.uniform(0, Lx)
        py = np.random.uniform(0, Ly)
        pz = 0.05  # Slightly above surface
        points.InsertNextPoint(px, py, pz)
        initial_positions.append([px, py, pz])

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(min(Lx, Ly) * 0.015)
    sphere.SetThetaResolution(8)
    sphere.SetPhiResolution(8)

    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetInputData(polyData)
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1.0, 0.9, 0.2)
    actor.GetProperty().SetOpacity(0.9)

    return points, glyph, actor, initial_positions


def vtk_3d_visualization(
    x: np.ndarray,
    y: np.ndarray,
    phi: np.ndarray,
    u: float,
    v: float,
    Gamma: float,
    Pe_x: float,
    Pe_y: float,
    phi_left: float,
    phi_right: float,
    phi_bottom: float,
    phi_top: float,
    title: str = "2D Convection-Diffusion Visualization",
):
    """
    Create an ambitious 3D visualization of the 2D convection-diffusion solution.

    Features (extending 1D version to 2D):
    - 3D surface with height representing scalar value
    - Color mapping for scalar field
    - Velocity glyph field showing flow direction
    - Contour lines on the surface
    - Animated particles showing convection-diffusion transport
    - Boundary condition annotations on all four edges
    - Governing equation display and Peclet numbers
    - Enhanced lighting

    :param x: X-coordinate vector (numpy.ndarray)
    :param y: Y-coordinate vector (numpy.ndarray)
    :param phi: Scalar field (numpy.ndarray)
    :param u: x-velocity component (float)
    :param v: y-velocity component (float)
    :param Gamma: Diffusion coefficient (float)
    :param Pe_x: Peclet number in x-direction (float)
    :param Pe_y: Peclet number in y-direction (float)
    :param phi_left: Left boundary value (float)
    :param phi_right: Right boundary value (float)
    :param phi_bottom: Bottom boundary value (float)
    :param phi_top: Top boundary value (float)
    :param title: Window title (str)
    """
    Lx = x[-1] - x[0]
    Ly = y[-1] - y[0]
    height_scale = min(Lx, Ly) * 0.5

    # ========== Create 3D surface ==========
    surface_polydata = create_2d_surface(x, y, phi, height_scale=height_scale)
    color_map = create_color_map(phi)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(surface_polydata)
    mapper.SetLookupTable(color_map)
    mapper.SetScalarModeToUsePointData()
    mapper.SetScalarRange(np.min(phi), np.max(phi))

    surface_actor = vtk.vtkActor()
    surface_actor.SetMapper(mapper)
    surface_actor.GetProperty().SetInterpolationToPhong()
    surface_actor.GetProperty().SetSpecular(0.4)
    surface_actor.GetProperty().SetSpecularPower(30)

    # ========== Contour lines on surface ==========
    contour_actor = create_contour_lines(surface_polydata, phi, n_contours=12)

    # ========== Velocity glyph field ==========
    velocity_glyphs = create_velocity_glyph_field_2d(x, y, u, v, n_glyphs=6)

    # ========== Particle system for animation ==========
    particle_points, particle_glyph, particle_actor, initial_positions = \
        create_particle_system_2d(Lx, Ly, n_particles=N_PARTICLES_2D)

    # ========== Boundary annotations ==========
    offset = 0.15 * max(Lx, Ly)
    annotations = [
        create_boundary_annotation_2d(f"phi={phi_left}", (-offset, Ly/2, 0), (0.3, 0.7, 1.0)),
        create_boundary_annotation_2d(f"phi={phi_right}", (Lx + offset, Ly/2, 0), (1.0, 0.5, 0.4)),
        create_boundary_annotation_2d(f"phi={phi_bottom}", (Lx/2, -offset, 0), (0.4, 0.9, 0.4)),
        create_boundary_annotation_2d(f"phi={phi_top}", (Lx/2, Ly + offset, 0), (1.0, 0.8, 0.3)),
    ]

    # ========== Physics label ==========
    physics_label = create_boundary_annotation_2d(
        f"Pe_x={Pe_x:.1f}  Pe_y={Pe_y:.1f}\nu={u} m/s  v={v} m/s\nGamma={Gamma}",
        (Lx/2, Ly + offset * 2, height_scale * 0.8),
        color=(0.9, 0.9, 0.4),
    )

    # ========== Scalar bar ==========
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Phi")
    scalar_bar.SetNumberOfLabels(6)
    scalar_bar.SetWidth(0.06)
    scalar_bar.SetHeight(0.35)
    scalar_bar.SetPosition(0.92, 0.32)

    # ========== Title annotation ==========
    title_actor = vtk.vtkTextActor()
    title_actor.SetInput(title)
    title_actor.GetTextProperty().SetFontSize(20)
    title_actor.GetTextProperty().SetColor(1, 1, 1)
    title_actor.GetTextProperty().SetBold(True)
    title_actor.SetPosition(15, 15)

    # ========== Governing equation annotation ==========
    equation_actor = vtk.vtkTextActor()
    equation_actor.SetInput("Governing Equation:\nrho*u*dphi/dx + rho*v*dphi/dy = Gamma*(d2phi/dx2 + d2phi/dy2)")
    equation_actor.GetTextProperty().SetFontSize(12)
    equation_actor.GetTextProperty().SetColor(0.8, 0.8, 0.8)
    equation_actor.GetTextProperty().SetFontFamilyToCourier()
    equation_actor.SetPosition(15, 80)

    # ========== Renderer ==========
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(contour_actor)
    renderer.AddActor(velocity_glyphs)
    renderer.AddActor(particle_actor)
    for ann in annotations:
        renderer.AddActor(ann)
    renderer.AddActor(physics_label)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)
    renderer.AddActor2D(equation_actor)

    # Gradient background
    renderer.SetBackground(0.08, 0.08, 0.12)
    renderer.SetBackground2(0.02, 0.02, 0.04)
    renderer.GradientBackgroundOn()

    # Enhanced lighting
    light1 = vtk.vtkLight()
    light1.SetPosition(Lx/2, Ly/2, Lx * 2)
    light1.SetFocalPoint(Lx/2, Ly/2, 0)
    light1.SetIntensity(0.8)
    renderer.AddLight(light1)

    light2 = vtk.vtkLight()
    light2.SetPosition(Lx * 2, Ly/2, Lx)
    light2.SetFocalPoint(Lx/2, Ly/2, 0)
    light2.SetIntensity(0.4)
    renderer.AddLight(light2)

    # ========== Render window ==========
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1400, 900)
    renderWindow.SetWindowName(title)

    # ========== Interactor ==========
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interact_style)

    # ========== Axes widget ==========
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0, 0, 0.18, 0.18)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # ========== Camera position ==========
    camera = renderer.GetActiveCamera()
    camera.SetPosition(Lx * 2.5, -Ly * 1.5, max(Lx, Ly) * 2)
    camera.SetFocalPoint(Lx/2, Ly/2, height_scale * 0.3)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    camera.Zoom(1.1)

    # ========== Animation callback for particles ==========
    class ParticleAnimator2D:
        def __init__(self, points, glyph, initial_pos, Lx, Ly, u, v, dt=0.02):
            self.points = points
            self.glyph = glyph
            self.positions = [list(p) for p in initial_pos]
            self.Lx = Lx
            self.Ly = Ly
            self.u = u
            self.v = v
            self.dt = dt
            self.time = 0

        def update(self, obj, event):
            self.time += self.dt
            for i, pos in enumerate(self.positions):
                # Convection: particle advected by flow velocity
                pos[0] += self.u * self.dt * CONVECTION_SCALE_FACTOR
                pos[1] += self.v * self.dt * CONVECTION_SCALE_FACTOR

                # Diffusion: random walk (using dimensionless scale for visualization)
                # Note: This is a visualization aid, not physically accurate diffusion
                diffusion_scale = min(self.Lx, self.Ly) * 0.01 * DIFFUSION_SCALE_FACTOR
                pos[0] += np.random.normal(0, diffusion_scale)
                pos[1] += np.random.normal(0, diffusion_scale)

                # Respawn at inlet when hitting outlet (consistent with Dirichlet BCs)
                # Particles represent tracer injection at boundaries
                if pos[0] > self.Lx:
                    pos[0] = 0  # Respawn at left inlet
                    pos[1] = np.random.uniform(0, self.Ly)
                if pos[0] < 0:
                    pos[0] = self.Lx  # Respawn at right
                    pos[1] = np.random.uniform(0, self.Ly)
                if pos[1] > self.Ly:
                    pos[1] = 0  # Respawn at bottom
                    pos[0] = np.random.uniform(0, self.Lx)
                if pos[1] < 0:
                    pos[1] = self.Ly  # Respawn at top
                    pos[0] = np.random.uniform(0, self.Lx)

                self.points.SetPoint(i, pos[0], pos[1], pos[2])

            self.points.Modified()
            self.glyph.Update()
            renderWindow.Render()

    animator = ParticleAnimator2D(
        particle_points, particle_glyph, initial_positions, Lx, Ly, u, v
    )

    # Create timer for animation
    renderWindowInteractor.AddObserver("TimerEvent", animator.update)
    timer_id = renderWindowInteractor.CreateRepeatingTimer(ANIMATION_TIMER_MS)

    renderWindow.Render()
    renderWindowInteractor.Start()

    # Cleanup timer
    renderWindowInteractor.DestroyTimer(timer_id)


def main():
    """
    Main function demonstrating 2D convection-diffusion solver.

    This is the logical extension of the 1D solver to two dimensions.
    """
    # Physical parameters (extending 1D case to 2D)
    Lx = 1.0        # Domain length in x (m)
    Ly = 1.0        # Domain length in y (m)
    nx = 50         # Number of grid divisions in x
    ny = 50         # Number of grid divisions in y
    rho = 1.0       # Fluid density (kg/m³)
    Gamma = 0.1     # Diffusion coefficient (kg/m·s)
    u = 2.0         # x-velocity (m/s)
    v = 1.0         # y-velocity (m/s)

    # Boundary conditions (Dirichlet on all edges)
    phi_left = 0.0      # x = 0
    phi_right = 1.0     # x = Lx
    phi_bottom = 0.0    # y = 0
    phi_top = 1.0       # y = Ly

    Rtol = 1e-6     # Convergence tolerance

    # Print header
    print("=" * 70)
    print("2D Steady-State Convection-Diffusion Solver")
    print("(Logical Extension of 1D Solver)")
    print("=" * 70)
    print(f"Domain: {Lx}m x {Ly}m, Grid: ({nx+1}) x ({ny+1})")
    print(f"Velocity: u = {u} m/s, v = {v} m/s")
    print(f"Diffusion coefficient: Gamma = {Gamma} kg/m·s")
    print("-" * 70)
    print("Boundary Conditions:")
    print(f"  Left (x=0):   phi = {phi_left}")
    print(f"  Right (x=Lx): phi = {phi_right}")
    print(f"  Bottom (y=0): phi = {phi_bottom}")
    print(f"  Top (y=Ly):   phi = {phi_top}")
    print("-" * 70)
    print("Solving...")

    # Solve the 2D convection-diffusion equation
    x, y, phi, iterations, Pe_x, Pe_y = solve_convection_diffusion_2d(
        nx, ny, Lx, Ly, Gamma, rho, u, v,
        phi_left, phi_right, phi_bottom, phi_top, Rtol
    )

    # Print results summary
    print("-" * 70)
    print(f"Converged in {iterations} iterations")
    print(f"Peclet numbers: Pe_x = {Pe_x:.2f}, Pe_y = {Pe_y:.2f}")
    print(f"Phi range: [{np.min(phi):.4f}, {np.max(phi):.4f}]")
    print(f"Phi at center: {phi[ny//2, nx//2]:.4f}")
    print("=" * 70)

    # Plot results using matplotlib
    plot_solution_2d(x, y, phi, u, v, Pe_x, Pe_y)

    # 3D VTK visualization with all ambitious features
    vtk_3d_visualization(
        x, y, phi, u, v, Gamma, Pe_x, Pe_y,
        phi_left, phi_right, phi_bottom, phi_top
    )


if __name__ == "__main__":
    main()
