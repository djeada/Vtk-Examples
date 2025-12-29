"""
1D Steady-State Convection-Diffusion Solver

This module provides a solver for the 1D steady-state convection-diffusion equation,
which models the transport of a scalar quantity (e.g., temperature, concentration)
under the combined effects of convection (advection by fluid flow) and diffusion.

Workflow:
1. Initialization: Set up the computational grid, material properties (diffusivity Gamma),
   and flow characteristics (velocity u, density rho).
2. Matrix Assembly: Construct the tridiagonal matrix using either Upwind Difference Scheme
   (UDS) or Central Difference Scheme (CDS) for the convection term.
3. Boundary Conditions: Apply Dirichlet boundary conditions at both ends.
4. Solver: Use the Tridiagonal Matrix Algorithm (TDMA) to solve the linear system.
5. Post-Processing: Compare numerical results with the analytical solution for validation.

Physics:
The steady-state convection-diffusion equation in 1D is:
    d/dx(rho * u * phi) = d/dx(Gamma * dphi/dx)

where:
- phi: transported scalar (temperature, concentration, etc.)
- rho: fluid density (kg/m³)
- u: fluid velocity (m/s)
- Gamma: diffusion coefficient (kg/m·s or W/m·K for heat)

For constant properties, this simplifies to:
    rho * u * dphi/dx = Gamma * d²phi/dx²

The Peclet number Pe = (rho * u * L) / Gamma characterizes the relative importance
of convection to diffusion. High Pe indicates convection-dominated flow.

Discretization Schemes:
1. Upwind Difference Scheme (UDS):
   - Uses upwind approximation for convection: stable for all Pe
   - First-order accurate, introduces numerical diffusion
   - Coefficients: a_P = a_E + a_W + (F_e - F_w)
     where a_E = D_e, a_W = D_w + F_w (for positive flow)

2. Central Difference Scheme (CDS):
   - Uses central difference for convection: second-order accurate
   - Can become unstable for |Pe| > 2 (cell Peclet number)
   - Coefficients: a_E = D_e - F_e/2, a_W = D_w + F_w/2

Analytical Solution:
For constant properties with boundary conditions phi(0) = phi_0, phi(L) = phi_L:
    phi(x) = phi_0 + (phi_L - phi_0) * (exp(Pe*x/L) - 1) / (exp(Pe) - 1)
"""

import matplotlib.pyplot as plt
import numpy as np
import vtk

# ========== Animation and Visualization Constants ==========
ANIMATION_TIMER_MS = 50  # Animation update interval in milliseconds
PARTICLE_RADIUS_FACTOR = 0.8  # Particles confined to this fraction of rod radius
CONVECTION_SCALE_FACTOR = 0.5  # Scale factor for particle convection animation
DIFFUSION_SCALE_FACTOR = 0.3  # Scale factor for particle diffusion animation
LATERAL_DIFFUSION_STD = 0.002  # Standard deviation for lateral particle motion


def tdma_solve(
    a: np.ndarray, ad: np.ndarray, au: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """
    Solves a tridiagonal matrix equation using the Thomas algorithm (TDMA).

    The system has the form: ad[i]*x[i-1] + a[i]*x[i] + au[i]*x[i+1] = b[i]

    :param a: Main diagonal elements (numpy.ndarray)
    :param ad: Sub-diagonal (lower) elements (numpy.ndarray)
    :param au: Super-diagonal (upper) elements (numpy.ndarray)
    :param b: Right-hand side vector (numpy.ndarray)
    :return: Solution vector (numpy.ndarray)
    """
    n = len(b)
    # Make copies to avoid modifying original arrays
    a = a.copy()
    b = b.copy()
    x = np.zeros(n)

    # Forward elimination
    for i in range(1, n):
        f = ad[i] / a[i - 1]
        a[i] -= f * au[i - 1]
        b[i] -= f * b[i - 1]

    # Back substitution
    x[-1] = b[-1] / a[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (b[i] - au[i] * x[i + 1]) / a[i]

    return x


def setup_matrix(
    n: int, dx: float, Gamma: float, F: float, use_uds: bool
) -> tuple:
    """
    Sets up the coefficient matrix for the TDMA solver based on the chosen scheme.

    The discretized convection-diffusion equation for interior nodes is:
        a_W * phi[i-1] + a_P * phi[i] + a_E * phi[i+1] = 0

    where:
    - D = Gamma/dx (diffusion conductance)
    - F = rho*u (convective mass flux)

    For UDS (Upwind Difference Scheme) with positive flow (F > 0):
        a_W = D + F,  a_E = D,  a_P = a_W + a_E = 2D + F

    For CDS (Central Difference Scheme):
        a_W = D + F/2,  a_E = D - F/2,  a_P = 2D

    Note: CDS can become unstable when cell Peclet number |Pe_cell| = |F*dx/Gamma| > 2

    :param n: Number of grid divisions (int)
    :param dx: Grid spacing (float)
    :param Gamma: Diffusion coefficient (float, kg/m·s or W/m·K)
    :param F: Convective mass flux rho*u (float, kg/m²·s)
    :param use_uds: True for Upwind scheme, False for Central scheme (bool)
    :return: Tuple of (A, Ad, Au) - main, lower, and upper diagonals
    """
    num_nodes = n + 1
    D = Gamma / dx  # Diffusion conductance

    A = np.zeros(num_nodes)   # Main diagonal (a_P)
    Au = np.zeros(num_nodes)  # Upper diagonal (coefficient for phi[i+1])
    Ad = np.zeros(num_nodes)  # Lower diagonal (coefficient for phi[i-1])

    if use_uds:
        # Upwind Difference Scheme (first-order, unconditionally stable)
        # For positive flow: a_W = D + F, a_E = D
        a_W = D + max(F, 0)
        a_E = D + max(-F, 0)
        a_P = a_W + a_E

        A.fill(a_P)
        Au.fill(-a_E)   # Upper: -a_E (coefficient for phi[i+1])
        Ad.fill(-a_W)   # Lower: -a_W (coefficient for phi[i-1])
    else:
        # Central Difference Scheme (second-order, conditionally stable)
        a_W = D + F / 2
        a_E = D - F / 2
        a_P = a_W + a_E  # = 2D

        A.fill(a_P)
        Au.fill(-a_E)
        Ad.fill(-a_W)

    return A, Ad, Au


def analytical_solution(
    x: np.ndarray, L: float, Pe: float, phi_0: float, phi_L: float
) -> np.ndarray:
    """
    Computes the analytical solution for the 1D steady-state convection-diffusion equation.

    For the equation: rho*u*dphi/dx = Gamma*d²phi/dx²
    With boundary conditions: phi(0) = phi_0, phi(L) = phi_L

    The analytical solution is:
        phi(x) = phi_0 + (phi_L - phi_0) * (exp(Pe*x/L) - 1) / (exp(Pe) - 1)

    For Pe -> 0 (pure diffusion), this reduces to linear profile.
    For large Pe (convection-dominated), the profile becomes exponential.

    :param x: Position array (numpy.ndarray)
    :param L: Domain length (float)
    :param Pe: Peclet number = rho*u*L/Gamma (float)
    :param phi_0: Value at x=0 (float)
    :param phi_L: Value at x=L (float)
    :return: Solution array (numpy.ndarray)
    """
    if abs(Pe) < 1e-10:
        # Pure diffusion case: linear profile
        return phi_0 + (phi_L - phi_0) * x / L
    else:
        # General case: exponential profile
        return phi_0 + (phi_L - phi_0) * (np.exp(Pe * x / L) - 1) / (np.exp(Pe) - 1)


def plot_solution(
    x: np.ndarray,
    phi_numerical: np.ndarray,
    phi_analytical: np.ndarray,
    Pe: float,
    scheme_name: str,
):
    """
    Plots the numerical and analytical solutions.

    :param x: Grid positions (numpy.ndarray)
    :param phi_numerical: Numerical solution (numpy.ndarray)
    :param phi_analytical: Analytical solution (numpy.ndarray)
    :param Pe: Peclet number (float)
    :param scheme_name: Name of the discretization scheme (str)
    """
    plt.figure(figsize=(10, 6))
    plt.plot(x, phi_analytical, "r-", linewidth=2, label="Analytical Solution")
    plt.plot(
        x, phi_numerical, "b--", linewidth=2, marker="o", markersize=3,
        markevery=max(1, len(x) // 20), label=f"Numerical ({scheme_name})"
    )
    plt.xlabel("Position x (m)", fontsize=12)
    plt.ylabel("Phi (Transported Scalar)", fontsize=12)
    plt.title(f"1D Convection-Diffusion: Pe = {Pe:.1f}", fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
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


def create_cylinder_rod(
    x: np.ndarray, phi: np.ndarray, radius: float = 0.05, resolution: int = 20
) -> tuple:
    """
    Create a 3D cylindrical rod representation with scalar field coloring.

    :param x: X-coordinate vector (numpy.ndarray)
    :param phi: Scalar field values (numpy.ndarray)
    :param radius: Radius of the cylinder (float)
    :param resolution: Number of segments around the cylinder circumference (int)
    :return: Tuple of (polydata, phi_values) for VTK
    """
    points = vtk.vtkPoints()
    phi_values = vtk.vtkFloatArray()
    phi_values.SetName("Phi")
    cells = vtk.vtkCellArray()

    # Create cylinder surface points for each x position
    for i, (xi, phii) in enumerate(zip(x, phi)):
        for j in range(resolution):
            theta = 2.0 * np.pi * j / resolution
            y = radius * np.cos(theta)
            z = radius * np.sin(theta)
            points.InsertNextPoint(xi, y, z)
            phi_values.InsertNextValue(phii)

    # Create quad cells connecting adjacent rings
    for i in range(len(x) - 1):
        for j in range(resolution):
            quad = vtk.vtkQuad()
            p0 = i * resolution + j
            p1 = i * resolution + (j + 1) % resolution
            p2 = (i + 1) * resolution + (j + 1) % resolution
            p3 = (i + 1) * resolution + j

            quad.GetPointIds().SetId(0, p0)
            quad.GetPointIds().SetId(1, p1)
            quad.GetPointIds().SetId(2, p2)
            quad.GetPointIds().SetId(3, p3)
            cells.InsertNextCell(quad)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetPolys(cells)
    polyData.GetPointData().SetScalars(phi_values)

    return polyData, phi_values


def create_end_cap(
    x_pos: float, phi_val: float, radius: float = 0.05, resolution: int = 20
) -> vtk.vtkPolyData:
    """
    Create an end cap (disk) for the cylinder.

    :param x_pos: X position of the cap (float)
    :param phi_val: Scalar value for the cap (float)
    :param radius: Radius of the cap (float)
    :param resolution: Number of segments (int)
    :return: vtkPolyData for the end cap
    """
    disk = vtk.vtkDiskSource()
    disk.SetInnerRadius(0)
    disk.SetOuterRadius(radius)
    disk.SetRadialResolution(1)
    disk.SetCircumferentialResolution(resolution)
    disk.Update()

    transform = vtk.vtkTransform()
    transform.Translate(x_pos, 0, 0)
    transform.RotateY(90)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(disk.GetOutput())
    transformFilter.Update()

    output = transformFilter.GetOutput()
    phi_array = vtk.vtkFloatArray()
    phi_array.SetName("Phi")
    for _ in range(output.GetNumberOfPoints()):
        phi_array.InsertNextValue(phi_val)
    output.GetPointData().SetScalars(phi_array)

    return output


def create_flow_arrow(
    start_pos: tuple, direction: tuple, scale: float = 0.15
) -> vtk.vtkActor:
    """
    Create an arrow to represent flow direction.

    :param start_pos: Starting position (x, y, z)
    :param direction: Direction vector (dx, dy, dz)
    :param scale: Scale of the arrow (float)
    :return: vtkActor for the arrow
    """
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(20)
    arrow.SetShaftResolution(20)

    transform = vtk.vtkTransform()
    transform.Translate(start_pos)

    # Rotate arrow if flow is in negative direction
    if direction[0] < 0:
        transform.RotateZ(180)

    transform.Scale(scale, scale * 0.3, scale * 0.3)

    transformFilter = vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputConnection(arrow.GetOutputPort())

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(transformFilter.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.2, 0.6, 1.0)  # Blue color for flow

    return actor


def create_velocity_glyph_field(
    L: float, radius: float, u: float, n_glyphs: int = 8
) -> vtk.vtkActor:
    """
    Create a field of velocity arrows showing flow direction along the domain.

    :param L: Domain length (float)
    :param radius: Rod radius (float)
    :param u: Flow velocity (float)
    :param n_glyphs: Number of glyph arrows (int)
    :return: vtkActor with velocity glyphs
    """
    # Create points for glyph positions
    points = vtk.vtkPoints()
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetName("Velocity")

    flow_dir = 1 if u >= 0 else -1
    y_offset = radius * 2.0

    for i in range(n_glyphs):
        x_pos = L * (i + 0.5) / n_glyphs
        points.InsertNextPoint(x_pos, y_offset, 0)
        vectors.InsertNextTuple3(flow_dir * abs(u), 0, 0)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.GetPointData().SetVectors(vectors)

    # Create arrow glyph
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetShaftResolution(16)

    glyph = vtk.vtkGlyph3D()
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputData(polyData)
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    # Use a default scale when velocity is zero or very small
    base_scale = 0.08 * L
    glyph.SetScaleFactor(base_scale / max(abs(u), 0.1))
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.3, 0.7, 1.0)
    actor.GetProperty().SetOpacity(0.8)

    return actor


def create_profile_curve_3d(
    x: np.ndarray, phi: np.ndarray, y_offset: float, color: tuple, line_width: float = 3.0
) -> vtk.vtkActor:
    """
    Create a 3D line showing the phi profile floating above/below the rod.

    :param x: X-coordinate vector (numpy.ndarray)
    :param phi: Scalar field values (numpy.ndarray)
    :param y_offset: Y position offset (float)
    :param color: Line color (r, g, b)
    :param line_width: Width of the line (float)
    :return: vtkActor for the profile curve
    """
    points = vtk.vtkPoints()
    phi_min, phi_max = np.min(phi), np.max(phi)
    phi_range = phi_max - phi_min if phi_max > phi_min else 1.0
    L = x[-1] - x[0]
    scale = 0.3 * L / phi_range

    for i, (xi, phii) in enumerate(zip(x, phi)):
        z = (phii - phi_min) * scale
        points.InsertNextPoint(xi, y_offset, z)

    # Create polyline
    polyline = vtk.vtkPolyLine()
    polyline.GetPointIds().SetNumberOfIds(len(x))
    for i in range(len(x)):
        polyline.GetPointIds().SetId(i, i)

    cells = vtk.vtkCellArray()
    cells.InsertNextCell(polyline)

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)
    polyData.SetLines(cells)

    # Create tube filter for nicer appearance
    tube = vtk.vtkTubeFilter()
    tube.SetInputData(polyData)
    tube.SetRadius(0.005 * L)
    tube.SetNumberOfSides(12)
    tube.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetSpecular(0.5)
    actor.GetProperty().SetSpecularPower(30)

    return actor


def create_contour_bands(
    rod_polydata: vtk.vtkPolyData, phi: np.ndarray, n_bands: int = 10
) -> vtk.vtkActor:
    """
    Create contour bands on the rod surface.

    :param rod_polydata: The rod polydata (vtkPolyData)
    :param phi: Scalar field values (numpy.ndarray)
    :param n_bands: Number of contour bands (int)
    :return: vtkActor with contour lines
    """
    phi_min, phi_max = np.min(phi), np.max(phi)

    # Create banded contour filter
    contour = vtk.vtkBandedPolyDataContourFilter()
    contour.SetInputData(rod_polydata)
    contour.SetNumberOfContours(n_bands)
    contour.GenerateValues(n_bands, phi_min, phi_max)
    contour.SetGenerateContourEdges(True)
    contour.Update()

    # Create mapper for contour edges
    edge_mapper = vtk.vtkPolyDataMapper()
    edge_mapper.SetInputConnection(contour.GetOutputPort(1))
    edge_mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(edge_mapper)
    actor.GetProperty().SetColor(0.1, 0.1, 0.1)
    actor.GetProperty().SetLineWidth(2.0)

    return actor


def create_particle_system(
    L: float, radius: float, n_particles: int = 20
) -> tuple:
    """
    Create a particle system for animation.

    :param L: Domain length (float)
    :param radius: Rod radius (float)
    :param n_particles: Number of particles (int)
    :return: Tuple of (points, actor, initial_positions)
    """
    points = vtk.vtkPoints()
    initial_positions = []

    for i in range(n_particles):
        x = np.random.uniform(0, L)
        y = np.random.uniform(-radius * PARTICLE_RADIUS_FACTOR, radius * PARTICLE_RADIUS_FACTOR)
        z = np.random.uniform(-radius * PARTICLE_RADIUS_FACTOR, radius * PARTICLE_RADIUS_FACTOR)
        # Ensure inside cylinder (r^2 < (radius * factor)^2)
        max_r_squared = (radius * PARTICLE_RADIUS_FACTOR) ** 2
        while y * y + z * z > max_r_squared:
            y = np.random.uniform(-radius * PARTICLE_RADIUS_FACTOR, radius * PARTICLE_RADIUS_FACTOR)
            z = np.random.uniform(-radius * PARTICLE_RADIUS_FACTOR, radius * PARTICLE_RADIUS_FACTOR)
        points.InsertNextPoint(x, y, z)
        initial_positions.append([x, y, z])

    polyData = vtk.vtkPolyData()
    polyData.SetPoints(points)

    # Create sphere glyphs for particles
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius * 0.08)
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
    actor.GetProperty().SetColor(1.0, 0.9, 0.2)  # Yellow particles
    actor.GetProperty().SetOpacity(0.9)

    return points, glyph, actor, initial_positions


def create_equation_annotation(L: float, radius: float) -> vtk.vtkTextActor:
    """
    Create a text annotation showing the governing equation.

    :param L: Domain length (float)
    :param radius: Rod radius (float)
    :return: vtkTextActor
    """
    text_actor = vtk.vtkTextActor()
    text_actor.SetInput("Governing Equation:\nrho*u*(dphi/dx) = Gamma*(d2phi/dx2)")
    text_actor.GetTextProperty().SetFontSize(14)
    text_actor.GetTextProperty().SetColor(0.8, 0.8, 0.8)
    text_actor.GetTextProperty().SetFontFamilyToCourier()
    text_actor.SetPosition(10, 100)

    return text_actor


def create_legend_actor(
    numerical_color: tuple, analytical_color: tuple
) -> vtk.vtkLegendBoxActor:
    """
    Create a legend showing numerical vs analytical curves.

    :param numerical_color: Color for numerical solution (r, g, b)
    :param analytical_color: Color for analytical solution (r, g, b)
    :return: vtkLegendBoxActor
    """
    legend = vtk.vtkLegendBoxActor()
    legend.SetNumberOfEntries(2)

    # Create small colored boxes for legend
    sphere = vtk.vtkSphereSource()
    sphere.Update()

    legend.SetEntry(0, sphere.GetOutput(), "Numerical", numerical_color)
    legend.SetEntry(1, sphere.GetOutput(), "Analytical", analytical_color)

    legend.SetPosition(0.02, 0.85)
    legend.SetWidth(0.15)
    legend.SetHeight(0.12)
    legend.GetProperty().SetColor(0.9, 0.9, 0.9)

    return legend


def create_boundary_annotation(
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
    text_actor.GetTextProperty().SetFontSize(14)
    text_actor.GetTextProperty().SetColor(color)
    text_actor.GetTextProperty().SetJustificationToCentered()
    text_actor.GetTextProperty().SetBold(True)

    return text_actor


def vtk_3d_visualization(
    x: np.ndarray,
    phi_numerical: np.ndarray,
    phi_analytical: np.ndarray,
    phi_0: float,
    phi_L: float,
    Pe: float,
    u: float,
    Gamma: float,
    scheme_name: str,
    title: str = "1D Convection-Diffusion Visualization",
):
    """
    Create an ambitious 3D visualization of the convection-diffusion solution.

    Features:
    - 3D cylindrical rod with scalar field coloring and contour bands
    - Velocity glyph field showing flow direction
    - 3D profile curves for numerical and analytical solutions
    - Animated particles showing convection-diffusion transport
    - Boundary condition annotations
    - Governing equation display
    - Legend comparing numerical vs analytical
    - Split viewport with profile comparison

    :param x: X-coordinate vector (numpy.ndarray)
    :param phi_numerical: Numerical solution (numpy.ndarray)
    :param phi_analytical: Analytical solution (numpy.ndarray)
    :param phi_0: Left boundary value (float)
    :param phi_L: Right boundary value (float)
    :param Pe: Peclet number (float)
    :param u: Flow velocity (float)
    :param Gamma: Diffusion coefficient (float)
    :param scheme_name: Name of the discretization scheme (str)
    :param title: Window title (str)
    """
    L = x[-1] - x[0]
    radius = 0.08 * L

    # ========== Create 3D rod with scalar coloring ==========
    rod_polydata, _ = create_cylinder_rod(x, phi_numerical, radius=radius)

    # Create end caps
    left_cap = create_end_cap(x[0], phi_numerical[0], radius=radius)
    right_cap = create_end_cap(x[-1], phi_numerical[-1], radius=radius)

    # Combine rod and caps
    append_filter = vtk.vtkAppendPolyData()
    append_filter.AddInputData(rod_polydata)
    append_filter.AddInputData(left_cap)
    append_filter.AddInputData(right_cap)
    append_filter.Update()

    # Color map
    phi_all = np.concatenate([phi_numerical, phi_analytical])
    color_map = create_color_map(phi_all)

    # Rod mapper and actor
    rod_mapper = vtk.vtkPolyDataMapper()
    rod_mapper.SetInputData(append_filter.GetOutput())
    rod_mapper.SetLookupTable(color_map)
    rod_mapper.SetScalarModeToUsePointData()
    rod_mapper.SetScalarRange(np.min(phi_all), np.max(phi_all))

    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(rod_mapper)
    rod_actor.GetProperty().SetInterpolationToPhong()
    rod_actor.GetProperty().SetSpecular(0.4)
    rod_actor.GetProperty().SetSpecularPower(30)

    # ========== Contour bands on rod ==========
    contour_actor = create_contour_bands(rod_polydata, phi_numerical, n_bands=8)

    # ========== Velocity glyph field ==========
    velocity_glyphs = create_velocity_glyph_field(L, radius, u, n_glyphs=10)

    # ========== 3D Profile curves ==========
    # Numerical solution curve (blue)
    numerical_curve = create_profile_curve_3d(
        x, phi_numerical, y_offset=-radius * 3.5, color=(0.2, 0.6, 1.0), line_width=4.0
    )
    # Analytical solution curve (red)
    analytical_curve = create_profile_curve_3d(
        x, phi_analytical, y_offset=-radius * 3.5, color=(1.0, 0.3, 0.3), line_width=4.0
    )

    # ========== Particle system for animation ==========
    particle_points, particle_glyph, particle_actor, initial_positions = \
        create_particle_system(L, radius, n_particles=30)

    # ========== Boundary annotations ==========
    annotation_offset = radius * 3
    left_text = f"Inlet\nphi = {phi_0}"
    right_text = f"Outlet\nphi = {phi_L}"

    left_annotation = create_boundary_annotation(
        left_text, (x[0] - 0.15 * L, 0, annotation_offset), color=(0.3, 0.7, 1.0)
    )
    right_annotation = create_boundary_annotation(
        right_text, (x[-1] + 0.15 * L, 0, annotation_offset), color=(1.0, 0.5, 0.4)
    )

    # ========== Flow and physics labels ==========
    flow_direction = 1 if u >= 0 else -1
    arrow_y_offset = radius * 3.5

    flow_label = create_boundary_annotation(
        f"Flow --> u = {u} m/s" if u >= 0 else f"<-- Flow u = {u} m/s",
        (L / 2, arrow_y_offset, 0),
        color=(0.3, 0.8, 1.0),
    )

    physics_label = create_boundary_annotation(
        f"Pe = {Pe:.1f} | Gamma = {Gamma}\n{scheme_name} Scheme",
        (L / 2, -radius * 5.5, 0),
        color=(0.9, 0.9, 0.4),
    )

    # Profile legend label
    profile_legend = create_boundary_annotation(
        "-- Numerical (blue)\n-- Analytical (red)",
        (L / 2, -radius * 7, 0),
        color=(0.8, 0.8, 0.8),
    )

    # ========== Scalar bar ==========
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Phi")
    scalar_bar.SetNumberOfLabels(6)
    scalar_bar.SetWidth(0.06)
    scalar_bar.SetHeight(0.35)
    scalar_bar.SetPosition(0.92, 0.32)
    scalar_bar.GetTitleTextProperty().SetFontSize(12)
    scalar_bar.GetLabelTextProperty().SetFontSize(10)

    # ========== Title and equation annotations ==========
    title_actor = vtk.vtkTextActor()
    title_actor.SetInput(title)
    title_actor.GetTextProperty().SetFontSize(22)
    title_actor.GetTextProperty().SetColor(1, 1, 1)
    title_actor.GetTextProperty().SetBold(True)
    title_actor.SetPosition(15, 15)

    equation_actor = create_equation_annotation(L, radius)

    # ========== Error metrics annotation ==========
    max_error = np.max(np.abs(phi_numerical - phi_analytical))
    error_text = vtk.vtkTextActor()
    error_text.SetInput(f"Max Error: {max_error:.4e}")
    error_text.GetTextProperty().SetFontSize(14)
    error_text.GetTextProperty().SetColor(0.9, 0.7, 0.3)
    error_text.SetPosition(15, 70)

    # ========== Main renderer (3D view) ==========
    renderer = vtk.vtkRenderer()
    renderer.SetViewport(0.0, 0.0, 1.0, 1.0)

    # Add all actors
    renderer.AddActor(rod_actor)
    renderer.AddActor(contour_actor)
    renderer.AddActor(velocity_glyphs)
    renderer.AddActor(numerical_curve)
    renderer.AddActor(analytical_curve)
    renderer.AddActor(particle_actor)
    renderer.AddActor(left_annotation)
    renderer.AddActor(right_annotation)
    renderer.AddActor(flow_label)
    renderer.AddActor(physics_label)
    renderer.AddActor(profile_legend)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)
    renderer.AddActor2D(equation_actor)
    renderer.AddActor2D(error_text)

    # Gradient background
    renderer.SetBackground(0.08, 0.08, 0.12)
    renderer.SetBackground2(0.02, 0.02, 0.04)
    renderer.GradientBackgroundOn()

    # Enhanced lighting
    light1 = vtk.vtkLight()
    light1.SetPosition(L / 2, L, L)
    light1.SetFocalPoint(L / 2, 0, 0)
    light1.SetIntensity(0.8)
    renderer.AddLight(light1)

    light2 = vtk.vtkLight()
    light2.SetPosition(L / 2, -L, L / 2)
    light2.SetFocalPoint(L / 2, 0, 0)
    light2.SetIntensity(0.4)
    renderer.AddLight(light2)

    # ========== Render window ==========
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1400, 800)
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
    camera.SetPosition(L / 2, -L * 1.0, L * 0.7)
    camera.SetFocalPoint(L / 2, 0, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    camera.Zoom(1.1)

    # ========== Animation callback for particles ==========
    class ParticleAnimator:
        def __init__(self, points, glyph, initial_pos, L, u, Gamma, dt=0.02):
            self.points = points
            self.glyph = glyph
            self.positions = [list(p) for p in initial_pos]
            self.L = L
            self.u = u
            self.Gamma = Gamma
            self.dt = dt
            self.time = 0

        def update(self, obj, event):
            self.time += self.dt
            for i, pos in enumerate(self.positions):
                # Convection: particle advected by flow velocity (scaled for visualization)
                pos[0] += self.u * self.dt * CONVECTION_SCALE_FACTOR
                # Diffusion: random walk in x-direction based on diffusion coefficient
                diffusion_std = np.sqrt(2 * self.Gamma * self.dt) * DIFFUSION_SCALE_FACTOR
                pos[0] += np.random.normal(0, diffusion_std)
                # Lateral diffusion (small random motion in y and z)
                pos[1] += np.random.normal(0, LATERAL_DIFFUSION_STD)
                pos[2] += np.random.normal(0, LATERAL_DIFFUSION_STD)

                # Periodic boundary
                if pos[0] > self.L:
                    pos[0] = 0
                if pos[0] < 0:
                    pos[0] = self.L

                self.points.SetPoint(i, pos[0], pos[1], pos[2])

            self.points.Modified()
            self.glyph.Update()
            renderWindow.Render()

    animator = ParticleAnimator(particle_points, particle_glyph, initial_positions, L, u, Gamma)

    # Create timer for animation
    renderWindowInteractor.AddObserver("TimerEvent", animator.update)
    timer_id = renderWindowInteractor.CreateRepeatingTimer(ANIMATION_TIMER_MS)

    renderWindow.Render()
    renderWindowInteractor.Start()

    # Cleanup timer
    renderWindowInteractor.DestroyTimer(timer_id)


def main():
    """
    Main function demonstrating 1D convection-diffusion solver.

    Solves the steady-state convection-diffusion equation with Dirichlet BCs
    and compares with the analytical solution.
    """
    # Physical parameters
    L = 1.0         # Domain length (m)
    n = 100         # Number of grid divisions
    rho = 1.0       # Fluid density (kg/m³)
    Gamma = 0.1     # Diffusion coefficient (kg/m·s)
    u = 2.5         # Flow velocity (m/s)
    phi_0 = 0.0     # Left boundary value
    phi_L = 1.0     # Right boundary value
    use_uds = True  # Use Upwind Difference Scheme

    # Derived quantities
    dx = L / n
    F = rho * u                  # Convective mass flux
    Pe = rho * u * L / Gamma     # Peclet number
    Pe_cell = abs(F * dx / Gamma)  # Cell Peclet number

    # Scheme name for display
    scheme_name = "UDS" if use_uds else "CDS"

    # Check CDS stability
    if not use_uds and Pe_cell > 2:
        print(f"WARNING: Cell Peclet number = {Pe_cell:.2f} > 2")
        print("CDS may produce oscillatory/unstable solution.")
        print("Consider using UDS or refining the mesh.")

    # Create position array
    x = np.linspace(0, L, n + 1)

    # Setup coefficient matrix
    A, Ad, Au = setup_matrix(n, dx, Gamma, F, use_uds)

    # Apply Dirichlet boundary conditions
    # Left BC: phi[0] = phi_0
    A[0] = 1.0
    Au[0] = 0.0
    Ad[0] = 0.0

    # Right BC: phi[n] = phi_L
    A[-1] = 1.0
    Au[-1] = 0.0
    Ad[-1] = 0.0

    # Right-hand side vector
    b = np.zeros(n + 1)
    b[0] = phi_0
    b[-1] = phi_L

    # Solve the system
    phi_numerical = tdma_solve(A, Ad, Au, b)

    # Compute analytical solution
    phi_analytical = analytical_solution(x, L, Pe, phi_0, phi_L)

    # Plot results
    plot_solution(x, phi_numerical, phi_analytical, Pe, scheme_name)

    # 3D VTK visualization with all the ambitious features
    vtk_3d_visualization(
        x, phi_numerical, phi_analytical, phi_0, phi_L, Pe, u, Gamma, scheme_name
    )

    # Print results summary
    print("=" * 65)
    print("1D Steady-State Convection-Diffusion Solver")
    print("=" * 65)
    print(f"Domain length: L = {L} m, Grid points: {n + 1}")
    print(f"Flow velocity: u = {u} m/s")
    print(f"Diffusion coefficient: Gamma = {Gamma} kg/m·s")
    print(f"Peclet number: Pe = {Pe:.2f}")
    print(f"Cell Peclet number: Pe_cell = {Pe_cell:.4f}")
    print(f"Scheme: {scheme_name}")
    print("-" * 65)
    print(f"Left BC: phi(0) = {phi_0}")
    print(f"Right BC: phi(L) = {phi_L}")
    print("-" * 65)
    print(f"Numerical solution at x=0:  {phi_numerical[0]:.6f}")
    print(f"Numerical solution at x=L:  {phi_numerical[-1]:.6f}")
    print(f"Analytical solution at x=0: {phi_analytical[0]:.6f}")
    print(f"Analytical solution at x=L: {phi_analytical[-1]:.6f}")
    print("-" * 65)
    max_error = np.max(np.abs(phi_numerical - phi_analytical))
    rms_error = np.sqrt(np.mean((phi_numerical - phi_analytical) ** 2))
    print(f"Maximum error: {max_error:.6e}")
    print(f"RMS error: {rms_error:.6e}")
    print("=" * 65)


if __name__ == "__main__":
    main()
