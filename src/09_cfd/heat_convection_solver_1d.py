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
    phi: np.ndarray,
    phi_0: float,
    phi_L: float,
    Pe: float,
    u: float,
    scheme_name: str,
    title: str = "1D Convection-Diffusion",
):
    """
    Create an advanced 3D visualization of the convection-diffusion solution.

    Features:
    - 3D cylindrical representation with scalar field coloring
    - End caps for the domain
    - Flow direction arrows
    - Boundary condition annotations
    - Peclet number display

    :param x: X-coordinate vector (numpy.ndarray)
    :param phi: Scalar field values (numpy.ndarray)
    :param phi_0: Left boundary value (float)
    :param phi_L: Right boundary value (float)
    :param Pe: Peclet number (float)
    :param u: Flow velocity (float)
    :param scheme_name: Name of the discretization scheme (str)
    :param title: Window title (str)
    """
    L = x[-1] - x[0]
    radius = 0.08 * L

    # Create 3D cylinder rod
    rod_polydata, _ = create_cylinder_rod(x, phi, radius=radius)

    # Create end caps
    left_cap = create_end_cap(x[0], phi[0], radius=radius)
    right_cap = create_end_cap(x[-1], phi[-1], radius=radius)

    # Combine rod and caps
    append_filter = vtk.vtkAppendPolyData()
    append_filter.AddInputData(rod_polydata)
    append_filter.AddInputData(left_cap)
    append_filter.AddInputData(right_cap)
    append_filter.Update()

    # Color map
    color_map = create_color_map(phi)

    # Rod mapper and actor
    rod_mapper = vtk.vtkPolyDataMapper()
    rod_mapper.SetInputData(append_filter.GetOutput())
    rod_mapper.SetLookupTable(color_map)
    rod_mapper.SetScalarModeToUsePointData()
    rod_mapper.SetScalarRange(np.min(phi), np.max(phi))

    rod_actor = vtk.vtkActor()
    rod_actor.SetMapper(rod_mapper)
    rod_actor.GetProperty().SetInterpolationToPhong()
    rod_actor.GetProperty().SetSpecular(0.3)
    rod_actor.GetProperty().SetSpecularPower(20)

    # Create flow direction arrow
    flow_direction = 1 if u >= 0 else -1
    arrow_y_offset = radius * 2.5
    flow_arrow = create_flow_arrow(
        start_pos=(L / 2 - 0.1 * L * flow_direction, arrow_y_offset, 0),
        direction=(flow_direction, 0, 0),
        scale=0.2 * L,
    )

    # Create boundary annotations
    annotation_offset = radius * 3
    left_text = f"Inlet\nφ = {phi_0}"
    right_text = f"Outlet\nφ = {phi_L}"

    left_annotation = create_boundary_annotation(
        left_text, (x[0] - 0.12 * L, 0, annotation_offset), color=(0.3, 0.6, 1.0)
    )
    right_annotation = create_boundary_annotation(
        right_text, (x[-1] + 0.12 * L, 0, annotation_offset), color=(1.0, 0.4, 0.4)
    )

    # Flow and Peclet number labels
    flow_label = create_boundary_annotation(
        f"Flow → (u = {u} m/s)" if u >= 0 else f"← Flow (u = {u} m/s)",
        (L / 2, arrow_y_offset + radius, 0),
        color=(0.2, 0.7, 1.0),
    )
    peclet_label = create_boundary_annotation(
        f"Pe = {Pe:.1f}\n({scheme_name})",
        (L / 2, -radius * 2.5, 0),
        color=(0.9, 0.9, 0.3),
    )

    # Scalar bar
    scalar_bar = vtk.vtkScalarBarActor()
    scalar_bar.SetLookupTable(color_map)
    scalar_bar.SetTitle("Phi (Scalar)")
    scalar_bar.SetNumberOfLabels(5)
    scalar_bar.SetWidth(0.08)
    scalar_bar.SetHeight(0.4)
    scalar_bar.SetPosition(0.9, 0.3)

    # Title annotation
    title_actor = vtk.vtkTextActor()
    title_actor.SetInput(title)
    title_actor.GetTextProperty().SetFontSize(18)
    title_actor.GetTextProperty().SetColor(1, 1, 1)
    title_actor.GetTextProperty().SetBold(True)
    title_actor.SetPosition(10, 10)

    # Renderer setup
    renderer = vtk.vtkRenderer()
    renderer.AddActor(rod_actor)
    renderer.AddActor(flow_arrow)
    renderer.AddActor(left_annotation)
    renderer.AddActor(right_annotation)
    renderer.AddActor(flow_label)
    renderer.AddActor(peclet_label)
    renderer.AddActor2D(scalar_bar)
    renderer.AddActor2D(title_actor)

    # Gradient background
    renderer.SetBackground(0.1, 0.1, 0.15)
    renderer.SetBackground2(0.02, 0.02, 0.05)
    renderer.GradientBackgroundOn()

    # Lighting
    light = vtk.vtkLight()
    light.SetPosition(L / 2, L, L)
    light.SetFocalPoint(L / 2, 0, 0)
    light.SetIntensity(1.0)
    renderer.AddLight(light)

    # Render window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(1200, 700)
    renderWindow.SetWindowName(title)

    # Interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    interact_style = vtk.vtkInteractorStyleTrackballCamera()
    renderWindowInteractor.SetInteractorStyle(interact_style)

    # Axes widget
    axes = vtk.vtkAxesActor()
    axes_widget = vtk.vtkOrientationMarkerWidget()
    axes_widget.SetOrientationMarker(axes)
    axes_widget.SetInteractor(renderWindowInteractor)
    axes_widget.SetViewport(0, 0, 0.2, 0.2)
    axes_widget.SetEnabled(1)
    axes_widget.InteractiveOff()

    # Camera position
    camera = renderer.GetActiveCamera()
    camera.SetPosition(L / 2, -L * 0.8, L * 0.6)
    camera.SetFocalPoint(L / 2, 0, 0)
    camera.SetViewUp(0, 0, 1)
    renderer.ResetCamera()
    camera.Zoom(1.2)

    renderWindow.Render()
    renderWindowInteractor.Start()


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

    # 3D VTK visualization
    vtk_3d_visualization(x, phi_numerical, phi_0, phi_L, Pe, u, scheme_name)

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
