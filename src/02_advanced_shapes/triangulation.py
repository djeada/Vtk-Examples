"""
CFD Mesh Generation using Delaunay Triangulation with VTK

This module demonstrates mesh generation for computational fluid dynamics using
Delaunay triangulation. Mesh generation is a critical preprocessing step in CFD,
where the computational domain is divided into discrete elements for numerical
solution of governing equations.

Physical Context:
-----------------
We create a mesh for a thermal analysis problem: a plate with a circular hole
(representing a cooling channel or structural cutout). This is a common geometry
in heat transfer and stress analysis applications.

The mesh quality directly affects:
- Solution accuracy (element shape affects interpolation error)
- Convergence rate (poorly shaped elements slow iterative solvers)
- Numerical stability (extreme aspect ratios cause ill-conditioning)

Mathematical Background:
------------------------
Delaunay Triangulation:
For a set of points P in R², the Delaunay triangulation DT(P) satisfies the
"empty circumcircle" property: no point in P lies inside the circumcircle of
any triangle in DT(P).

Key properties:
1. Maximizes minimum angle: DT avoids "sliver" triangles, improving mesh quality
2. Unique (for general position): Given points, DT is uniquely determined
3. Dual of Voronoi diagram: Each Delaunay vertex corresponds to a Voronoi cell

For 3D, Delaunay tetrahedralization extends these concepts to tetrahedra and
circumspheres.

Quality Metrics:
- Aspect Ratio: ratio of longest to shortest edge (ideally close to 1)
- Minimum Angle: should be > 20-30° to avoid numerical issues
- Skewness: deviation from equilateral (0 = perfect, 1 = degenerate)

Mesh Generation Workflow:
-------------------------
1. Define domain boundaries (outer boundary + internal features)
2. Generate boundary points with appropriate spacing
3. Add interior points for resolution control
4. Apply Delaunay triangulation
5. Remove triangles outside the domain
6. Refine/smooth as needed for quality

CFD Mesh Requirements:
---------------------
For finite element/volume methods:
- Boundary layer resolution near walls (for viscous flows)
- Higher resolution in regions of expected gradients
- Smooth size transitions (avoid abrupt changes)
- Element alignment with flow direction (when possible)

CFD Educational Value:
----------------------
This example demonstrates:
- Computational domain definition for CFD
- Point distribution strategies for mesh generation
- Delaunay triangulation for unstructured meshes
- Mesh quality considerations
- VTK mesh visualization techniques
- Boundary handling in mesh generation

References:
-----------
1. Thompson, J.F., et al. "Handbook of Grid Generation"
2. Frey, P.J., George, P.L. "Mesh Generation"
3. Shewchuk, J.R. "Delaunay Refinement Algorithms"
"""

import numpy as np
import vtk


class CFDMeshGenerator:
    """
    Mesh generator for CFD applications using Delaunay triangulation.

    Creates an unstructured triangular mesh for a domain with optional
    internal features (holes, refinement regions).

    Attributes:
        domain_size: Physical size of the square domain (m).
        hole_radius: Radius of central circular hole (m).
        num_boundary_points: Points on outer boundary.
        num_hole_points: Points on hole boundary.
        num_interior_points: Additional interior points.
    """

    def __init__(
        self,
        domain_size: float = 1.0,
        hole_radius: float = 0.2,
        num_boundary_points: int = 40,
        num_hole_points: int = 30,
        num_interior_points: int = 100,
    ):
        """
        Initialize the mesh generator.

        Args:
            domain_size: Size of square domain (m).
            hole_radius: Radius of central hole (m).
            num_boundary_points: Number of points on outer boundary.
            num_hole_points: Number of points on hole boundary.
            num_interior_points: Number of interior seed points.
        """
        self.domain_size = domain_size
        self.hole_radius = hole_radius
        self.num_boundary_points = num_boundary_points
        self.num_hole_points = num_hole_points
        self.num_interior_points = num_interior_points

        self.points = vtk.vtkPoints()
        self.triangle_filter = None
        self.mesh_polydata = None

        self._generate_mesh()

    def _generate_boundary_points(self) -> None:
        """
        Generate points along the outer square boundary.

        Distributes points evenly along the four edges of the domain.
        """
        L = self.domain_size
        n = self.num_boundary_points // 4

        # Bottom edge (y = 0)
        for i in range(n):
            x = i * L / n
            self.points.InsertNextPoint(x, 0.0, 0.0)

        # Right edge (x = L)
        for i in range(n):
            y = i * L / n
            self.points.InsertNextPoint(L, y, 0.0)

        # Top edge (y = L)
        for i in range(n):
            x = L - i * L / n
            self.points.InsertNextPoint(x, L, 0.0)

        # Left edge (x = 0)
        for i in range(n):
            y = L - i * L / n
            self.points.InsertNextPoint(0.0, y, 0.0)

    def _generate_hole_points(self) -> None:
        """
        Generate points along the circular hole boundary.

        Uses uniform angular spacing around the hole perimeter.
        """
        center_x = self.domain_size / 2.0
        center_y = self.domain_size / 2.0
        r = self.hole_radius

        for i in range(self.num_hole_points):
            theta = 2.0 * np.pi * i / self.num_hole_points
            x = center_x + r * np.cos(theta)
            y = center_y + r * np.sin(theta)
            self.points.InsertNextPoint(x, y, 0.0)

    def _generate_interior_points(self) -> None:
        """
        Generate interior points for mesh resolution.

        Uses stratified random sampling to avoid clustering while
        maintaining irregular spacing (beneficial for Delaunay quality).
        Points inside the hole are rejected.
        """
        np.random.seed(42)  # Reproducibility
        center_x = self.domain_size / 2.0
        center_y = self.domain_size / 2.0

        generated = 0
        while generated < self.num_interior_points:
            x = np.random.uniform(0.05, self.domain_size - 0.05)
            y = np.random.uniform(0.05, self.domain_size - 0.05)

            # Reject points inside the hole (with small margin)
            dist_to_center = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
            if dist_to_center > self.hole_radius + 0.02:
                self.points.InsertNextPoint(x, y, 0.0)
                generated += 1

    def _generate_refinement_points(self) -> None:
        """
        Generate additional points near the hole for refinement.

        In CFD, regions near geometric features often require finer meshes
        to capture gradients (boundary layers, stress concentrations).
        """
        center_x = self.domain_size / 2.0
        center_y = self.domain_size / 2.0

        # Two refinement rings around the hole
        for ring_factor in [1.3, 1.6]:
            r = self.hole_radius * ring_factor
            num_points = int(self.num_hole_points * ring_factor)
            for i in range(num_points):
                theta = 2.0 * np.pi * i / num_points + 0.1 * ring_factor
                x = center_x + r * np.cos(theta)
                y = center_y + r * np.sin(theta)
                if 0 < x < self.domain_size and 0 < y < self.domain_size:
                    self.points.InsertNextPoint(x, y, 0.0)

    def _perform_triangulation(self) -> None:
        """
        Perform Delaunay triangulation on the generated points.

        Uses VTK's vtkDelaunay2D for 2D triangulation, which is appropriate
        for planar CFD domains.
        """
        input_polydata = vtk.vtkPolyData()
        input_polydata.SetPoints(self.points)

        delaunay = vtk.vtkDelaunay2D()
        delaunay.SetInputData(input_polydata)
        delaunay.SetTolerance(0.0001)
        delaunay.Update()

        self.triangle_filter = vtk.vtkTriangleFilter()
        self.triangle_filter.SetInputConnection(delaunay.GetOutputPort())
        self.triangle_filter.Update()

        self.mesh_polydata = self.triangle_filter.GetOutput()

    def _remove_hole_triangles(self) -> None:
        """
        Remove triangles that lie inside the circular hole.

        This post-processing step ensures the mesh respects internal boundaries.
        """
        center_x = self.domain_size / 2.0
        center_y = self.domain_size / 2.0

        # Get the triangulated mesh
        mesh = self.triangle_filter.GetOutput()

        # Create a new polydata for filtered triangles
        new_points = vtk.vtkPoints()
        new_cells = vtk.vtkCellArray()

        # Copy all points
        for i in range(mesh.GetNumberOfPoints()):
            new_points.InsertNextPoint(mesh.GetPoint(i))

        # Filter triangles based on centroid location
        for i in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(i)
            if cell.GetCellType() == vtk.VTK_TRIANGLE:
                # Compute triangle centroid
                pts = cell.GetPoints()
                cx = sum(pts.GetPoint(j)[0] for j in range(3)) / 3.0
                cy = sum(pts.GetPoint(j)[1] for j in range(3)) / 3.0

                # Check if centroid is outside the hole
                dist = np.sqrt((cx - center_x) ** 2 + (cy - center_y) ** 2)
                if dist > self.hole_radius:
                    new_cells.InsertNextCell(cell)

        # Create filtered mesh
        filtered_mesh = vtk.vtkPolyData()
        filtered_mesh.SetPoints(new_points)
        filtered_mesh.SetPolys(new_cells)

        # Clean up unused points
        clean = vtk.vtkCleanPolyData()
        clean.SetInputData(filtered_mesh)
        clean.Update()

        self.mesh_polydata = clean.GetOutput()

    def _generate_mesh(self) -> None:
        """
        Generate the complete mesh.

        Workflow:
        1. Generate boundary points
        2. Generate hole points
        3. Generate refinement points
        4. Generate interior points
        5. Perform Delaunay triangulation
        6. Remove triangles inside the hole
        """
        self._generate_boundary_points()
        self._generate_hole_points()
        self._generate_refinement_points()
        self._generate_interior_points()
        self._perform_triangulation()
        self._remove_hole_triangles()

    def compute_mesh_statistics(self) -> dict:
        """
        Compute mesh quality statistics.

        Returns:
            Dictionary with mesh quality metrics.
        """
        mesh = self.mesh_polydata
        num_points = mesh.GetNumberOfPoints()
        num_cells = mesh.GetNumberOfCells()

        # Compute aspect ratios
        aspect_ratios = []
        min_angles = []

        for i in range(num_cells):
            cell = mesh.GetCell(i)
            if cell.GetCellType() == vtk.VTK_TRIANGLE:
                pts = [np.array(cell.GetPoints().GetPoint(j)) for j in range(3)]

                # Edge lengths
                edges = [
                    np.linalg.norm(pts[1] - pts[0]),
                    np.linalg.norm(pts[2] - pts[1]),
                    np.linalg.norm(pts[0] - pts[2]),
                ]
                max_edge = max(edges)
                min_edge = min(edges)

                if min_edge > 0:
                    aspect_ratios.append(max_edge / min_edge)

                    # Compute angles using law of cosines
                    angles = []
                    for j in range(3):
                        a, b, c = edges[j], edges[(j + 1) % 3], edges[(j + 2) % 3]
                        cos_angle = (a**2 + b**2 - c**2) / (2 * a * b + 1e-10)
                        cos_angle = np.clip(cos_angle, -1, 1)
                        angles.append(np.degrees(np.arccos(cos_angle)))
                    min_angles.append(min(angles))

        return {
            "num_points": num_points,
            "num_triangles": num_cells,
            "avg_aspect_ratio": np.mean(aspect_ratios) if aspect_ratios else 0,
            "max_aspect_ratio": max(aspect_ratios) if aspect_ratios else 0,
            "min_angle_deg": min(min_angles) if min_angles else 0,
            "avg_min_angle_deg": np.mean(min_angles) if min_angles else 0,
        }

    def get_actor(self) -> vtk.vtkActor:
        """
        Create a VTK actor for mesh visualization.

        Returns:
            vtkActor with wireframe representation.
        """
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.mesh_polydata)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetRepresentationToWireframe()
        actor.GetProperty().SetLineWidth(1.5)
        actor.GetProperty().SetColor(0.2, 0.6, 1.0)

        return actor

    def get_surface_actor(self) -> vtk.vtkActor:
        """
        Create a VTK actor with surface representation.

        Returns:
            vtkActor with surface shading.
        """
        # Compute normals for shading
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputData(self.mesh_polydata)
        normals.ComputePointNormalsOn()
        normals.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(normals.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.8, 0.8, 0.9)
        actor.GetProperty().SetOpacity(0.8)

        return actor


def visualize_mesh(mesh_generator: CFDMeshGenerator) -> None:
    """
    Visualize the CFD mesh with quality information.

    Displays:
    - Wireframe mesh overlay
    - Surface representation
    - Mesh statistics text

    Args:
        mesh_generator: The CFDMeshGenerator with computed mesh.
    """
    # Get mesh statistics
    stats = mesh_generator.compute_mesh_statistics()
    print("\nMesh Statistics:")
    print(f"  Number of points: {stats['num_points']}")
    print(f"  Number of triangles: {stats['num_triangles']}")
    print(f"  Average aspect ratio: {stats['avg_aspect_ratio']:.2f}")
    print(f"  Maximum aspect ratio: {stats['max_aspect_ratio']:.2f}")
    print(f"  Minimum angle: {stats['min_angle_deg']:.1f}°")
    print(f"  Average minimum angle: {stats['avg_min_angle_deg']:.1f}°")

    # Create actors
    wireframe_actor = mesh_generator.get_actor()
    surface_actor = mesh_generator.get_surface_actor()

    # Create text actor for statistics
    text = f"Triangles: {stats['num_triangles']}\n"
    text += f"Min angle: {stats['min_angle_deg']:.1f}°"

    text_actor = vtk.vtkTextActor()
    text_actor.SetInput(text)
    text_actor.GetTextProperty().SetFontSize(16)
    text_actor.GetTextProperty().SetColor(1, 1, 1)
    text_actor.SetPosition(10, 10)

    # Create renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(surface_actor)
    renderer.AddActor(wireframe_actor)
    renderer.AddActor2D(text_actor)
    renderer.SetBackground(0.1, 0.1, 0.15)

    # Create render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(1024, 768)
    render_window.SetWindowName("CFD Mesh - Delaunay Triangulation")

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

    # Position camera for 2D view
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0.5, 0.5, 3.0)
    camera.SetFocalPoint(0.5, 0.5, 0.0)
    camera.SetViewUp(0, 1, 0)
    renderer.ResetCamera()

    # Start visualization
    interactor.Initialize()
    render_window.Render()
    interactor.Start()


def main():
    """
    Main function to demonstrate CFD mesh generation and visualization.

    Creates a mesh for a plate with a circular hole, a common geometry
    in thermal and structural analysis.
    """
    print("=" * 60)
    print("CFD Mesh Generation - Delaunay Triangulation")
    print("=" * 60)
    print()
    print("Generating mesh for plate with circular hole...")
    print("This geometry is common in heat transfer analysis.")
    print()

    # Create mesh generator
    mesh_gen = CFDMeshGenerator(
        domain_size=1.0,
        hole_radius=0.15,
        num_boundary_points=48,
        num_hole_points=32,
        num_interior_points=150,
    )

    print("Mesh generation complete.")
    print("Use mouse to rotate, zoom, and pan.")
    print()

    # Visualize
    visualize_mesh(mesh_gen)


if __name__ == "__main__":
    main()
