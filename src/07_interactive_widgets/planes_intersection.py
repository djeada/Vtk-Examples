"""
This module provides a set of tools to visualize the intersection of planes in 3D space using VTK (Visualization Toolkit).
It facilitates the creation and rendering of planes and their intersection lines. The module is capable of handling multiple planes
and their intersections, providing an intuitive visual representation of complex geometrical relationships.

Workflow:
1. Initialization: An instance of PlaneIntersectionVisualizer is created with a list of plane parameters.
2. Plane Creation: For each set of plane parameters, a VTK plane actor is created.
3. Intersection Computation: The module computes the intersection lines of all possible plane pairs.
4. Clipping Intersections: The intersection lines are clipped to lie within the bounds of the planes.
5. Visualization: The planes and their intersection lines are rendered in a 3D space using VTK.

Mathematics:

- Plane Representation:
  Each plane is defined by three points: an origin (O) and two other points (P1, P2) on the plane.
  A plane can be represented in the form ax + by + cz = d, where [a, b, c] is the plane's normal vector.

- Normal Calculation:
  The normal (N) of the plane is calculated using the cross product of vectors formed by the plane points.
  N = (P1 - O) × (P2 - O), where × denotes the cross product.

- Intersection Line:
  The intersection line of two planes is found by solving the system of equations formed by the normals (N1, N2) of the planes and a direction vector (D).
  The direction vector D is the cross product of the normals: D = N1 × N2.
  The line of intersection can be represented by a point (P) and a direction (D): Line = P + tD, where t is a parameter.

  To find P, we solve the system of equations:
    N1 · P = N1 · O1
    N2 · P = N2 · O2
    D · P = 0
  where O1 and O2 are the origins of the two planes.

- Line-Plane Intersection:
  The intersection of a line with a plane is computed using the parametric equation of the line and the plane equation.
  Let L(t) = L0 + tV be the parametric equation of the line, where L0 is a point on the line, V is the direction vector of the line, and t is a parameter.
  The plane equation is N · X = d, where N is the normal of the plane, and X is any point on the plane.

  To find the intersection, we substitute the line equation into the plane equation:
    N · (L0 + tV) = d
  Solving for t gives the intersection point on the line.

- Clipping:
  Intersection lines are clipped to the bounds of the planes by checking if intersection points lie on the plane edges.
  This is done by ensuring the intersection point falls within the segment defined by the edge endpoints.

Usage Example:
To use the module, define the parameters for the planes you wish to visualize and create an instance of PlaneIntersectionVisualizer.
Call the visualize method to render the planes and their intersections in 3D space.
"""
from itertools import combinations

import numpy as np
import vtk


class PlaneIntersectionVisualizer:
    def __init__(self, plane_parameters):
        self.plane_parameters = plane_parameters
        self.plane_actors = []
        self.intersection_actors = []

    @staticmethod
    def create_plane_actor(origin, point1, point2, color):
        plane = vtk.vtkPlaneSource()
        plane.SetOrigin(origin)
        plane.SetPoint1(point1)
        plane.SetPoint2(point2)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(plane.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        return actor

    @staticmethod
    def create_line_actor(p1, p2, color, line_width=4):
        line = vtk.vtkLineSource()
        line.SetPoint1(p1)
        line.SetPoint2(p2)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(line.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetLineWidth(line_width)

        return actor

    def compute_and_clip_intersections(self):
        for (p1_params, p2_params) in combinations(self.plane_parameters, 2):
            line_point, line_dir = self.compute_plane_intersection(
                *p1_params[:3], *p2_params[:3]
            )

            for plane_params in [p1_params, p2_params]:
                p1, p2 = self.clip_line_to_plane_bounds(
                    line_point, line_dir, *plane_params[:3]
                )
                if p1 is not None and p2 is not None:
                    intersection_actor = self.create_line_actor(
                        p1, p2, [1, 1, 0]
                    )  # Yellow color for intersections
                    self.intersection_actors.append(intersection_actor)
                    break

    @staticmethod
    def compute_plane_intersection(
        p1_origin, p1_point1, p1_point2, p2_origin, p2_point1, p2_point2
    ):
        # Calculate normals
        normal1 = np.cross(
            np.array(p1_point1) - np.array(p1_origin),
            np.array(p1_point2) - np.array(p1_origin),
        )
        normal2 = np.cross(
            np.array(p2_point1) - np.array(p2_origin),
            np.array(p2_point2) - np.array(p2_origin),
        )

        A = np.array([normal1, normal2, np.cross(normal1, normal2)])
        b = np.array([np.dot(normal1, p1_origin), np.dot(normal2, p2_origin), 0])
        point_on_line = np.linalg.solve(A, b)
        direction = np.cross(normal1, normal2)

        return point_on_line, direction

    @staticmethod
    def clip_line_to_plane_bounds(
        line_point, line_dir, plane_origin, plane_point1, plane_point2
    ):
        plane_origin, plane_point1, plane_point2 = map(
            np.array, [plane_origin, plane_point1, plane_point2]
        )
        line_point, line_dir = np.array(line_point), np.array(line_dir)

        plane_normal = np.cross(
            plane_point1 - plane_origin, plane_point2 - plane_origin
        )

        edges = [
            (plane_origin, plane_point1),
            (plane_point1, plane_point2),
            (plane_point2, plane_origin),
        ]
        intersection_points = []
        for edge_start, edge_end in edges:
            edge_dir = edge_end - edge_start
            edge_normal = np.cross(plane_normal, edge_dir)

            intersection = PlaneIntersectionVisualizer.line_plane_intersection(
                line_point, line_dir, edge_start, edge_normal
            )
            if (
                intersection is not None
                and PlaneIntersectionVisualizer.is_point_on_edge(
                    intersection, edge_start, edge_end
                )
            ):
                intersection_points.append(intersection)

        if len(intersection_points) == 2:
            return intersection_points[0], intersection_points[1]
        else:
            return None, None

    @staticmethod
    def line_plane_intersection(line_point, line_dir, plane_point, plane_normal):
        denominator = np.dot(plane_normal, line_dir)
        if np.isclose(denominator, 0):
            return None

        d = np.dot(plane_normal, (plane_point - line_point)) / denominator
        intersection = line_point + d * line_dir

        if np.any(np.isnan(intersection)) or np.any(np.isinf(intersection)):
            return None

        return intersection

    @staticmethod
    def is_point_on_edge(point, edge_start, edge_end, tolerance=1e-6):
        edge_vector = edge_end - edge_start
        point_vector = point - edge_start
        cross_product = np.cross(edge_vector, point_vector)
        dot_product = np.dot(point_vector, edge_vector)
        edge_length_squared = np.dot(edge_vector, edge_vector)
        return (
            np.linalg.norm(cross_product) <= tolerance
            and 0 <= dot_product <= edge_length_squared
        )

    def visualize(self):
        self.plane_actors = [
            self.create_plane_actor(*params) for params in self.plane_parameters
        ]
        self.compute_and_clip_intersections()

        renderer = vtk.vtkRenderer()
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window_interactor = vtk.vtkRenderWindowInteractor()
        render_window_interactor.SetRenderWindow(render_window)

        for actor in self.plane_actors + self.intersection_actors:
            renderer.AddActor(actor)

        renderer.SetBackground(0.1, 0.2, 0.4)
        render_window.Render()
        render_window_interactor.Start()


# Example usage
plane_parameters = [
    ([0, -1, 0], [1, -1, 0], [0, 1, 0], [1, 0, 0]),
    ([0, 0, -1], [1, 0, -1], [0, 0, 1], [0, 1, 0]),
    ([0, -1, 0], [1, -1, 0], [0, 1, 1], [1, 0, 1]),
]

visualizer = PlaneIntersectionVisualizer(plane_parameters)
visualizer.visualize()
