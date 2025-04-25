"""
This module demonstrates the visualization of a dipole field in a three-dimensional space using the Visualization Toolkit (VTK). A dipole field is a vector field that is typically created by two equal but opposite charges separated by a distance.

Workflow:
1. Generate a structured grid of points in a 3D space.
2. Calculate the dipole field vectors at each point in the grid.
3. Combine these points and vectors into a vtkStructuredGrid object.
4. Create seed points for streamline initiation.
5. Generate streamlines representing the flow lines of the dipole field using vtkStreamTracer.
6. Visualize these streamlines using vtkPolyDataMapper, vtkActor, and vtkRenderer.

Mathematics:
- The dipole field at a point in space is calculated using the formula for the electric field due to point charges. The field due to each charge is \vec{E} = \frac{q \vec{r}}{4 \pi \epsilon_0 r^3}, where \vec{r} is the position vector relative to the charge, q is the charge magnitude, and \epsilon_0 is the vacuum permittivity.
- The net field at any point is the vector sum of the fields due to both charges.
- Streamlines are generated using vtkStreamTracer, which integrates the vector field to trace the paths of particles in the flow.
"""

import numpy as np
import vtk


class DipoleFieldVisualization:
    def __init__(self, num_seed_points=30):
        self.num_seed_points = num_seed_points
        self.grid = self.create_structured_grid()
        self.seeds = self.create_seeds()
        self.streamer = self.create_streamlines()

    def create_structured_grid(self):
        points, vectors = self.generate_points_and_vectors()
        grid = vtk.vtkStructuredGrid()
        grid.SetDimensions(21, 21, 21)
        grid.SetPoints(points)
        grid.GetPointData().SetVectors(vectors)
        return grid

    def generate_points_and_vectors(self):
        points = vtk.vtkPoints()
        vectors = vtk.vtkDoubleArray()
        vectors.SetNumberOfComponents(3)

        for i in np.linspace(-10, 10, 21):
            for j in np.linspace(-10, 10, 21):
                for k in np.linspace(-10, 10, 21):
                    points.InsertNextPoint(i, j, k)
                    x, y, z = i, j, k
                    r1 = np.sqrt((x - 1) ** 2 + y**2 + z**2) + 1e-2
                    r2 = np.sqrt((x + 1) ** 2 + y**2 + z**2) + 1e-2
                    vectors.InsertNextTuple3(
                        (x - 1) / r1**3 - (x + 1) / r2**3,
                        y / r1**3 - y / r2**3,
                        z / r1**3 - z / r2**3,
                    )

        return points, vectors

    def create_seeds(self):
        seeds = vtk.vtkPointSource()
        seeds.SetCenter(0, 0, 0)
        seeds.SetRadius(10)
        seeds.SetNumberOfPoints(self.num_seed_points)
        return seeds

    def create_streamlines(self):
        streamer = vtk.vtkStreamTracer()
        streamer.SetInputData(self.grid)
        streamer.SetSourceConnection(self.seeds.GetOutputPort())
        streamer.SetMaximumPropagation(500)
        streamer.SetIntegrationStepUnit(vtk.vtkStreamTracer.LENGTH_UNIT)
        streamer.SetInitialIntegrationStep(0.1)
        streamer.SetIntegrationDirectionToBoth()
        streamer.SetComputeVorticity(False)
        return streamer

    def visualize(self):
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self.streamer.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        renderer = vtk.vtkRenderer()
        renderer.AddActor(actor)
        renderer.SetBackground(0, 0, 0)

        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)

        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(render_window)
        interactor.Initialize()
        interactor.Start()


if __name__ == "__main__":
    visualization = DipoleFieldVisualization()
    visualization.visualize()
