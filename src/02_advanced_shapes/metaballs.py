"""Merging implicit blobs rendered as an animated isosurface."""
import math

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy


def make_field(size=48):
    image = vtk.vtkImageData()
    image.SetDimensions(size, size, size)
    image.SetOrigin(-2, -2, -2)
    image.SetSpacing(*(4 / (size - 1),) * 3)
    scalars = numpy_to_vtk(np.zeros(size**3, dtype=np.float32), deep=True)
    scalars.SetName("density")
    image.GetPointData().SetScalars(scalars)
    axis = np.linspace(-2, 2, size)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")

    def update(time):
        centers = [
            (-0.8 + 0.48 * math.cos(time), 0.30 * math.sin(time), 0),
            (0.8 - 0.48 * math.cos(time), 0.15, 0.25 * math.sin(time)),
            (0, -0.7 + 0.4 * math.sin(time * 0.7), -0.3),
        ]
        field = np.zeros_like(x, dtype=np.float32)
        for cx, cy, cz in centers:
            field += np.exp(-2.4 * ((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2))
        vtk_to_numpy(scalars)[:] = field.ravel(order="F")
        scalars.Modified()
        image.Modified()

    update(0)
    return image, update


def build_scene(size=48):
    image, update = make_field(size)
    contour = vtk.vtkFlyingEdges3D()
    contour.SetInputData(image)
    contour.SetValue(0, 0.52)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(contour.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.22, 0.84, 0.93)
    actor.GetProperty().SetSpecular(0.65)
    actor.GetProperty().SetSpecularPower(35)
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0.015, 0.025, 0.055)
    frame = [0]

    def step():
        frame[0] += 1
        update(frame[0] * 0.045)
        actor.RotateY(0.35)

    return renderer, step


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Merging Metaballs")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Zoom(1.3)

    def tick(obj, event):
        step()
        window.Render()

    interactor.AddObserver("TimerEvent", tick)
    window.Render()
    interactor.Initialize()
    interactor.CreateRepeatingTimer(40)
    interactor.Start()


if __name__ == "__main__":
    main()
