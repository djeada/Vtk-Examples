"""Browse cross-sections through a synthetic volume (not medical MRI)."""
import math

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk


def make_volume(size=64):
    axis = np.linspace(-1, 1, size)
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    shell = (x / 0.83) ** 2 + (y / 0.9) ** 2 + (z / 0.75) ** 2
    center = (x / 0.45) ** 2 + (y / 0.6) ** 2 + (z / 0.45) ** 2
    pocket = ((x - 0.28) / 0.2) ** 2 + ((y + 0.18) / 0.25) ** 2 + (z / 0.3) ** 2
    density = np.where(shell < 1, 100, 0).astype(np.float32)
    density[center < 1] = 205
    density[pocket < 1] = 40
    image = vtk.vtkImageData()
    image.SetDimensions(size, size, size)
    image.SetSpacing(*(2.0 / (size - 1),) * 3)
    image.SetOrigin(-1, -1, -1)
    pixels = numpy_to_vtk(density.ravel(order="F"), deep=True)
    pixels.SetName("synthetic_intensity")
    image.GetPointData().SetScalars(pixels)
    return image


def build_scene(size=64):
    image = make_volume(size)
    mapper = vtk.vtkImageSliceMapper()
    mapper.SetInputData(image)
    mapper.SetOrientationToZ()
    mapper.SetSliceNumber(size // 2)
    slice_actor = vtk.vtkImageSlice()
    slice_actor.SetMapper(mapper)
    slice_actor.GetProperty().SetColorWindow(220)
    slice_actor.GetProperty().SetColorLevel(110)
    outline = vtk.vtkOutlineFilter()
    outline.SetInputData(image)
    outline_mapper = vtk.vtkPolyDataMapper()
    outline_mapper.SetInputConnection(outline.GetOutputPort())
    box = vtk.vtkActor()
    box.SetMapper(outline_mapper)
    box.GetProperty().SetColor(0.20, 0.78, 0.95)
    renderer = vtk.vtkRenderer()
    renderer.AddViewProp(slice_actor)
    renderer.AddActor(box)
    renderer.SetBackground(0.015, 0.022, 0.045)
    renderer.GetActiveCamera().SetPosition(3, -4, 3)
    renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
    renderer.GetActiveCamera().SetViewUp(0, 0, 1)
    frame = [0]

    def step():
        frame[0] += 1
        position = (math.sin(frame[0] * 0.027) + 1) / 2
        mapper.SetSliceNumber(round(position * (size - 1)))

    return renderer, step


def main():
    renderer, step = build_scene()
    window = vtk.vtkRenderWindow()
    window.SetWindowName("Synthetic Volume Scanner")
    window.SetSize(900, 700)
    window.AddRenderer(renderer)
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)

    def tick(obj, event):
        step()
        window.Render()

    interactor.AddObserver("TimerEvent", tick)
    window.Render()
    interactor.Initialize()
    interactor.CreateRepeatingTimer(50)
    interactor.Start()


if __name__ == "__main__":
    main()
