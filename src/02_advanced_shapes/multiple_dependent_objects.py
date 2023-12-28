import vtk


def create_cube_mapper(center: tuple, length: float) -> vtk.vtkPolyDataMapper:
    cube_source = vtk.vtkCubeSource()
    cube_source.SetCenter(center)
    cube_source.SetXLength(length)
    cube_source.SetYLength(length)
    cube_source.SetZLength(length)

    cube_mapper = vtk.vtkPolyDataMapper()
    cube_mapper.SetInputConnection(cube_source.GetOutputPort())

    return cube_mapper


def create_cone_mapper(
    center: tuple, height: float, radius: float, resolution: int = 30
) -> vtk.vtkPolyDataMapper:
    cone_source = vtk.vtkConeSource()
    cone_source.SetCenter(center)
    cone_source.SetHeight(height)
    cone_source.SetRadius(radius)
    cone_source.SetResolution(resolution)

    cone_mapper = vtk.vtkPolyDataMapper()
    cone_mapper.SetInputConnection(cone_source.GetOutputPort())

    return cone_mapper


if __name__ == "__main__":
    cube_mapper = create_cube_mapper(center=(0, 0, 0), length=1.0)
    cone_mapper = create_cone_mapper(center=(0, 1.5, 0), height=1.5, radius=0.5)

    cube_actor = vtk.vtkActor()
    cube_actor.SetMapper(cube_mapper)

    cone_actor = vtk.vtkActor()
    cone_actor.SetMapper(cone_mapper)

    assembly = vtk.vtkAssembly()
    assembly.AddPart(cube_actor)
    assembly.AddPart(cone_actor)

    renderer = vtk.vtkRenderer()
    renderer.AddActor(assembly)
    renderer.SetBackground(0, 0, 0)

    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(800, 600)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    interactor.Initialize()

    interactor.Start()
