import vtk

# create a sphere source
sphereSource = vtk.vtkSphereSource()
sphereSource.SetCenter(0, 0, 0)
sphereSource.SetRadius(1)

# create a cube source
cubeSource = vtk.vtkCubeSource()
cubeSource.SetCenter(2, 0, 0)
cubeSource.SetXLength(1)
cubeSource.SetYLength(1)
cubeSource.SetZLength(1)

# create mappers
sphereMapper = vtk.vtkPolyDataMapper()
sphereMapper.SetInputConnection(sphereSource.GetOutputPort())

cubeMapper = vtk.vtkPolyDataMapper()
cubeMapper.SetInputConnection(cubeSource.GetOutputPort())

# create actors
sphereActor = vtk.vtkActor()
sphereActor.SetMapper(sphereMapper)

cubeActor = vtk.vtkActor()
cubeActor.SetMapper(cubeMapper)

# create a renderer and add both actors
renderer = vtk.vtkRenderer()
renderer.AddActor(sphereActor)
renderer.AddActor(cubeActor)

# create a render window and add the renderer
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# create an interactor and set it up to work with the render window
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

# initialize the interactor and start the rendering loop
interactor.Initialize()
interactor.Start()
