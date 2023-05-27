import vtk

# create a cube source
cubeSource = vtk.vtkCubeSource()
cubeSource.SetCenter(-2, 0, 0)  # Move the cube to the left
cubeSource.SetXLength(1)
cubeSource.SetYLength(1)
cubeSource.SetZLength(1)

# create a mapper for the original cube
cubeMapper = vtk.vtkPolyDataMapper()
cubeMapper.SetInputConnection(cubeSource.GetOutputPort())

# create an actor for the original cube
cubeActor = vtk.vtkActor()
cubeActor.SetMapper(cubeMapper)

# create an elevation filter
elevationFilter = vtk.vtkElevationFilter()
elevationFilter.SetInputConnection(cubeSource.GetOutputPort())
elevationFilter.SetLowPoint(-1, -1, -1)
elevationFilter.SetHighPoint(1, 1, 1)

# create a mapper for the cube after elevation filter
elevationMapper = vtk.vtkPolyDataMapper()
elevationMapper.SetInputConnection(elevationFilter.GetOutputPort())
elevationMapper.SetScalarRange(
    0, 1
)  # To use the scalar data generated from the elevation filter

# create an actor for the cube after elevation filter
elevationActor = vtk.vtkActor()
elevationActor.SetMapper(elevationMapper)
elevationActor.SetPosition(2, 0, 0)  # Move this cube to the center

# create a warp filter
warpFilter = vtk.vtkWarpScalar()
warpFilter.SetInputConnection(elevationFilter.GetOutputPort())
warpFilter.SetScaleFactor(0.5)  # Modify this to increase or decrease the deformation

# create a mapper for the cube after warp filter
warpMapper = vtk.vtkPolyDataMapper()
warpMapper.SetInputConnection(warpFilter.GetOutputPort())
warpMapper.SetScalarRange(
    0, 1
)  # To use the scalar data generated from the elevation filter

# create an actor for the cube after warp filter
warpActor = vtk.vtkActor()
warpActor.SetMapper(warpMapper)
warpActor.SetPosition(4, 0, 0)  # Move this cube to the right

# create a renderer and add the actors
renderer = vtk.vtkRenderer()
renderer.AddActor(cubeActor)
renderer.AddActor(elevationActor)
renderer.AddActor(warpActor)

# create a render window and add the renderer
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# create an interactor and set it up to work with the render window
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

# initialize the interactor and start the rendering loop
interactor.Initialize()
interactor.Start()
