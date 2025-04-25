import numpy as np
import vtk

from src.common.simple_pipeline import VisualisationPipeline


def create_rectilinear_grid(x_coords, y_coords, z_coords):
    """
    Create and return a rectilinear grid.
    """
    rgrid = vtk.vtkRectilinearGrid()
    rgrid.SetDimensions(len(x_coords), len(y_coords), len(z_coords))

    # Set X coordinates
    x_array = vtk.vtkFloatArray()
    for x in x_coords:
        x_array.InsertNextValue(x)
    rgrid.SetXCoordinates(x_array)

    # Set Y coordinates
    y_array = vtk.vtkFloatArray()
    for y in y_coords:
        y_array.InsertNextValue(y)
    rgrid.SetYCoordinates(y_array)

    # Set Z coordinates
    z_array = vtk.vtkFloatArray()
    for z in z_coords:
        z_array.InsertNextValue(z)
    rgrid.SetZCoordinates(z_array)

    return rgrid


def attach_fields_to_grid(rgrid):
    """
    Attach scalar and vector fields to the rectilinear grid.
    """
    # Attach a scalar field
    scalars = vtk.vtkFloatArray()
    scalars.SetName("Scalars")
    num_points = rgrid.GetNumberOfPoints()
    for i in range(num_points):
        scalars.InsertNextValue(np.random.rand())
    rgrid.GetPointData().SetScalars(scalars)

    # Attach a vector field
    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)
    vectors.SetName("Vectors")
    for i in range(num_points):
        vectors.InsertNextTuple3(np.random.rand(), np.random.rand(), np.random.rand())
    rgrid.GetPointData().SetVectors(vectors)


def visualize_rectilinear_grid(rgrid):
    """
    Visualize the rectilinear grid with scalar and vector fields.
    """
    geometry_filter = vtk.vtkRectilinearGridGeometryFilter()
    geometry_filter.SetInputData(rgrid)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(geometry_filter.GetOutputPort())
    mapper.SetScalarModeToUsePointData()
    mapper.SetColorModeToMapScalars()
    mapper.SelectColorArray("Scalars")

    pipeline = VisualisationPipeline(mappers=[mapper], point_size=30)
    pipeline.run()


def main():
    x_coords = np.linspace(0, 10, 11)
    y_coords = np.linspace(0, 5, 6)
    z_coords = np.linspace(0, 3, 4)
    rgrid = create_rectilinear_grid(x_coords, y_coords, z_coords)
    attach_fields_to_grid(rgrid)
    visualize_rectilinear_grid(rgrid)


if __name__ == "__main__":
    main()
