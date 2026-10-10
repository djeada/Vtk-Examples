import matplotlib.pyplot as plt
import numpy as np
import vtkmodules.all as vtk
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # Registers the "3d" projection.

# Create a grid of points for the Rosenbrock Function
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)
Z = (1 - X) ** 2 + 100 * (Y - X**2) ** 2  # Rosenbrock Function


def create_surface_from_grid(x, y, z):
    """Create quad cells from matching 2D coordinate arrays."""
    x, y, z = (np.asarray(array) for array in (x, y, z))
    if x.ndim != 2 or x.shape != y.shape or x.shape != z.shape:
        raise ValueError("Coordinate arrays must have the same 2D shape.")
    rows, columns = x.shape
    points = vtk.vtkPoints()
    polygons = vtk.vtkCellArray()

    for i in range(rows):
        for j in range(columns):
            points.InsertNextPoint(x[i][j], y[i][j], z[i][j])
            if i < rows - 1 and j < columns - 1:
                polygons.InsertNextCell(4)
                polygons.InsertCellPoint(i * columns + j)
                polygons.InsertCellPoint((i + 1) * columns + j)
                polygons.InsertCellPoint((i + 1) * columns + (j + 1))
                polygons.InsertCellPoint(i * columns + (j + 1))

    grid = vtk.vtkPolyData()
    grid.SetPoints(points)
    grid.SetPolys(polygons)
    return grid


def main():
    # Create a surface from the grid data
    surface = create_surface_from_grid(X, Y, Z)
    print(
        f"VTK surface: {surface.GetNumberOfPoints()} points, "
        f"{surface.GetNumberOfCells()} quads"
    )

    # Create a Matplotlib 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot the surface
    ax.plot_surface(
        X, Y, Z, cmap="viridis", rstride=1, cstride=1, linewidth=0, antialiased=False
    )

    # Customize the plot (add labels, title, etc.)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("Rosenbrock Function")

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
