import matplotlib.pyplot as plt
import numpy as np
import vtkmodules.all as vtk
from mpl_toolkits.mplot3d import Axes3D

# Create a grid of points for the Rosenbrock Function
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
X, Y = np.meshgrid(x, y)
Z = (1 - X) ** 2 + 100 * (Y - X ** 2) ** 2  # Rosenbrock Function


def create_surface_from_grid(x, y, z):
    # Create a surface from the grid data
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()

    for i in range(len(x)):
        for j in range(len(y)):
            points.InsertNextPoint(x[i][j], y[i][j], z[i][j])
            if i < len(x) - 1 and j < len(y) - 1:
                vertices.InsertNextCell(4)
                vertices.InsertCellPoint(i * len(y) + j)
                vertices.InsertCellPoint((i + 1) * len(y) + j)
                vertices.InsertCellPoint((i + 1) * len(y) + (j + 1))
                vertices.InsertCellPoint(i * len(y) + (j + 1))

    grid = vtk.vtkPolyData()
    grid.SetPoints(points)
    grid.SetPolys(vertices)
    return grid


def main():
    # Create a surface from the grid data
    surface = create_surface_from_grid(X, Y, Z)

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
