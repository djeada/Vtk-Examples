import matplotlib.pyplot as plt
import vtkmodules.all as vtk
from vtkmodules.util import numpy_support


def create_sphere(radius=1.0, theta_resolution=40, phi_resolution=40):
    # Create a sphere
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(theta_resolution)
    sphere.SetPhiResolution(phi_resolution)
    sphere.Update()
    return sphere


def convert_vtk_to_numpy(vtk_data):
    # Convert vtk to numpy
    numpy_array = numpy_support.vtk_to_numpy(vtk_data)
    return numpy_array


def main():
    # Create a sphere
    sphere = create_sphere()

    # Convert vtk to numpy
    vtk_array = sphere.GetOutput().GetPoints().GetData()
    numpy_array = convert_vtk_to_numpy(vtk_array)

    # Split the numpy array into x, y, z components for 3D plotting
    x, y, z = numpy_array[:, 0], numpy_array[:, 1], numpy_array[:, 2]

    # Plot the sphere using matplotlib
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, color="b", alpha=0.6, edgecolors="w", s=20)
    plt.show()


if __name__ == "__main__":
    main()
