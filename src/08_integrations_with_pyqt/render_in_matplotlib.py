import matplotlib.pyplot as plt
import vtk
from vtkmodules.util import numpy_support

# Create a sphere
sphere = vtk.vtkSphereSource()
sphere.SetRadius(1.0)
sphere.SetThetaResolution(40)
sphere.SetPhiResolution(40)
sphere.Update()

# Convert vtk to numpy
vtk_array = sphere.GetOutput().GetPoints().GetData()
numpy_array = numpy_support.vtk_to_numpy(vtk_array)

# Split the numpy array into x, y, z components for 3D plotting
x = numpy_array[:, 0]
y = numpy_array[:, 1]
z = numpy_array[:, 2]

# Plot the sphere using matplotlib
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.scatter(x, y, z, color="b", alpha=0.6, edgecolors="w", s=20)
plt.show()
