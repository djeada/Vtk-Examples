# Ten bite-sized VTK projects for YouTube

These are small, visually striking **Python + VTK** projects to build on camera. Aim for a clear *minimum viable effect* before adding lighting, text, or cinematic polish. The result can be a 30–90-second Short or the central reveal of a 5–10-minute tutorial.

The times below are rough targets for a basic demo **once you know the rendering boilerplate**; animated, highly polished versions will take longer. All ten can use generated data, so no external 3D assets are required.

## The project shortlist

| # | Project / video hook | First visible result | Useful VTK building blocks | MVP target |
| --- | --- | --- | --- | --- |
| 1 | **Lava planet** — "I made a planet in Python" | A crater-like sphere with fiery height colors | `vtkSphereSource`, point scalars, `vtkWarpScalar`, `vtkPolyDataNormals` | 10–20 min |
| 2 | **Merging metaballs** — "These blobs become one object" | Smooth blobs fuse and separate | `vtkImageData`, `vtkFlyingEdges3D` or `vtkContourFilter` | 20–30 min |
| 3 | **Black-hole gravity well** — "I bent a 3D grid" | A wireframe surface dips toward a dark center | `vtkPoints`, `vtkPolyData`, wireframe `vtkActor` | 10–20 min |
| 4 | **Tornado ribbons** — "A tornado made from curves" | Helices spiral around a narrowing funnel | `vtkPolyLine`, `vtkTubeFilter`, `vtkPolyDataMapper` | 10–20 min |
| 5 | **Live 3D slicer** — "Watch me cut through a sphere" | A moving plane reveals the intersection ring | `vtkPlane`, `vtkCutter`, optional `vtkImplicitPlaneWidget2` | 5–10 min |
| 6 | **Exploding mesh** — "I separated every triangle" | A solid sphere flies apart into facets | `vtkTriangleFilter`, `vtkPolyData`, `vtkPoints` | 15–25 min |
| 7 | **Procedural terrain** — "No heightmap, just noise" | A flat grid becomes mountains | `vtkPlaneSource`, height scalars, `vtkWarpScalar` | 10–20 min |
| 8 | **Faux MRI scanner** — "Scan a volume slice by slice" | A slice moves through a glowing synthetic volume | `vtkImageData`, `vtkImageReslice`, `vtkImageSlice` | 15–25 min |
| 9 | **Boolean geometry** — "Subtract a cylinder from a sphere" | A clean geometric hole appears | `vtkTriangleFilter`, `vtkBooleanOperationPolyDataFilter` | 10–20 min |
| 10 | **Recursive fractal sculpture** — "Geometry inside geometry" | Repeated miniature shapes form a 3D pattern | `vtkTransformPolyDataFilter`, `vtkAppendPolyData` | 15–25 min |

## How to build each effect

1. **Lava planet:** Start with a high-resolution sphere. Assign deterministic noise or a mathematical height function to its points, displace the surface, then map the heights to dark rock/orange/red. Slowly orbit the camera. This is *procedural terrain*, not a lava-flow simulation.
2. **Merging metaballs:** Fill a small 3D image grid with the sum of several Gaussian fields centered on moving points. Extract an isosurface at one threshold and update the centers; the surface merges naturally. Begin with a low grid resolution to keep frame times manageable.
3. **Black-hole gravity well:** Construct a rectangular mesh whose height follows a steep radial depression, then render its edges. Move a small sphere around the center for scale. This is a *stylized gravity-well diagram*, not a general-relativistic or light-bending simulation.
4. **Tornado ribbons:** Generate several helices whose radius and height change along the path. Turn each polyline into a colored tube and rotate its phase per frame. This visualizes a designed spiral, not simulated fluid dynamics.
5. **Live 3D slicer:** Feed a sphere to a cutter driven by a `vtkPlane`. Draw the resulting contour in a contrasting color over a translucent copy of the sphere. Slide the plane origin to make the circular cross-section grow and shrink.
6. **Exploding mesh:** Triangulate a sphere. Copy each triangle with its **own three points** (so neighboring faces are no longer shared), then move its points outward along the triangle-centroid direction. Interpolate the offset from zero to one for the reveal.
7. **Procedural terrain:** Generate a plane grid, compute one height per vertex from layered sines or noise, and warp in the plane-normal direction. Recompute normals for smooth lighting; add a height-based color ramp.
8. **Faux MRI scanner:** Make a synthetic 3D voxel field (for example, nested ellipsoids with different intensities). Display one resliced plane at a time and animate the slice offset. Label it **synthetic volume**: this is an imaging visualization, not real patient data.
9. **Boolean geometry:** Triangulate two **closed, intersecting** surfaces and subtract one from the other. Boolean polygon filters can be sensitive to non-manifold geometry or near-coplanar faces; if the operation fails, slightly change the overlap or use implicit-function sampling plus contour extraction instead.
10. **Recursive fractal sculpture:** Repeatedly create scaled and translated instances of one simple source at child locations. Append the transformed meshes, stop after a few levels, and animate the camera. Cap recursion depth because instance count grows exponentially.

## First build: a tiny slice-through-sphere demo

Run this from the repository root after installing the [core requirements](../README.md#getting-started). The circle at the cut plane is the exact surface intersection computed by `vtkCutter`. You can animate `plane.SetOrigin(x, 0, 0)` and re-render to turn it into a Short.

```python
import vtk

sphere = vtk.vtkSphereSource()
sphere.SetThetaResolution(64)
sphere.SetPhiResolution(64)

plane = vtk.vtkPlane()
plane.SetOrigin(0.25, 0, 0)
plane.SetNormal(1, 0, 0)

cutter = vtk.vtkCutter()
cutter.SetCutFunction(plane)
cutter.SetInputConnection(sphere.GetOutputPort())

shell_mapper = vtk.vtkPolyDataMapper()
shell_mapper.SetInputConnection(sphere.GetOutputPort())
shell = vtk.vtkActor()
shell.SetMapper(shell_mapper)
shell.GetProperty().SetOpacity(0.25)

cut_mapper = vtk.vtkPolyDataMapper()
cut_mapper.SetInputConnection(cutter.GetOutputPort())
cut_actor = vtk.vtkActor()
cut_actor.SetMapper(cut_mapper)
cut_actor.GetProperty().SetColor(1, 0.25, 0.1)
cut_actor.GetProperty().SetLineWidth(5)

renderer = vtk.vtkRenderer()
renderer.AddActor(shell)
renderer.AddActor(cut_actor)
renderer.SetBackground(0.04, 0.06, 0.11)

window = vtk.vtkRenderWindow()
window.AddRenderer(renderer)
window.SetSize(900, 700)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)
renderer.ResetCamera()
window.Render()
interactor.Start()
```

## A simple YouTube shot list

1. **0–3 s:** Show the final effect first ("I made a tornado with VTK").
2. **3–10 s:** Show the primitive or field that everything starts from.
3. **10–40 s:** Show two or three transformations: source → filter → mapper → actor.
4. **40–60 s:** Reveal the animation, rotate the camera, and display one meaningful takeaway.
5. **For longer tutorials:** Explain the one important idea (e.g., an isosurface or a cut plane), not every line of boilerplate.

For recording and video export, start with [animations](06_animations.md) or the repository's [educational reel workflow](../docs/reels.md). The reel tool renders **existing runnable scenes**; ideas on this page become reel inputs only after they have been implemented as scene examples.

## Relevant examples already in this repository

- [Isosurface](../src/02_advanced_shapes/isosurface.py) — a natural starting point for metaballs.
- [Triangulation](../src/02_advanced_shapes/triangulation.py) — inspect face topology before an explosion.
- [Image data](../src/03_structures_and_datasets/image_data.py) — build a synthetic voxel volume.
- [Pipeline animation](../src/06_visualization_pipeline/pipeline_animation.py) — use repeatable scene updates.
- [Planes intersection](../src/07_interactive_widgets/planes_intersection.py) — explore cutting planes and interaction.

**Fastest first video:** the live slicer (#5). **Strongest visual reveal:** metaballs (#2) or the exploding mesh (#6).
