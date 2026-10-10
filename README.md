# VTK Examples

A collection of Python examples for VTK, covering geometry, data structures, file formats, rendering, interaction, UI integrations, and small CFD demonstrations. The repository also includes a small browser example using VTK.js.

[![Watch a VTK rendering demo](https://img.youtube.com/vi/Qgoh9NbNqdc/maxresdefault.jpg)](https://youtube.com/shorts/Qgoh9NbNqdc?feature=share)

Many examples follow the VTK rendering flow: **source → filter → mapper → actor → renderer**.

## Getting started

The Python examples require Python 3.10 or newer. From the repository root, create and activate a virtual environment, then install the dependencies:

```bash
git clone https://github.com/djeada/VTK-Examples.git
cd VTK-Examples
# Linux or macOS
python3 -m venv .venv
source .venv/bin/activate
# Windows PowerShell:
# py -m venv .venv
# .venv\Scripts\Activate.ps1

python -m pip install -r requirements-core.txt
```

The core requirements cover VTK, NumPy, Matplotlib, and the CFD examples. Install `requirements.txt` when you also need the Qt, pandas, Spark, or other integrations.

Run a Python example by its filename or relative path. For example:

```bash
python -m src --list
python -m src circle
python -m src 01_basic_shapes/circle
python -m src naca_airfoil_generator -- 0012 --export airfoil.vtp
```

Run commands from the repository root so examples can find the shared `data/` directory. Some examples open an interactive window; close the window to return to the terminal. Examples that read or write external file formats may need sample files in `data/` and can create output files there. File-writing examples may overwrite an existing output file.

## Example catalog

Use `python -m src --list` for the runnable catalog, or filter it with a fragment such as `python -m src --list heat`. Select an example by filename, unique filename fragment, or category/path. The links below point to the source files.

### Basic shapes

- [circle](src/01_basic_shapes/circle.py)
- [cone](src/01_basic_shapes/cone.py)
- [cube](src/01_basic_shapes/cube.py)
- [cylinder](src/01_basic_shapes/cylinder.py)
- [glyph](src/01_basic_shapes/glyph.py)
- [square](src/01_basic_shapes/square.py)
- [triangle](src/01_basic_shapes/triangle.py)

### Advanced shapes

- [custom filter sphere](src/02_advanced_shapes/custom_filter_sphere.py)
- [enclosing box](src/02_advanced_shapes/enclosing_box.py)
- [flow simulation visualization](src/02_advanced_shapes/flow_simulation_visualization.py)
- [isosurface](src/02_advanced_shapes/isosurface.py)
- [multiple dependent objects](src/02_advanced_shapes/multiple_dependent_objects.py)
- [multiple independent objects](src/02_advanced_shapes/multiple_independent_objects.py)
- [streamlines](src/02_advanced_shapes/streamlines.py)
- [triangulation](src/02_advanced_shapes/triangulation.py)
- [visualization techniques comparison](src/02_advanced_shapes/visualization_techniques_comparison.py)
- [volume rendering](src/02_advanced_shapes/volume_rendering.py)

### Structures and datasets

- [Ahmed body visualization](src/03_structures_and_datasets/ahmed_body_visualization.py)
- [cell types demo](src/03_structures_and_datasets/cell_types_demo.py)
- [cells](src/03_structures_and_datasets/cells.py)
- [connectivity demo](src/03_structures_and_datasets/connectivity_demo.py)
- [data types demo](src/03_structures_and_datasets/data_types_demo.py)
- [fields](src/03_structures_and_datasets/fields.py)
- [generate NACA 0012](src/03_structures_and_datasets/generate_naca0012.py)
- [image data](src/03_structures_and_datasets/image_data.py)
- [multiblock dataset](src/03_structures_and_datasets/multiblock_dataset.py)
- [NACA airfoil generator](src/03_structures_and_datasets/naca_airfoil_generator.py)
- [points](src/03_structures_and_datasets/points.py)
- [poly data](src/03_structures_and_datasets/poly_data.py)
- [rectilinear grid](src/03_structures_and_datasets/rectilinear_grid.py)
- [structured grid](src/03_structures_and_datasets/structured_grid.py)
- [structured mesh](src/03_structures_and_datasets/structured_mesh.py)
- [topology vs geometry](src/03_structures_and_datasets/topology_vs_geometry.py)
- [unstructured grid](src/03_structures_and_datasets/unstructured_grid.py)
- [unstructured mesh](src/03_structures_and_datasets/unstructured_mesh.py)

### Input and output

- [Exodus II](src/04_input_output/io_exodus_ii.py)
- [OBJ](src/04_input_output/io_obj.py)
- [OpenFOAM](src/04_input_output/io_openfoam.py)
- [PLY](src/04_input_output/io_ply.py)
- [STL](src/04_input_output/io_stl.py)
- [VTK](src/04_input_output/io_vtk.py)
- [VTK PolyData](src/04_input_output/io_vtk_polydata.py)
- [VTM](src/04_input_output/io_vtm.py)
- [VTP](src/04_input_output/io_vtp.py)
- [VTU](src/04_input_output/io_vtu.py)

### Data conversion

- [OBJ to VTP](src/05_data_conversion/obj_vtp.py)
- [PLY to OBJ](src/05_data_conversion/ply_obj.py)
- [PLY to STL](src/05_data_conversion/ply_stl.py)
- [STL to OBJ](src/05_data_conversion/stl_obj.py)
- [STL to VTK](src/05_data_conversion/stl_vtk.py)
- [STL to VTP](src/05_data_conversion/stl_vtp.py)
- [VTK to OBJ](src/05_data_conversion/vtk_obj.py)
- [VTK to VTM](src/05_data_conversion/vtk_vtm.py)
- [VTK to VTP](src/05_data_conversion/vtk_vtp.py)
- [VTK to VTU](src/05_data_conversion/vtk_vtu.py)

### Visualization pipeline

- [actor and mapper with multiple objects](src/06_visualization_pipeline/actor_mapper_multiple_objects.py)
- [text labels](src/06_visualization_pipeline/adding_text_labels.py)
- [camera movement](src/06_visualization_pipeline/camera_movement.py)
- [cycles decomposition](src/06_visualization_pipeline/claude_cycles_decomposition.py)
- [filters in action](src/06_visualization_pipeline/filters_in_action.py)
- [lighting and shadows](src/06_visualization_pipeline/lighting_and_shadows.py)
- [pipeline animation](src/06_visualization_pipeline/pipeline_animation.py)
- [scalar color mapping](src/06_visualization_pipeline/scalar_color_mapping.py)

### Interactive widgets

- [deformation cylinder](src/07_interactive_widgets/deformation_cylinder.py)
- [interactive flower mesh](src/07_interactive_widgets/interactive_flower_mesh.py)
- [Knuth Hamiltonian cycle](src/07_interactive_widgets/knuth_hamiltonian_cycle.py)
- [orientation marker](src/07_interactive_widgets/orientation_marker.py)
- [planes intersection](src/07_interactive_widgets/planes_intersection.py)
- [simple button](src/07_interactive_widgets/simple_button.py)
- [slider](src/07_interactive_widgets/slider.py)
- [solar eclipse simulation](src/07_interactive_widgets/solar_eclipse_simulation.py)

### UI integrations

- [geospatial data](src/08_integration_with_ui/geospatial_data.py)
- [Matplotlib sphere](src/08_integration_with_ui/matplotlib_sphere.py)
- [Matplotlib surface plot](src/08_integration_with_ui/matplotlib_surface_plot.py)
- [Qt sphere](src/08_integration_with_ui/qt_sphere.py)
- [pandas fallback for a Spark and VTK workflow](src/08_integration_with_ui/spark_vtk_workflow.py)

### CFD

See the [CFD model and visualization guide](src/09_cfd/README.md) for the equations, boundary conditions, checks, and limitations of these teaching examples.

- [convective end heat transfer, 1D](src/09_cfd/convective_end_heat_transfer_1d.py)
- [enhanced heat transfer, 1D](src/09_cfd/enhanced_heat_transfer_solver_1d.py)
- [enhanced heat transfer, 2D](src/09_cfd/enhanced_heat_transfer_solver_2d.py)
- [fixed end heat transfer, 1D](src/09_cfd/fixed_end_heat_transfer_1d.py)
- [fluid flow simulator](src/09_cfd/fluid_flow_simulator.py)
- [heat conduction, 2D](src/09_cfd/heat_conduction_solver_2d.py)
- [heat convection, 1D](src/09_cfd/heat_convection_solver_1d.py)
- [heat convection, 2D](src/09_cfd/heat_convection_solver_2d.py)
- [obstacle flow simulation](src/09_cfd/obstacle_flow_simulation.py)

### Educational science reels

Create one vertical video per selected example. The nine curated CFD stories
use numerical adapters for the example solvers; other examples use their
original VTK scene with a generic story or a supplied JSON story. By default,
reels render at 1080 × 1920 and 30 fps; `fast` uses 720 × 1280 at 24 fps.
Quality also selects different numerical grids for curated models, so fast and
studio results can differ slightly; the manifest records the model parameters.
The curated obstacle reel uses `Re=50`, while the standalone example defaults
to `Re=100`.

Install the reel dependencies and system FFmpeg (including `ffprobe`), then list
or render stories:

```bash
python -m pip install -r requirements-reels.txt
python -m src.reels --list
python -m src.reels obstacle_flow_simulation
python -m src.reels --all cfd
```

Each selected example produces its own MP4, thumbnail, contact sheet, JSON
manifest, and SRT captions under `reels/`. Linux systems without a display use
Xvfb automatically when installed. See the [reel guide](docs/reels.md) for
quality, audio, output, and custom-story options.

## VTK.js example

The browser example lives in [`vtk_js/basic`](vtk_js/basic/README.md). Follow its README to install the Node.js dependencies and start the development server; it prints the local URL to open.

## Notes

The [`notes/`](notes/) directory contains guides on [data types and structures](notes/01_data_types_and_structures.md), [filters and algorithms](notes/02_filters_and_algorithms.md), [input and output](notes/03_input_and_output.md), [visualization techniques](notes/04_visualization_techniques.md), [interactivity](notes/05_interactivity.md), [animations](notes/06_animations.md), [performance](notes/07_performance_optimization_and_parallelism.md), [tool integration](notes/08_integration_with_other_tools.md), and [custom filters](notes/09_custom_filters_and_algorithms.md).

## Contributing

Bug reports, example fixes, and documentation improvements are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the local workflow and project conventions.

## License

This project is licensed under the [MIT License](LICENSE).
