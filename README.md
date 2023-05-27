# VTK

Python examples for the Visualization Toolkit (VTK).

## Overview

* VTK is a powerful high-level 3D graphics and visualization library.
* Widely used for scientific data visualization, 3D graphics, and image processing.
* Provides hundreds of algorithms for various purposes.
* First released in 1993, actively developed and maintained by Kitware Inc.
* Open Source, written in C++, with wrappers for Python, Java, and Tcl.

## Getting Started

### Installation

To install VTK on various platforms (Windows, Linux, macOS), visit the official download page for platform-specific installation guides: https://vtk.org/download/

### Basic Workflow

1. Data input and representation: Select the appropriate data structure for your data.
2. Data manipulation and processing: Apply filters and algorithms to process the data.
3. Visualization and rendering: Use mappers, actors, and renderers to display the data.
4. Interaction and exporting: Add user interactivity and export the visualizations as needed.

## Repository Structure

This repository contains examples and tutorials to help users get started with VTK in Python. Each subdirectory contains a specific topic, along with corresponding code examples.

## Running Examples

To run the examples in this repository, follow the steps below:

### Setting Up the Environment

1. Make sure to have a virtual environment set up for your project. You can use `virtualenv` or `conda` to create a virtual environment.
2. Activate the virtual environment and install the necessary packages. For example:

```
pip install -r requirements.txt
```

### Running Individual Scripts

To run an individual example script, navigate to the appropriate folder within the `src/` directory and execute the Python script using the command line. For example:

```
cd src/01_basic_shapes
python circle.py
```

This will run the `circle.py` script, which creates a simple circle using VTK.

### Using an IDE

If you prefer to use an Integrated Development Environment (IDE) like Visual Studio Code or PyCharm, follow these steps:

1. Open the IDE and load the project folder.
2. Set up the Python interpreter for the project by selecting the virtual environment you created earlier.
3. Open the example script you would like to run.
4. Run the script using the IDE's built-in tools (e.g., the "Run" button or right-clicking on the script and selecting "Run").

With the IDE, you can also set breakpoints, step through the code, and debug the example scripts as needed.

## Basic Shapes

| Number | Description | Link |
| --- | --- | --- |
| 1 | Demonstrates the steps to create a perfectly round circle | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/circle.py) |
| 2 | Describes the process of generating a 3D cone | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cone.py) |
| 3 | Explore how to make a simple cube in 3D space | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cube.py) |
| 4 | Learn to draw a cylinder using vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cylinder.py) |
| 5 | Get started with glyph production in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/glyph.py) |
| 6 | Tutorial on creating a 2D square | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/square.py) |
| 7 | Guides you through the creation of a simple triangle | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/triangle.py) |

## Advanced Shapes

| Number | Description | Link |
| --- | --- | --- |
| 1 | A detailed example of creating a 3D box that encloses other objects | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/enclosing_box.py) |
| 2 | Shows how to construct an isosurface in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/isosurface.py) |
| 3 | Delves into creating scenes with multiple dependent objects | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/multiple_dependent_objects.py) |
| 4 | Explores creating complex scenes with multiple independent objects | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/multiple_independent_objects.py) |
| 5 | Introduction to the creation of streamlines in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/streamlines.py) |
| 6 | Teaches you how to use triangulation techniques in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/triangulation.py) |
| 7 | Explores the world of volume rendering in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/volume_rendering.py) |
| 8 | A comparison of different visualization techniques in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/visualization_techniques_comparison.py) |

## Structures and Datasets

| Number | Description | Link |
| --- | --- | --- |
| 1 | Shows you how to manipulate points in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/points.py) |
| 2 | Dives into the workings of cells within vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/cells.py) |
| 3 | An overview of handling fields in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/fields.py) |
| 4 | Detailed example of working with multiblock datasets | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/multiblock_dataset.py) |
| 5 | Covers the basics of poly data structures in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/poly_data.py) |
| 6 | Teaches you how to work with structured grids | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/structured_grid.py) |
| 7 | Guides you through the intricacies of unstructured grids | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/unstructured_grid.py) |

## Input and Output

| Number | Description | Link |
| --- | --- | --- |
| 1 | Walkthrough of handling Exodus II files | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_exodus_ii.py) |
| 2 | Introduction to working with OBJ files in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_obj.py) |
| 3 | Tutorial on interacting with OpenFOAM files | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_openfoam.py) |
| 4 | Learn to read and write STL files with vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_stl.py) |
| 5 | Overview of handling VTK file format | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtk.py) |
| 6 | Guides you through the usage of VTM file format | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtm.py) |
| 7 | Explores handling VTU file formats with vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtu.py) |

## Data Conversion

| Number | Description | Link |
| --- | --- | --- |
| 1 | Interface for data conversion utilities | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/converter_interface.py) |
| 2 | Teaches conversion between STL and OBJ formats | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/stl_obj.py) |
| 3 | Learn to convert between STL and VTK formats | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/stl_vtk.py) |
| 4 | Demonstrates conversion between VTK and OBJ formats | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_obj.py) |
| 5 | Shows you how to convert between VTK and VTM formats | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_vtm.py) |
| 6 | Guides you through the conversion between VTK and VTU formats | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_vtu.py) |

## Visualization Pipeline

| Number | Description | Link |
| --- | --- | --- |
| 1 | Shows handling multiple objects in the actor-mapper setup | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/actor_mapper_multiple_objects.py) |
| 2 | Teaches you how to add text labels in your visualization | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/adding_text_labels.py) |
| 3 | Walks you through creating camera movements | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/camera_movement.py) |
| 4 | Shows various filters in action in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/filters_in_action.py) |
| 5 | Guides you through creating lighting and shadows in your visualization | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/lighting_and_shadows.py) |
| 6 | Shows you how to animate your visualization pipeline | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/pipeline_animation.py) |
| 7 | Demonstrates scalar color mapping in vtk | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/scalar_color_mapping.py) |

## Interactive Widgets

| Number | Description | Link |
| --- | --- | --- |
| 1 | Guides you to use the orientation marker widget | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/orientation_marker.py) |
| 2 | Learn to create slicing planes in your visualization | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/slicing_planes.py) |
| 3 | Tutorial on creating and using sliders in your visualization | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/slider.py) |

## Integrations with external tools

| Number | Description | Link |
| --- | --- | --- |
| 1 | Learn to create a 3D window with PyQt | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/08_integrations_with_pyqt/3d_window.py) |
| 2 | Shows you how to render vtk in matplotlib | [Python](https://github.com/djeada/VTK-Examples/blob/main/src/08_integrations_with_pyqt/render_in_matplotlib.py) |

## References

* [PyScience: VTK](https://pyscience.wordpress.com/tag/vtk/)
* [Program Creek: VTK](https://www.programcreek.com/python/index/480/vtk)
* [Visualization Course 2014](https://www.cb.uu.se/~aht/Vis2014/)
* [Data Visualization Course, SJTU](https://www.cs.sjtu.edu.cn/~shengbin/course/datavis/)

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
