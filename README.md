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
