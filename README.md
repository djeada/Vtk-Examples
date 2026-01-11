<div align="center">

# 🎨 VTK Examples - Master 3D Visualization with Python

**A comprehensive collection of VTK (Visualization Toolkit) examples showcasing 3D graphics, scientific visualization, and computational fluid dynamics**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![VTK 9.4](https://img.shields.io/badge/VTK-9.4-green.svg)](https://vtk.org/)
[![Contributions Welcome](https://img.shields.io/badge/Contributions-Welcome-brightgreen.svg)](CONTRIBUTING.md)

[📖 Documentation](#-documentation) • [🚀 Quick Start](#-quick-start) • [💡 Examples](#-examples-gallery) • [🤝 Contributing](#-how-to-contribute)

---

![VTK Visualization Demo](https://github.com/djeada/VTK-Examples/assets/37275728/d843d55f-ab7a-4830-a762-34cc71041fdc)

</div>

## 🎬 See VTK in Action

Experience the power of VTK through these interactive demonstrations. Click on any thumbnail to watch the video on YouTube:

<div align="center">

| Visualization Demo | Animation Showcase |
|:------------------:|:------------------:|
| [![Watch Demo 1](https://img.youtube.com/vi/Qgoh9NbNqdc/maxresdefault.jpg)](https://youtube.com/shorts/Qgoh9NbNqdc?feature=share) | [![Watch Demo 2](https://img.youtube.com/vi/0jAN9Q-GGCk/maxresdefault.jpg)](https://youtube.com/shorts/0jAN9Q-GGCk?feature=share) |
| *3D Rendering & Interaction* | *Dynamic Visualization* |

</div>

---

## 📋 Table of Contents

- [✨ Why This Repository?](#-why-this-repository)
- [🔍 What is VTK?](#-what-is-vtk)
- [🚀 Quick Start](#-quick-start)
- [📁 Project Structure](#-project-structure)
- [💡 Examples Gallery](#-examples-gallery)
  - [🔷 Basic Shapes](#-basic-shapes)
  - [🔶 Advanced Shapes](#-advanced-shapes)
  - [📊 Structures and Datasets](#-structures-and-datasets)
  - [📥 Input and Output](#-input-and-output)
  - [🔄 Data Conversion](#-data-conversion)
  - [🎯 Visualization Pipeline](#-visualization-pipeline)
  - [🎛️ Interactive Widgets](#-interactive-widgets)
  - [🔗 Integrations](#-integrations-with-external-tools)
  - [🌊 Computational Fluid Dynamics](#-computational-fluid-dynamics-cfd)
- [🌐 VTK.js Web Examples](#-vtkjs-web-examples)
- [📖 Documentation](#-documentation)
- [🤝 How to Contribute](#-how-to-contribute)
- [📚 References](#-references)
- [📄 License](#-license)

---

## ✨ Why This Repository?

<div align="center">

| 🎓 **Learn by Example** | 🛠️ **Production Ready** | 📈 **Comprehensive Coverage** |
|:-----------------------:|:------------------------:|:-----------------------------:|
| 50+ working examples with clear explanations | Well-structured, reusable code patterns | From basics to advanced CFD simulations |

| 🔬 **Scientific Focus** | 🌐 **Web & Desktop** | 📚 **Rich Documentation** |
|:-----------------------:|:--------------------:|:-------------------------:|
| Perfect for research and engineering | Python + VTK.js examples | In-depth notes and tutorials |

</div>

Whether you're a **student** learning visualization concepts, a **researcher** analyzing scientific data, or an **engineer** building professional applications, this repository provides the building blocks you need.

### 🎯 What You'll Learn

- **Create stunning 3D visualizations** from simple shapes to complex scientific data
- **Master the VTK pipeline** architecture for efficient rendering
- **Handle various file formats** (STL, OBJ, VTK, VTU, VTM, and more)
- **Build interactive applications** with widgets and UI integration
- **Implement CFD simulations** for heat transfer and fluid flow
- **Deploy web visualizations** using VTK.js

---

## 🔍 What is VTK?

The **Visualization Toolkit (VTK)** is the gold standard for 3D computer graphics, image processing, and scientific visualization.

<table>
<tr>
<td width="50%">

### 🏆 Features

- **Extensive Algorithm Library** — Hundreds of visualization algorithms for any use case
- **Cross-Platform** — Windows, Linux, macOS support
- **Multi-Language** — C++, Python, Java, Tcl bindings
- **Open Source** — BSD license, free for commercial use
- **Active Community** — 30+ years of development by Kitware

</td>
<td width="50%">

### 🔧 Core Capabilities

- 3D surface and volume rendering
- Scalar, vector, and tensor visualization
- Image processing and analysis
- Mesh generation and manipulation
- Scientific data formats I/O
- Interactive 3D widgets

</td>
</tr>
</table>

### 📊 VTK Pipeline Architecture

Understanding VTK's pipeline is essential for effective visualization:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         VTK VISUALIZATION PIPELINE                      │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐           │
│  │  SOURCE  │───▶│  FILTER  │───▶│  MAPPER  │───▶│  ACTOR   │           │
│  │          │    │          │    │          │    │          │           │
│  │ Generate │    │ Process  │    │ Convert  │    │ Render   │           │
│  │   Data   │    │   Data   │    │ to Geo   │    │  Object  │           │
│  └──────────┘    └──────────┘    └──────────┘    └────┬─────┘           │
│                                                       │                 │
│                                                       ▼                 │
│                                               ┌──────────────┐          │
│                                               │   RENDERER   │          │
│                                               │              │          │
│                                               │   Display    │          │
│                                               └──────────────┘          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 🚀 Quick Start

Get up and running in minutes with these simple steps:

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/djeada/VTK-Examples.git
cd VTK-Examples

# 2. Create a virtual environment (recommended)
python -m venv venv

# 3. Activate the virtual environment
# On Linux/macOS:
source venv/bin/activate
# On Windows:
venv\Scripts\activate

# 4. Install dependencies
pip install -r requirements.txt
```

### Run Your First Example

```bash
# Navigate to the examples directory and run a script
cd src/01_basic_shapes
python circle.py
```

🎉 **Congratulations!** You should see a window displaying circles rendered with VTK.

### 💻 IDE Setup (Optional)

For an enhanced development experience with **VS Code** or **PyCharm**:

1. **Open** the project folder in your IDE
2. **Configure** the Python interpreter to use your virtual environment
3. **Run** any script using the built-in run button or debugger
4. **Debug** with breakpoints for step-through analysis

---

## 📁 Project Structure

```
VTK-Examples/
├── 📂 src/                          # Source code examples
│   ├── 📂 01_basic_shapes/          # Fundamental geometric primitives
│   ├── 📂 02_advanced_shapes/       # Complex geometries & techniques
│   ├── 📂 03_structures_and_datasets/ # VTK data structures
│   ├── 📂 04_input_output/          # File format handling
│   ├── 📂 05_data_conversion/       # Format conversion utilities
│   ├── 📂 06_visualization_pipeline/ # Rendering pipeline examples
│   ├── 📂 07_interactive_widgets/   # UI widgets & interaction
│   ├── 📂 08_integration_with_ui/   # Qt, Matplotlib integration
│   ├── 📂 09_cfd/                   # Computational Fluid Dynamics
│   └── 📂 common/                   # Shared utilities & helpers
│
├── 📂 notes/                        # In-depth documentation & tutorials
├── 📂 vtk_js/                       # VTK.js web visualization examples
├── 📂 data/                         # Sample data files
│   ├── 📂 stls/                     # STL mesh files
│   ├── 📂 objs/                     # OBJ model files
│   ├── 📂 vtks/                     # VTK data files
│   └── 📂 ...                       # Other format samples
│
├── 📄 requirements.txt              # Python dependencies
├── 📄 LICENSE                       # MIT License
└── 📄 README.md                     # This file
```

---

## 💡 Examples Gallery

### 🔷 Basic Shapes

Start your VTK journey with fundamental geometric primitives. These examples demonstrate the core concepts of creating and rendering 3D objects.

| # | Example | Description | Link |
|:-:|:--------|:------------|:----:|
| 1 | **Circle** | Create perfectly round 2D circles with customizable radius and resolution | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/circle.py) |
| 2 | **Cone** | Generate 3D cones with adjustable height, radius, and resolution | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cone.py) |
| 3 | **Cube** | Build solid cubes and boxes in 3D space | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cube.py) |
| 4 | **Cylinder** | Create cylindrical shapes with height and radius parameters | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/cylinder.py) |
| 5 | **Glyph** | Introduction to glyph-based visualization for representing data points | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/glyph.py) |
| 6 | **Square** | Render 2D squares and rectangles | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/square.py) |
| 7 | **Triangle** | Create triangular primitives, the building blocks of 3D meshes | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/01_basic_shapes/triangle.py) |

### 🔶 Advanced Shapes

Take your visualization skills to the next level with complex geometries, volume rendering, and multi-object scenes.

| # | Example | Description | Link |
|:-:|:--------|:------------|:----:|
| 1 | **Enclosing Box** | Create bounding boxes that enclose other objects—useful for spatial analysis | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/enclosing_box.py) |
| 2 | **Isosurface** | Extract surfaces of constant value from volumetric data | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/isosurface.py) |
| 3 | **Multiple Dependent Objects** | Build scenes with hierarchically linked objects | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/multiple_dependent_objects.py) |
| 4 | **Multiple Independent Objects** | Create complex scenes with multiple standalone objects | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/multiple_independent_objects.py) |
| 5 | **Streamlines** | Visualize flow fields and vector data with streamline rendering | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/streamlines.py) |
| 6 | **Triangulation** | Apply triangulation techniques for mesh generation | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/triangulation.py) |
| 7 | **Volume Rendering** | Render volumetric data with customizable transfer functions | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/volume_rendering.py) |
| 8 | **Visualization Comparison** | Compare different visualization techniques side-by-side | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/visualization_techniques_comparison.py) |
| 9 | **Flow Simulation** | Visualize computational fluid dynamics data | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/02_advanced_shapes/flow_simulation_visualization.py) |

### 📊 Structures and Datasets

Master VTK's data structures—the foundation of all visualization operations.

| # | Example | Description | Link |
|:-:|:--------|:------------|:----:|
| 1 | **Points** | Work with vtkPoints—the fundamental spatial data structure | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/points.py) |
| 2 | **Cells** | Understand cells and their role in defining mesh topology | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/cells.py) |
| 3 | **Fields** | Handle scalar, vector, and tensor field data | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/fields.py) |
| 4 | **Multiblock Dataset** | Organize complex data with hierarchical multiblock structures | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/multiblock_dataset.py) |
| 5 | **PolyData** | Work with polygonal surfaces—the most common VTK data type | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/poly_data.py) |
| 6 | **Structured Grid** | Handle regularly spaced 3D grid data | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/structured_grid.py) |
| 7 | **Unstructured Grid** | Work with arbitrary cell types and connectivity | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/unstructured_grid.py) |
| 8 | **Structured Mesh** | Create and manipulate structured mesh geometries | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/structured_mesh.py) |
| 9 | **Unstructured Mesh** | Build and process unstructured mesh geometries | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/unstructured_mesh.py) |
| 10 | **Cell Types Demo** | Interactive demo of all VTK cell types with combo box selection | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/cell_types_demo.py) |
| 11 | **Topology vs Geometry** | Interactive demo showing the difference between mesh connectivity and point positions | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/03_structures_and_datasets/topology_vs_geometry.py) |

### 📥 Input and Output

Learn to read and write various 3D file formats commonly used in engineering and scientific applications.

| # | Format | Description | Link |
|:-:|:-------|:------------|:----:|
| 1 | **OBJ** | Wavefront OBJ—widely used 3D model format | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_obj.py) |
| 2 | **STL** | Stereolithography format—standard for 3D printing | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_stl.py) |
| 3 | **VTK** | Native VTK format for polygonal data | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtk.py) |
| 4 | **VTM** | VTK multiblock format for composite datasets | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtm.py) |
| 5 | **VTU** | VTK unstructured grid format | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_vtu.py) |
| 6 | **Exodus II** | SANDIA National Labs format for FEM/CFD results | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_exodus_ii.py) |
| 7 | **OpenFOAM** | Read OpenFOAM CFD simulation results | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/04_input_output/io_openfoam.py) |

### 🔄 Data Conversion

Convert between popular 3D file formats with ease.

| # | Conversion | Description | Link |
|:-:|:-----------|:------------|:----:|
| 1 | **Converter Interface** | Base interface for building custom converters | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/converter_interface.py) |
| 2 | **STL ↔ VTK** | Convert between STL and VTK formats | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/stl_vtk.py) |
| 3 | **VTK ↔ OBJ** | Convert between VTK and OBJ formats | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_obj.py) |
| 4 | **STL ↔ OBJ** | Convert between STL and OBJ formats | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/stl_obj.py) |
| 5 | **VTK ↔ VTM** | Convert between VTK and VTM formats | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_vtm.py) |
| 6 | **VTK ↔ VTU** | Convert between VTK and VTU formats | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/05_data_conversion/vtk_vtu.py) |

### 🎯 Visualization Pipeline

Understand and master the VTK rendering pipeline for professional-quality visualizations.

| # | Example | Description | Link |
|:-:|:--------|:------------|:----:|
| 1 | **Actor-Mapper Setup** | Handle multiple objects with proper actor-mapper relationships | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/actor_mapper_multiple_objects.py) |
| 2 | **Text Labels** | Add informative text annotations to your visualizations | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/adding_text_labels.py) |
| 3 | **Scalar Color Mapping** | Map data values to colors using customizable color maps | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/scalar_color_mapping.py) |
| 4 | **Camera Movement** | Control camera position, orientation, and movement paths | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/camera_movement.py) |
| 5 | **Filters in Action** | Apply VTK filters for data processing and transformation | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/filters_in_action.py) |
| 6 | **Lighting & Shadows** | Create realistic lighting effects and shadow rendering | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/lighting_and_shadows.py) |
| 7 | **Pipeline Animation** | Animate your visualization pipeline over time | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/06_visualization_pipeline/pipeline_animation.py) |

### 🎛️ Interactive Widgets

Build interactive applications with VTK's powerful widget system.

| # | Widget | Description | Link |
|:-:|:-------|:------------|:----:|
| 1 | **Orientation Marker** | Add 3D axes indicator for spatial orientation | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/orientation_marker.py) |
| 2 | **Slider** | Create interactive sliders for parameter control | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/slider.py) |
| 3 | **Button** | Implement clickable buttons for user interaction | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/simple_button.py) |
| 4 | **Planes Intersection** | Interactively slice volumes with intersecting planes | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/07_interactive_widgets/planes_intersection.py) |

### 🔗 Integrations with External Tools

Combine VTK with popular frameworks for enhanced visualization capabilities.

| # | Integration | Description | Link |
|:-:|:------------|:------------|:----:|
| 1 | **Qt Integration** | Embed VTK visualizations in Qt desktop applications | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/08_integration_with_ui/qt_sphere.py) |
| 2 | **Matplotlib Sphere** | Combine VTK 3D rendering with Matplotlib | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/08_integration_with_ui/matplotlib_sphere.py) |
| 3 | **Matplotlib Surface** | Create surface plots using VTK and Matplotlib together | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/08_integration_with_ui/matplotlib_surface_plot.py) |

### 🌊 Computational Fluid Dynamics (CFD)

Simulate and visualize heat transfer and fluid flow phenomena with these physics-based examples.

| # | Simulation | Description | Link |
|:-:|:-----------|:------------|:----:|
| 1 | **1D Heat Convection** | Solve 1D heat convection problems with finite differences | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/heat_convection_solver_1d.py) |
| 2 | **1D Fixed-End Heat Transfer** | Heat transfer with fixed temperature boundary conditions | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/fixed_end_heat_transfer_1d.py) |
| 3 | **1D Convective-End Heat Transfer** | Heat transfer with convective boundary conditions | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/convective_end_heat_transfer_1d.py) |
| 4 | **1D Enhanced Heat Transfer** | Advanced 1D heat transfer with enhanced accuracy | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/enhanced_heat_transfer_solver_1d.py) |
| 5 | **2D Heat Conduction** | Solve 2D steady-state heat conduction problems | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/heat_conduction_solver_2d.py) |
| 6 | **2D Enhanced Heat Transfer** | Advanced 2D heat transfer solver | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/enhanced_heat_transfer_solver_2d.py) |
| 7 | **Fluid Flow Simulation** | Simulate basic fluid flow patterns | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/fluid_flow_simulator.py) |
| 8 | **Obstacle Flow** | Simulate fluid flow around obstacles | [📄 Code](https://github.com/djeada/VTK-Examples/blob/main/src/09_cfd/obstacle_flow_simulation.py) |

---

## 🌐 VTK.js Web Examples

Take your visualizations to the web with **VTK.js**—a JavaScript implementation of VTK for browser-based 3D rendering.

### Quick Start with VTK.js

```bash
# Navigate to the VTK.js examples
cd vtk_js/basic

# Install dependencies
npm install

# Build the project
npm run build

# Start the development server
npm start
```

Then open your browser to `http://localhost:8080` to see VTK running in the web!

📂 **Location:** [`vtk_js/basic/`](https://github.com/djeada/VTK-Examples/tree/main/vtk_js/basic)

---

## 📖 Documentation

Deepen your understanding with our comprehensive documentation covering VTK concepts, techniques, and best practices.

| Topic | Description | Link |
|:------|:------------|:----:|
| **Data Types & Structures** | vtkDataObject, points, cells, grids, and datasets | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/01_data_types_and_structures.md) |
| **Filters & Algorithms** | Data processing, transformation, and analysis filters | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/02_filters_and_algorithms.md) |
| **Input & Output** | File format handling and data serialization | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/03_input_and_output.md) |
| **Visualization Techniques** | Volume rendering, scalar mapping, vector visualization | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/04_visualization_techniques.md) |
| **Interactivity** | User interaction, picking, and selection | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/05_interactivity.md) |
| **Animations** | Creating smooth animations and time-series visualizations | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/06_animations.md) |
| **Performance Optimization** | Parallelism, LOD, and rendering optimization | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/07_performance_optimization_and_parallelism.md) |
| **Tool Integration** | Integrating VTK with Qt, Matplotlib, and other tools | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/08_integration_with_other_tools.md) |
| **Custom Filters** | Building custom filters and algorithms | [📚 Read](https://github.com/djeada/VTK-Examples/blob/main/notes/09_custom_filters_and_algorithms.md) |

---

## 🤝 How to Contribute

We welcome contributions from the community! Whether you're fixing a bug, adding new examples, or improving documentation, your help is appreciated.

### Ways to Contribute

- 🐛 **Report bugs** — Found an issue? [Open a bug report](https://github.com/djeada/VTK-Examples/issues/new)
- 💡 **Suggest features** — Have an idea? [Share it with us](https://github.com/djeada/VTK-Examples/issues/new)
- 📝 **Improve documentation** — Help make the docs even better
- 🔧 **Submit code** — Add new examples or enhance existing ones

### Contribution Workflow

```bash
# 1. Fork the repository on GitHub

# 2. Clone your fork
git clone https://github.com/YOUR_USERNAME/VTK-Examples.git
cd VTK-Examples

# 3. Create a feature branch
git checkout -b feature/AmazingFeature

# 4. Make your changes and commit
git add .
git commit -m "Add AmazingFeature: brief description"

# 5. Push to your fork
git push origin feature/AmazingFeature

# 6. Open a Pull Request on GitHub
```

### Code Guidelines

- Follow existing code style and patterns
- Add clear comments and docstrings
- Include example usage in new scripts
- Test your code before submitting

---

## 📚 References

Expand your VTK knowledge with these excellent resources:

### Official Resources
- 🌐 [VTK Official Website](https://vtk.org/) — Documentation, downloads, and news
- 📖 [VTK User's Guide](https://vtk.org/vtk-users-guide/) — Comprehensive official guide
- 📚 [VTK Examples](https://examples.vtk.org/) — Official VTK examples collection

### Community Resources
- 📝 [PyScience VTK Tutorials](https://pyscience.wordpress.com/tag/vtk/) — Python-focused VTK tutorials
- 💻 [Program Creek VTK Examples](https://www.programcreek.com/python/index/480/vtk) — Code snippets and examples
- 🎓 [Visualization Course 2014 (Uppsala)](https://www.cb.uu.se/~aht/Vis2014/) — Academic visualization course
- 🎓 [Data Visualization Course (SJTU)](https://www.cs.sjtu.edu.cn/~shengbin/course/datavis/) — University course materials

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

```
MIT License - Copyright (c) 2021 Adam Djellouli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software.
```

---

## ⭐ Star History

If you find this repository helpful, please consider giving it a star! Your support helps others discover these resources.

[![Star History Chart](https://api.star-history.com/svg?repos=djeada/Vtk-Examples&type=Date)](https://star-history.com/#djeada/Vtk-Examples&Date)

---

<div align="center">

**Made with ❤️ by [Adam Djellouli](https://github.com/djeada) and contributors**

[⬆ Back to Top](#-vtk-examples---master-3d-visualization-with-python)

</div>
