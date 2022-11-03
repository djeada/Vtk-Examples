# VTK
VTK examples in Python.

## Overview 

* Powerfull higl-level visualization library (with 3D graphics).
* Used for 3D graphics and image processing.
* Hundreds of implemented algorithms.
* It was first released in 1993 and is still being actively developed today.
* Open Source (mainteined by Kitware Inc).
* Written in C++ (there are wrappers for Python, Java and Tcl).

## Basic idea

Everything that you can do with VTK comes down to the following:

* Data (many different native types).
* Data manipulation methods (filters that can transform the data in many ways).
* IO (loading and exporting data to uncounted number of file formats).
* Visualization (methods of displaying the data on the screen). 

## VTK Library

The VTK library is enormous. How can we find anything useful there? Follow the guidelines below:

* VTK/Common -> core classes.
* VTK/IO -> reading and writing data files.
* VTK/Rendering -> a variety of rendering techinques.
* VTK/Graphics -> 3D data processing filter.
* VTK/Filtering -> pipeline for data processing.
* VTK/Imaging -> filters for image processing.
* VTK/Parallel -> support for parallel processing (MPI).

## Data types

### vtkDataObject 

* vtkImageData
* vtkRectlinearGrid (n blocks of the same size)
* vtkStructuredGrid
* vtkPolydata <- points
* vtkUnstructuredGrid <- combines any number of geometries

### vtkpolydata

### structued vs unstrucutred grid

* Unstrucutred grid can contain any kind of geometric primitive.
* You can mix different structures and put them in that container (let's say a line, a cube and two cones).
* All the parts will share the data.

### multiblock dataset


### PyVTK

* non graphic library.
* Used to convert data from calculations run in Python to a common format
* Besides different formats (Points, Cell, Field) you can also give names to your data (which will be shown in paraview.
* Fields (additional info for the cells).

## Filters and algorithms

### vtkAlgorithm 

* Source
  - Procedural source
  - Reader source  
* Filter

Source -> Data object

Data object -> Filter -> Data object


## IO

### VTK file format

ASCII or binary;

The VTK file format consists of five parts; the first three are mandatory, and the last two are optional.

    # vtk DataFile Version m.n
        The header line. "m.n" is the version number. 
    a title
        The title may be up to 256 characters, terminated by a new line. 
    ASCII or BINARY
        indicates the format used for subsequent data. 
    DATASET type
        values for type are "STRUCTURED_POINTS" or "STRUCTURED_GRID" or "UNSTRUCTURED_GRID" or "POLYDATA" or "RECTILINEAR_GRID"; Depending on the type chosen, there will be further lines of keywords and values to define the data. 
    POINT_DATA n
        The number of data items n of each type must match the number of points in the dataset. 

When describing the cells in terms of point indices, the points must be indexed starting at 0.

All point data must use 3 coordinates. Even if your model has 2D geometry, each point must have (X,Y,Z) coordinates. Simply set Z = 0 if you want to work in 2D. 


* File version and identifier
  - \# vtkDataFile Version x,x
* Header: comment, 256 char max
* Type: ASCII or BINARY
* Structure: specified the structure of the domain:
  - STRUCTURED_POINTS
  - STRUCTURED_GRID
  - UNSTRUCTURED_GRID
  - RECTLINEAR_GRID
  - POLYDATA
  - FIELD
  + Information related to the type of the strucutre (dimensions, spacing, coordinates).
  
![cell_types](https://user-images.githubusercontent.com/37275728/194726441-8361df2e-10c3-4792-a16f-4679d5178d05.png)

* Attributes: the values stored at the grid points
  - supported types: scalars, vectors, 3x3 tensors, normals, texture coordinates, and fields.
  - attributesType = {SCALARS, VECTORS, TENSORS, NORMALS, TEXTURE_COORDINATES}
  - dataName = Any name you like
  - dataType = {char, short, int, long, float, double}
  - attributeValues = the values stored at the grid points
  
## Visualization pipeline

Source -> Filter -> Mapper -> Actor -> Renderer

    main() {  
      create geometry;
      create a mapper;
      give the geometry to the mapper;
      create an actor;
      give the mapper to the actor;
      create a renderer;
      give the actor to the renderer;
      create a window;
      give the renderer to the window;
    }

### Mappers

* Mappers convert data into graphical primitives or write to a file (writer).
* Mapper require one or more input data object.
* Converts pure data to renderable geometry.

### Actors

* Actors represent graphical data or objects
* Object properties (color, shading type, etc.), geometry, transformations
* They work together with vtkLight and vtkCamera to make a scene
* The scene is rendered to an image by a renderer.

### Render Window

* The class, vtkRenderWindow ties the entire rendering process together.
* It manages all the platform dependnt window managment issues and hide the details from the user.
* It also stores graphics specific information such as window size, position, title, frame buffer depth, etc.


## Refrences

* https://pyscience.wordpress.com/tag/vtk/
* https://www.programcreek.com/python/index/480/vtk
* https://www.cb.uu.se/~aht/Vis2014/
* https://www.cs.sjtu.edu.cn/~shengbin/course/datavis/

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
