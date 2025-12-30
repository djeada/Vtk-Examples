## Data Types and Structures

VTK is built to carry *real-world* 2D/3D data all the way from “numbers in memory” to “something you can see and reason about.” That means it needs data types that don’t just store values, but also store **where** those values live in space and **how** they connect. If you pick the right structure early, everything downstream, filters, rendering, memory usage, speed, gets smoother. Pick the wrong one, and you’ll spend time converting, losing information, or wondering why performance tanks.

VTK uses 3D geometries, including points, lines, polygons, and volumes. It handles images and volumetric data for 2D and 3D visualization. It works with scalar, vector, and tensor fields for complex data representation. Supports structured and unstructured grid types for various spatial data layouts. Includes support for time-series data and hierarchical datasets.

The mental model to keep in your pocket is simple: **data in VTK is more than values**. It’s values + geometry + topology + metadata. That combination is what lets VTK do “smart visualization” instead of just drawing triangles.

### vtkDataObject

Before you meet the “fancy” datasets, it helps to understand the trunk of the tree. `vtkDataObject` is VTK’s universal container: it’s the “everything starts here” base class. You care because once you know what *every* data object guarantees, you can navigate the VTK ecosystem without feeling like you’re memorizing random class names.

**Do:** treat `vtkDataObject` as the shared contract. When you’re writing general pipelines or utilities, thinking at this level keeps your code flexible.
**Don’t:** assume every `vtkDataObject` is directly renderable or has points/cells, some are composite containers, some are specialized.

I. Base class for all data objects in VTK

II. Stores data and metadata including

* *Geometric Information* encompasses shape, position, and size of the data elements. Essential for defining the spatial representation of the data.
* *Topological Information* concerns the connectivity of data elements. It reveals the relationships between elements, like how points are linked to form lines or polygons in a mesh.

III. Common subclasses include

1. **vtkImageData**: Handles regular, rectilinear grid data.
2. **vtkRectilinearGrid**: Manages grid data with varying spacing.
3. **vtkStructuredGrid**: Deals with structured data points in 3D space.
4. **vtkPolyData**: Specializes in representing polygonal data.
5. **vtkUnstructuredGrid**: Ideal for representing complex, irregular grid data.

When people say “VTK dataset,” they’re often talking about nodes lower in this hierarchy, especially `vtkDataSet` and anything that has **points + cells**. The diagram below is useful because it tells you what you can safely expect. If you’re holding a `vtkPointSet`, for example, you know it has points. If you’re holding a `vtkCompositeDataSet`, you know you’re dealing with *multiple* datasets grouped together.

```
|
|--- vtkDataObject
     |
     |--- vtkDataSet
          |
          |--- vtkPointSet
          |    |
          |    |--- vtkPolyData
          |    |
          |    |--- vtkStructuredPoints (vtkImageData)
          |    |
          |    |--- vtkStructuredGrid
          |    |
          |    |--- vtkUnstructuredGrid
          |
          |--- vtkRectilinearGrid
          |
          |--- vtkCompositeDataSet
               |
               |--- vtkHierarchicalBoxDataSet
               |
               |--- vtkMultiBlockDataSet
               |
               |--- vtkHierarchicalDataSet
               |
               |--- vtkOverlappingAMR
               |
               |--- vtkNonOverlappingAMR
```

### vtkImageData

`vtkImageData` is the “pixel/voxel brain” of VTK. It’s what you reach for when your data lives on a **uniform grid**, like images (2D) or volumes (3D). You care because this structure is fast and compact: connectivity is *implicit*, and a lot of VTK’s imaging and volume rendering pipeline is built around it.

**Do:** use `vtkImageData` for CT/MRI, 3D textures, scalar fields sampled uniformly, and anything that feels like “voxels.”
**Don’t:** force irregular spacing into it. The moment spacing isn’t uniform, you’ll either distort the meaning of your data or end up wishing you chose a rectilinear or unstructured grid.

* Represents a regular grid with fixed topology and uniform spacing between points.
* Ideal for volumetric data like images or 3D scalar fields.
* CT and MRI scans, where vtkImageData stores pixel values at each grid point for 3D visualization.

![image\_data](https://github.com/djeada/Vtk-Examples/assets/37275728/1dd504e1-9532-4863-bd88-29f03ca6e449)

### vtkRectilinearGrid

Think of `vtkRectilinearGrid` as “structured, but with stretchy spacing.” Topology is still a neat grid, but the points along each axis can be unevenly spaced. This matters when your data’s resolution changes by design, common in scientific models where you want more detail in certain regions without paying the cost everywhere else.

**Do:** use this when you still want the simplicity of a structured grid, but your sampling is non-uniform.
**Don’t:** confuse this with curvilinear grids, rectilinear still means axes-aligned lines; spacing changes, not direction.

* A regular grid with fixed topology but non-uniform spacing between points.
* Suitable for data with varying resolution, like in climate or terrain elevation data.
* Climate models with varying altitude resolution, or terrain data with resolution changing with slope.

![rectilinear\_grid](https://github.com/djeada/Vtk-Examples/assets/37275728/1fd2f697-73aa-40aa-a319-302f37751b63)

### vtkPolyData

If `vtkImageData` is “voxels,” `vtkPolyData` is “surfaces.” This is the workhorse for meshes you’d actually model or export: triangle surfaces, polygon soups, lines, contours, and anything you’d describe as *geometry you can touch*. You care because most interactive 3D models and many visualization outputs end up as `vtkPolyData`.

**Do:** use `vtkPolyData` for surfaces, contour results, boundaries, and models that are primarily skins/shells.
**Don’t:** store volumetric cells here, `vtkPolyData` isn’t meant for full 3D cell types like tetrahedra or hexahedra.

* **Description**: Represents a dataset comprising points, vertices, lines, polygons, and triangle strips.
* **Applications**: Ideal for surface meshes and 3D models.
* **Examples**: Models of 3D objects (e.g., vehicles, buildings), terrain surfaces (mountains, valleys).

![poly\_data](https://github.com/djeada/VTK-Examples/assets/37275728/4c642d3b-ecd0-4397-9bdb-2e7a89095416)

### vtkStructuredGrid

A `vtkStructuredGrid` is where structure meets flexibility. The grid topology is still orderly (i,j,k indexing), but the points can warp through space, making it perfect for curvilinear coordinate systems and simulation meshes that bend around objects.

**Do:** use this when you have a logical 3D grid layout but the geometry isn’t axis-aligned.
**Don’t:** use it when connectivity isn’t grid-like. If neighbors aren’t predictable via indices, you’re heading into unstructured territory.

* Curvilinear grid maintaining fixed topology.
* Best for data on curvilinear coordinate systems.
* Fluid flow simulations around objects, finite volume method simulations with a fixed cell division.

![structured\_grid](https://github.com/djeada/VTK-Examples/assets/37275728/c55c5b5a-e6ae-45bc-aa3a-e591eda64e1d)

### vtkUnstructuredGrid

`vtkUnstructuredGrid` is VTK’s “anything goes” dataset. It can mix cell types and represent very complex geometries. The tradeoff is that you have to explicitly store connectivity, which costs memory and computation, but you gain freedom.

This is the one you choose when the world stops being neat: adaptive meshes, FEM simulations, irregular domains, complex anatomy, fractured geology, you name it.

**Do:** use it when the mesh is irregular, adaptive, or mixed-cell.
**Don’t:** default to it “just in case.” If your data is structured, structured grids will usually be faster and lighter.

* An irregular grid with flexible topology, capable of containing various geometric primitives.
* Suitable for complex geometries and adaptive meshes.
* Finite element method simulations with complex, dynamic domains, models of intricate 3D geometries (human brain, aircraft wings).

![unstructured\_grid](https://github.com/djeada/VTK-Examples/assets/37275728/aa34d289-5cc4-4611-bae3-368b3ac49ac4)

### Structured vs Unstructured Grids

This is one of the most practical forks in the road. If you remember nothing else, remember this: **structured grids buy you speed and simplicity**, unstructured grids buy you **geometric freedom**. Your “best” choice is usually the simplest structure that still represents your data truthfully.

**Do:** pick structured when the indexing is predictable and the domain is regular-ish.
**Don’t:** pick unstructured just because it sounds more powerful, it’s powerful, but you pay for it.

|                     | **Structured Grids**                                   | **Unstructured Grids**                                     |
| ------------------- | ------------------------------------------------------ | ---------------------------------------------------------- |
| **Characteristics** | Regular grids with fixed topology                      | Irregular grids with flexible topology                     |
| **Examples**        | vtkImageData, vtkRectilinearGrid, vtkStructuredGrid    | vtkUnstructuredGrid                                        |
| **Advantages**      | Easier to work with, more memory-efficient             | Highly versatile in handling complex shapes                |
| **Disadvantages**   | Limited flexibility in representing complex geometries | Less memory-efficient and can be computationally intensive |

### Multiblock Dataset

Once your datasets get big, or naturally split into parts, you’ll want a container that preserves that structure. `vtkMultiBlockDataSet` lets you keep data organized the way *you* think about it: domains, components, LOD levels, timesteps, regions, etc. This matters because organization isn’t just neatness, it affects pipeline clarity, debugging, and how easily you can process or render subsets.

**Do:** use multiblock when your data is naturally “many datasets,” especially for multi-domain simulations or grouped assets.
**Don’t:** flatten everything into one giant unstructured grid unless you truly need to, keeping blocks separate can make processing and rendering more manageable.

* A collection of multiple datasets, organized hierarchically.
* Handles complex, hierarchical data structures effectively.
* Multi-domain simulations (each domain as a separate dataset), multi-resolution data (datasets at different detail levels).

![multiblock\_dataset](https://github.com/djeada/VTK-Examples/assets/37275728/7a373425-41b2-4a62-bf4b-7f40d40ba04a)

A `vtkMultiBlockDataSet` organizes datasets (or blocks) hierarchically. Each block can be a vtkDataSet subclass or another composite dataset, allowing versatile dataset organization. For instance, a vtkMultiBlockDataSet might contain blocks like vtkPolyData, vtkStructuredGrid, and even another vtkMultiBlockDataSet. This structure is invaluable in large-scale simulations where data is segmented into blocks representing different simulation areas or system components, enabling VTK to manage complex, multi-part datasets coherently.

### Choosing the Right Data Structure

This is where everything clicks: your data structure choice is not “just a container decision.” It’s a commitment that affects performance, tooling, and what operations feel natural. If you match the structure to the data’s *true nature*, VTK feels effortless. If you mismatch it, you’ll constantly translate, patch, and second-guess.

**Do:** choose the least complex option that preserves meaning (geometry + topology + resolution).
**Don’t:** optimize prematurely by picking a complicated type “for future-proofing.” Future-proofing usually means choosing something accurate and convertible, not something maximal.

* Consider the geometry, topology, and resolution of your data
* The choice of data structure affects memory usage, processing time, and ease of use
* Some data structures can be converted to others, but not always without loss of information

| VTK Data Object               | Purpose and Use Cases                                                                     | Complexity                                                                                      | Flexibility                                                                                  | Efficiency                                                                                                      |
| ----------------------------- | ----------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- |
| PolyData                      | Ideal for representing 3D surfaces and complex shapes like sculptures or terrain models.  | Moderate: Comprised of polygons and polylines, allowing for manipulation and detailing.         | High: Excellent for complex, arbitrary shapes including intricate surface details.           | High: Efficiently stores data specific to the shape, omitting unnecessary information.                          |
| Structured Grids              | Suitable for 3D grid-based data like volumetric images or regular spatial data.           | Low: Simple, regular structure with implicit connectivity, making it straightforward to handle. | Low: Confined to regular grid shapes, limiting its application to uniformly structured data. | High for regular data: Efficiently handles uniformly spaced data but less suitable for irregular datasets.      |
| Unstructured Grids            | Used for complex, irregular shapes like biological structures or geological formations.   | High: Requires detailed definition of point connectivity, accommodating intricate geometries.   | High: Extremely versatile, capable of representing any shape, be it regular or irregular.    | Low to Moderate: Efficiency varies with the shape's complexity; generally less efficient than structured grids. |
| Structured Points (ImageData) | Perfect for regularly spaced data such as medical imaging (CT, MRI scans) or 3D textures. | Low: Simple and regular, with implicit connectivity for easy navigation.                        | Low: Restricted to regular grid shapes, ideal for uniformly spaced data.                     | Very High: Optimal for regular data thanks to implicit connectivity and uniform structure.                      |
| Rectilinear Grids             | Appropriate for data with variable resolution, like atmospheric or oceanic datasets.      | Moderate: Allows variable spacing while maintaining a regular grid structure.                   | Moderate: More adaptable than regular grids but less so than unstructured grids.             | High: More efficient than unstructured grids for data that aligns with its variable but regular structure.      |
