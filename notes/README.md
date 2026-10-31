# VTK Notes

These notes explain VTK's data model, pipelines, visualization techniques, and integrations. Start with the foundations, then follow the topics relevant to your project.

## Reading order

1. [Data types and structures](https://github.com/djeada/Vtk-Examples/blob/main/notes/01_data_types_and_structures.md) — geometry, topology, and dataset types.
2. [Filters and algorithms](https://github.com/djeada/Vtk-Examples/blob/main/notes/02_filters_and_algorithms.md) — sources, connectivity, and pipeline transformations.
3. [Input and output](https://github.com/djeada/Vtk-Examples/blob/main/notes/03_input_and_output.md) — choosing file formats, readers, and writers.
4. [Visualization techniques](https://github.com/djeada/Vtk-Examples/blob/main/notes/04_visualization_techniques.md) — volume rendering, flow, glyphs, and contours.
5. [Interactivity](https://github.com/djeada/Vtk-Examples/blob/main/notes/05_interactivity.md) — picking, cutting, widgets, and events.
6. [Animations](https://github.com/djeada/Vtk-Examples/blob/main/notes/06_animations.md) — keyframes, temporal data, and export.
7. [Performance and parallelism](https://github.com/djeada/Vtk-Examples/blob/main/notes/07_performance_optimization_and_parallelism.md) — level of detail, culling, and parallel processing.
8. [Integration with other tools](https://github.com/djeada/Vtk-Examples/blob/main/notes/08_integration_with_other_tools.md) — ParaView, VisIt, ITK, Python, and HPC.
9. [Custom filters and algorithms](https://github.com/djeada/Vtk-Examples/blob/main/notes/09_custom_filters_and_algorithms.md) — extending VTK's pipeline.
10. [End-to-end case study](https://github.com/djeada/Vtk-Examples/blob/main/notes/10_end_to_end_case_study.md) — build, annotate, export, reload and validate a mesh.

## Using these notes

The guides combine concepts with illustrative code snippets. Some snippets assume an existing VTK pipeline, dataset, or optional integration; they are not all standalone programs. For runnable Python examples and installation instructions, see the [project README](https://github.com/djeada/Vtk-Examples/blob/main/README.md) and [example catalog](https://github.com/djeada/Vtk-Examples/tree/main/src).

**Suggested paths:** For a first visualization, read chapters 1–4 and try the end-to-end case study; for interactive applications, continue with 5–6; for large datasets or new pipeline components, use 7–9.
