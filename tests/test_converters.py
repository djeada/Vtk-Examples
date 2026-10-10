import importlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import vtk


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


def converter(module, name):
    return getattr(importlib.import_module(f"src.05_data_conversion.{module}"), name)()


def read_data(path):
    readers = {
        ".obj": vtk.vtkOBJReader,
        ".stl": vtk.vtkSTLReader,
        ".ply": vtk.vtkPLYReader,
        ".vtp": vtk.vtkXMLPolyDataReader,
        ".vtu": vtk.vtkXMLUnstructuredGridReader,
        ".vtk": vtk.vtkGenericDataObjectReader,
        ".vtm": vtk.vtkXMLMultiBlockDataReader,
    }
    reader = readers[path.suffix]()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


class ConverterTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.output = Path(self.directory.name)

    def test_all_reader_writer_pairs_round_trip_geometry(self):
        cases = (
            (
                "stl_obj",
                "StlToObjConverter",
                "ObjToStlConverter",
                "stls/cube.stl",
                ".obj",
            ),
            (
                "obj_vtp",
                "ObjToVtpConverter",
                "VtpToObjConverter",
                "objs/cube.obj",
                ".vtp",
            ),
            (
                "ply_obj",
                "ObjToPlyConverter",
                "PlyToObjConverter",
                "objs/cube.obj",
                ".ply",
            ),
            (
                "ply_stl",
                "StlToPlyConverter",
                "PlyToStlConverter",
                "stls/cube.stl",
                ".ply",
            ),
            (
                "stl_vtk",
                "StlToVtkConverter",
                "VtkToStlConverter",
                "stls/cube.stl",
                ".vtk",
            ),
            (
                "stl_vtp",
                "StlToVtpConverter",
                "VtpToStlConverter",
                "stls/cube.stl",
                ".vtp",
            ),
            (
                "vtk_obj",
                "VtkToObjConverter",
                "ObjToVtkConverter",
                "vtks/cube_polydata.vtk",
                ".obj",
            ),
            (
                "vtk_vtp",
                "VtkToVtpConverter",
                "VtpToVtkConverter",
                "vtks/cube_polydata.vtk",
                ".vtp",
            ),
            (
                "vtk_vtu",
                "VtkToVtuConverter",
                "VtuToVtkConverter",
                "vtks/grid_of_triangles.vtk",
                ".vtu",
            ),
            (
                "vtk_vtm",
                "VtkToVtmConverter",
                "VtmToVtkConverter",
                "vtks/grid_of_triangles.vtk",
                ".vtm",
            ),
        )
        for module, forward, backward, relative_source, extension in cases:
            with self.subTest(module=module):
                source = DATA_DIR / relative_source
                intermediate = self.output / (module + extension)
                restored = self.output / (module + "_restored" + source.suffix)
                converter(module, forward).convert(source, intermediate)
                converter(module, backward).convert(intermediate, restored)
                original = read_data(source)
                result = read_data(restored)
                self.assertEqual(
                    result.GetNumberOfPoints(), original.GetNumberOfPoints()
                )
                # Text formats may round coordinate precision during export.
                for actual, expected in zip(result.GetBounds(), original.GetBounds()):
                    self.assertAlmostEqual(actual, expected, places=6)
                self.assertGreater(result.GetNumberOfCells(), 0)

    def test_legacy_volume_can_be_exported_as_a_surface(self):
        volume = vtk.vtkImageData()
        volume.SetDimensions(2, 2, 2)
        append = vtk.vtkAppendFilter()
        append.SetInputData(volume)
        append.Update()
        source = self.output / "volume.vtk"
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(str(source))
        writer.SetInputData(append.GetOutput())
        self.assertEqual(writer.Write(), 1)
        for module, name, suffix in (
            ("stl_vtk", "VtkToStlConverter", ".stl"),
            ("vtk_obj", "VtkToObjConverter", ".obj"),
            ("vtk_vtp", "VtkToVtpConverter", ".vtp"),
        ):
            with self.subTest(suffix=suffix):
                target = self.output / ("surface" + suffix)
                converter(module, name).convert(source, target)
                surface = read_data(target)
                self.assertIsInstance(surface, vtk.vtkPolyData)
                self.assertGreaterEqual(surface.GetNumberOfPolys(), 6)
                self.assertEqual(surface.GetBounds(), volume.GetBounds())

    def test_vtu_rejects_polydata_before_writing(self):
        target = self.output / "invalid.vtu"
        with self.assertRaisesRegex(
            ValueError, "Expected vtkUnstructuredGrid.*vtkPolyData"
        ):
            converter("vtk_vtu", "VtkToVtuConverter").convert(
                DATA_DIR / "vtks/cube_polydata.vtk", target
            )
        self.assertFalse(target.exists())

    def test_invalid_paths_and_extensions(self):
        instance = converter("stl_obj", "StlToObjConverter")
        with self.assertRaises(ValueError):
            instance.convert("", self.output / "x.obj")
        with self.assertRaisesRegex(ValueError, "extensions"):
            instance.convert(DATA_DIR / "stls/cube.stl", self.output / "x.vtu")
        with self.assertRaises(FileNotFoundError):
            instance.convert(self.output / "missing.stl", self.output / "x.obj")
        with self.assertRaisesRegex(FileNotFoundError, "Output directory"):
            instance.convert(
                DATA_DIR / "stls/cube.stl", self.output / "missing" / "x.obj"
            )

    def test_corrupt_input_does_not_report_success(self):
        source = self.output / "corrupt.vtp"
        source.write_text("not a VTK file")
        target = self.output / "invalid.obj"
        with self.assertRaises((OSError, ValueError)):
            converter("obj_vtp", "VtpToObjConverter").convert(source, target)
        self.assertFalse(target.exists())

    def test_empty_stl_is_rejected(self):
        source = self.output / "empty.stl"
        source.write_text("solid empty\nendsolid empty\n")
        with self.assertRaises((OSError, ValueError)):
            converter("stl_obj", "StlToObjConverter").convert(
                source, self.output / "empty.obj"
            )

    def test_line_data_cannot_be_written_as_stl(self):
        with self.assertRaisesRegex(ValueError, "STL requires a surface"):
            converter("stl_vtk", "VtkToStlConverter").convert(
                DATA_DIR / "vtks/polydata_line.vtk", self.output / "line.stl"
            )

    def test_writer_failure_raises_even_if_output_already_exists(self):
        class FailedWriter(vtk.vtkOBJWriter):
            def Write(self):
                return 0

        instance = converter("stl_obj", "StlToObjConverter")
        instance.writer_type = FailedWriter
        target = self.output / "existing.obj"
        target.write_text("old output")
        with self.assertRaisesRegex(OSError, "Failed to process"):
            instance.convert(DATA_DIR / "stls/cube.stl", target)

    def test_missing_first_multiblock_is_rejected(self):
        instance = converter("vtk_vtm", "VtmToVtkConverter")
        multiblock = vtk.vtkMultiBlockDataSet()
        with self.assertRaisesRegex(ValueError, "No blocks"):
            instance.prepare_data(multiblock)
        multiblock.SetNumberOfBlocks(1)
        with self.assertRaisesRegex(ValueError, "first VTM block is empty"):
            instance.prepare_data(multiblock)

    def test_multiblock_export_warns_when_dropping_blocks(self):
        multiblock = vtk.vtkMultiBlockDataSet()
        for index in range(2):
            sphere = vtk.vtkSphereSource()
            sphere.Update()
            multiblock.SetBlock(index, sphere.GetOutput())
        with self.assertLogs(level="WARNING") as logs:
            result = converter("vtk_vtm", "VtmToVtkConverter").prepare_data(multiblock)
        self.assertIn("Only the first block", logs.output[0])
        self.assertGreater(result.GetNumberOfPoints(), 0)


if __name__ == "__main__":
    unittest.main()
