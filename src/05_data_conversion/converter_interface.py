from abc import ABC, abstractmethod
import logging
from pathlib import Path

import vtk


class Converter(ABC):
    """Abstract base class for file format converters."""

    @abstractmethod
    def convert(self, input_filename: str, output_filename: str) -> None:
        """Convert a file from one format to another.

        Args:
            input_filename: Path to the input file.
            output_filename: Path to the output file.
        """
        pass


class VTKConverter(Converter):
    """Shared validation and error handling for the concrete reader/writer pairs.

    VTK often emits ErrorEvent instead of raising a Python exception. Check both
    those events and error codes so a failed conversion cannot report success.
    Subclasses declare their formats and override prepare_data when needed.
    """

    input_extension = ""
    output_extension = ""
    reader_type = None
    writer_type = None
    output_type = "vtkDataSet"
    extract_surface = False

    @staticmethod
    def _execute(algorithm, operation, filename):
        errors = []

        @vtk.calldata_type(vtk.VTK_STRING)
        def on_error(caller, event, message):
            errors.append(message.strip() if message else "VTK pipeline error")

        # Type and pipeline errors can be emitted by the executive itself.
        objects = (algorithm, algorithm.GetExecutive())
        observers = [obj.AddObserver("ErrorEvent", on_error) for obj in objects]
        try:
            result = operation()
        finally:
            for obj, observer in zip(objects, observers):
                obj.RemoveObserver(observer)
        if errors or algorithm.GetErrorCode() or result == 0:
            detail = (
                errors[0]
                if errors
                else (
                    f"VTK error code {algorithm.GetErrorCode()}"
                    if algorithm.GetErrorCode()
                    else "VTK operation returned failure"
                )
            )
            raise OSError(f"Failed to process {filename}: {detail}")

    def prepare_data(self, data):
        """Adapt input geometry to the target format before writing it."""
        if self.extract_surface and not data.IsA("vtkPolyData"):
            if not data.IsA("vtkDataSet"):
                raise ValueError(
                    f"Cannot extract a surface from {data.GetClassName()}."
                )
            surface = vtk.vtkDataSetSurfaceFilter()
            surface.SetInputData(data)
            self._execute(surface, surface.Update, "surface geometry")
            data = surface.GetOutput()
        return data

    def convert(self, input_filename, output_filename):
        """Convert str or pathlib paths; raise on invalid data or failed I/O.

        The output directory must already exist. Existing output files are
        overwritten by VTK, matching the original example behavior.
        """
        if not input_filename or not output_filename:
            raise ValueError("Input and output filenames must be provided.")
        source = Path(input_filename)
        destination = Path(output_filename)
        if (
            source.suffix.lower() != self.input_extension
            or destination.suffix.lower() != self.output_extension
        ):
            raise ValueError(
                f"Invalid file extensions. Expected '{self.input_extension}' and '{self.output_extension}'."
            )
        if not source.is_file():
            raise FileNotFoundError(f"Input file does not exist: {source}")
        if not destination.parent.is_dir():
            raise FileNotFoundError(
                f"Output directory does not exist: {destination.parent}"
            )
        if destination.is_dir():
            raise IsADirectoryError(f"Output path is a directory: {destination}")

        reader = self.reader_type()
        reader.SetFileName(str(source))
        self._execute(reader, reader.Update, source)
        data = reader.GetOutput()
        if data is None:
            raise ValueError(f"Input contains no dataset: {source}")
        data = self.prepare_data(data)
        if data is None or not data.IsA(self.output_type):
            actual = data.GetClassName() if data is not None else "no data"
            raise ValueError(
                f"Expected {self.output_type} for {destination.suffix}, got {actual}."
            )
        if data.IsA("vtkDataSet") and data.GetNumberOfPoints() == 0:
            raise ValueError(f"Input contains no points: {source}")
        if (
            self.output_extension == ".stl"
            and data.GetNumberOfPolys() + data.GetNumberOfStrips() == 0
        ):
            raise ValueError("STL requires a surface with polygons or triangle strips.")

        writer = self.writer_type()
        writer.SetFileName(str(destination))
        writer.SetInputData(data)
        self._execute(writer, writer.Write, destination)
        if not destination.is_file() or destination.stat().st_size == 0:
            raise OSError(f"VTK did not write any data to {destination}")
        logging.info("Conversion successful: %s to %s", source, destination)
