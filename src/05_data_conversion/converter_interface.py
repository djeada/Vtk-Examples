from abc import ABC, abstractmethod


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
