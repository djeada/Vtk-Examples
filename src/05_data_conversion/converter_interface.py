from abc import ABC, abstractmethod


class Converter(ABC):
    @abstractmethod
    def convert(self, input_filename: str, output_filename: str):
        pass
