import sys

import vtk
from PyQt6 import QtWidgets
from PyQt6.QtWidgets import QWidget, QMainWindow
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer


class CustomActor(vtkActor):
    def __init__(self, source=vtkSphereSource()):
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        self.SetMapper(mapper)
        self.SetPosition(0.2, 0.2, 0.6)
        self.GetProperty().SetColor(0.1, 0.8, 0.5)


class VtkDisplay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QtWidgets.QVBoxLayout(self)
        self.vtk_display = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.vtk_display)

        self.renderer = vtkRenderer()
        self.vtk_display.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_display.GetRenderWindow().GetInteractor()

        self.renderer.AddActor(CustomActor())

        self.renderer.ResetCamera()
        self.interactor.Initialize()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    window = QMainWindow()
    widget = VtkDisplay()
    window.setCentralWidget(widget)
    window.show()
    sys.exit(app.exec())
