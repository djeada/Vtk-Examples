import sys

import vtk
from PyQt6 import QtWidgets, QtCore
from PyQt6.QtWidgets import (
    QWidget,
    QMainWindow,
    QVBoxLayout,
    QSlider,
    QPushButton,
    QColorDialog,
)
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkFiltersSources import vtkSphereSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer


class CustomActor(vtkActor):
    def __init__(self, source=vtkSphereSource()):
        self.source = source
        self.color = (0.1, 0.8, 0.5)  # Default color
        self.update_mapper()

    def update_mapper(self):
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(self.source.GetOutputPort())
        self.SetMapper(mapper)
        self.SetPosition(0.2, 0.2, 0.6)
        self.GetProperty().SetColor(self.color)

    def update_radius(self, radius):
        self.source.SetRadius(radius)
        self.source.Update()
        self.update_mapper()

    def update_color(self, color):
        self.color = (color.red() / 255, color.green() / 255, color.blue() / 255)
        self.GetProperty().SetColor(self.color)


class VtkDisplay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        self.vtk_display = QVTKRenderWindowInteractor(self)
        layout.addWidget(self.vtk_display)

        self.slider = QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setMinimum(1)
        self.slider.setMaximum(10)
        self.slider.setValue(1)
        layout.addWidget(self.slider)

        self.color_button = QPushButton("Change color")
        layout.addWidget(self.color_button)

        self.renderer = vtkRenderer()
        self.vtk_display.GetRenderWindow().AddRenderer(self.renderer)
        self.interactor = self.vtk_display.GetRenderWindow().GetInteractor()

        self.actor = CustomActor()
        self.renderer.AddActor(self.actor)

        self.slider.valueChanged.connect(self.change_radius)
        self.color_button.clicked.connect(self.change_color)

        self.renderer.ResetCamera()
        self.interactor.Initialize()

    def change_radius(self, value):
        self.actor.update_radius(value)
        self.vtk_display.GetRenderWindow().Render()

    def change_color(self):
        color = QColorDialog.getColor()
        if color.isValid():
            self.actor.update_color(color)
            self.vtk_display.GetRenderWindow().Render()


if __name__ == "__main__":
    app = QtWidgets.QApplication([])

    window = QMainWindow()
    widget = VtkDisplay()
    window.setCentralWidget(widget)
    window.show()
    sys.exit(app.exec())
