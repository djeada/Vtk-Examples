import vtk


class Button2D:
    def __init__(self, interactor, label="Button", pos=(100, 100), size=(100, 40)):
        self.interactor = interactor
        self.label = label
        self.pos = pos
        self.size = size

        # Create the 2D button representation
        self.button_actor = vtk.vtkTexturedButtonRepresentation2D()
        self.button_actor.SetNumberOfStates(1)
        self.button_actor.SetPlaceFactor(1)
        self.button_actor.PlaceWidget(
            [pos[0], pos[0] + size[0], pos[1], pos[1] + size[1], 0, 0]
        )

        # Create a texture for the button
        self.create_button_texture()

        # Create the button widget
        self.button_widget = vtk.vtkButtonWidget()
        self.button_widget.SetInteractor(interactor)
        self.button_widget.SetRepresentation(self.button_actor)
        self.button_widget.On()

        # Create the text label
        self.create_button_label()

    def create_button_texture(self):
        # Create an image for the button texture
        image = vtk.vtkImageData()
        image.SetDimensions(self.size[0], self.size[1], 1)
        image.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 3)

        # Fill the image with a color
        for y in range(self.size[1]):
            for x in range(self.size[0]):
                pixel = [180, 180, 180]  # Light gray color
                image.SetScalarComponentFromFloat(x, y, 0, 0, pixel[0])
                image.SetScalarComponentFromFloat(x, y, 0, 1, pixel[1])
                image.SetScalarComponentFromFloat(x, y, 0, 2, pixel[2])

        # Set the image directly as the button texture
        self.button_actor.SetButtonTexture(0, image)

    def create_button_label(self):
        # Create a text actor
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(self.label)

        # Adjust font size and alignment
        text_property = text_actor.GetTextProperty()
        text_property.SetFontSize(
            min(self.size) // 2
        )  # Adjust font size based on button size
        text_property.SetColor(0, 0, 0)  # Black text
        text_property.SetJustificationToCentered()
        text_property.SetVerticalJustificationToCentered()

        # Estimate the vertical centering position
        approx_text_height = text_property.GetFontSize() * 0.5
        y_position = self.pos[1] + (self.size[1] - approx_text_height) / 2

        # Position the text in the center of the button
        text_actor.SetPosition(self.pos[0] + self.size[0] / 2, y_position)

        # Add the text actor to the renderer
        self.interactor.GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor(
            text_actor
        )

    def add_click_callback(self, callback):
        # Add click event listener
        def on_button_press(widget, event_string):
            callback()

        self.button_widget.AddObserver("StateChangedEvent", on_button_press)
