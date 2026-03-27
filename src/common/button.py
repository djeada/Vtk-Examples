import vtk


class Button2D:
    def __init__(
        self,
        interactor,
        label="Button",
        pos=(100, 100),
        size=(100, 40),
        background_color=(180, 180, 180),
        text_color=(0.0, 0.0, 0.0),
        font_size=None,
    ):
        self.interactor = interactor
        self.label = label
        self.pos = pos
        self.size = size
        self.background_color = background_color
        self.text_color = text_color
        self.font_size = font_size
        self.text_actor = None

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
                pixel = self.background_color
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
            self.font_size if self.font_size is not None else min(self.size) // 2
        )
        text_property.SetColor(*self.text_color)
        text_property.SetJustificationToCentered()
        text_property.SetVerticalJustificationToCentered()
        text_property.SetBold(True)

        # Estimate the vertical centering position
        approx_text_height = text_property.GetFontSize() * 0.5
        y_position = self.pos[1] + (self.size[1] - approx_text_height) / 2

        # Position the text in the center of the button
        text_actor.SetPosition(self.pos[0] + self.size[0] / 2, y_position)
        self.text_actor = text_actor

        # Add the text actor to the renderer
        self.interactor.GetRenderWindow().GetRenderers().GetFirstRenderer().AddActor2D(
            text_actor
        )

    def add_click_callback(self, callback):
        # Add click event listener
        def on_button_press(widget, event_string):
            callback()

        self.button_widget.AddObserver("StateChangedEvent", on_button_press)
