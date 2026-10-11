"""Cinematic vertical reel layout: a floating, glowing scene with kinetic type.

Everything is designed on a 1080 × 1920 grid and scaled to the output size.
Static artwork (backdrop, legend, chart, captions) is rasterized once; each
frame only composites the live scene, its bloom, and animated text sprites.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


INK = (244, 247, 252)
SOFT = (196, 207, 223)
MUTED = (128, 145, 170)
PANEL = (16, 24, 40)

_WEIGHTS = {
    "regular": ("Inter-Regular.otf", False),
    "medium": ("Inter-Medium.otf", False),
    "semibold": ("Inter-SemiBold.otf", True),
    "bold": ("InterDisplay-Bold.otf", True),
    "black": ("InterDisplay-ExtraBold.otf", True),
}

# Layout on the 1080 × 1920 design grid.
MARGIN = 72
PROGRESS_Y = 54
BRAND_Y = 92
CHIP_Y = 156
HEADLINE_Y = 236
SCENE_BOX = (0, 640, 1080, 1400)
CARD_BOX = (MARGIN, 1530, 1080 - MARGIN, 1752)
NOTE_Y = 1784


def scene_size(width):
    """Pixel size of the scene box for a reel of the given width."""
    scale = width / 1080
    return (
        round((SCENE_BOX[2] - SCENE_BOX[0]) * scale),
        round((SCENE_BOX[3] - SCENE_BOX[1]) * scale),
    )


@lru_cache(maxsize=96)
def font(size, weight="regular"):
    if weight is True:
        weight = "bold"
    elif weight is False:
        weight = "regular"
    name, bold = _WEIGHTS[weight]
    candidates = [Path("/usr/share/fonts/opentype/inter") / name]
    from matplotlib import get_data_path

    candidates.append(
        Path(get_data_path())
        / "fonts"
        / "ttf"
        / ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")
    )
    for candidate in candidates:
        if candidate.is_file():
            return ImageFont.truetype(str(candidate), size=max(1, round(size)))
    return ImageFont.load_default(size=max(1, round(size)))


def wrap(text, face, width):
    lines = []
    line = ""
    for word in text.split():
        trial = f"{line} {word}".strip()
        if line and face.getlength(trial) > width:
            lines.append(line)
            line = word
        else:
            line = trial
        # Long words in custom stories must fit the mobile caption column too.
        while face.getlength(line) > width:
            count = 1
            while count < len(line) and face.getlength(line[: count + 1]) <= width:
                count += 1
            lines.append(line[:count])
            line = line[count:]
    if line:
        lines.append(line)
    return lines


def fit_lines(text, size, weight, width, max_lines, minimum=10):
    """Wrap ``text``, shrinking the type until it fits in ``max_lines``."""
    face = font(size, weight)
    lines = wrap(text, face, width)
    while len(lines) > max_lines and size > minimum:
        size -= 1
        face = font(size, weight)
        lines = wrap(text, face, width)
    return lines, face, size


def _ease_out(value):
    value = min(1.0, max(0.0, value))
    return 1 - (1 - value) ** 3


def _smootherstep(value):
    value = min(1.0, max(0.0, value))
    return value * value * value * (value * (value * 6 - 15) + 10)


class Sprite:
    """Premultiplied RGBA pixels placed at a design position."""

    def __init__(self, image, x, y):
        pixels = np.asarray(image.convert("RGBA"), dtype=np.float32) / 255
        self.alpha = pixels[..., 3]
        self.rgb = pixels[..., :3] * self.alpha[..., None]
        self.x, self.y = int(round(x)), int(round(y))

    def draw(self, frame, opacity=1.0, dx=0, dy=0):
        if opacity <= 0.003:
            return
        x, y = self.x + int(round(dx)), self.y + int(round(dy))
        height, width = self.alpha.shape
        top, left = max(0, y), max(0, x)
        bottom = min(frame.shape[0], y + height)
        right = min(frame.shape[1], x + width)
        if bottom <= top or right <= left:
            return
        alpha = self.alpha[top - y : bottom - y, left - x : right - x, None]
        rgb = self.rgb[top - y : bottom - y, left - x : right - x]
        region = frame[top:bottom, left:right]
        region *= 1 - alpha * opacity
        region += rgb * opacity


def _text_image(text, face, color, tracking=0.0):
    ascent, descent = face.getmetrics()
    pad = max(2, face.size // 8)
    if tracking:
        width = sum(face.getlength(char) for char in text) + tracking * max(0, len(text) - 1)
    else:
        width = face.getlength(text)
    image = Image.new("RGBA", (int(width) + 2 * pad, ascent + descent + 2 * pad), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    if tracking:
        x = pad
        for char in text:
            draw.text((x, pad), char, font=face, fill=color)
            x += face.getlength(char) + tracking
    else:
        draw.text((pad, pad), text, font=face, fill=color)
    return image, pad


def _gradient_text(text, face, top, bottom):
    image, pad = _text_image(text, face, (255, 255, 255))
    height = image.height
    ramp = np.linspace(0, 1, height)[:, None, None]
    color = np.array(top, np.float32) * (1 - ramp) + np.array(bottom, np.float32) * ramp
    pixels = np.asarray(image, dtype=np.float32)
    pixels[..., :3] = np.broadcast_to(color, pixels[..., :3].shape)
    return Image.fromarray(pixels.astype(np.uint8), "RGBA"), pad


def _lighten(color, amount):
    return tuple(int(c + (255 - c) * amount) for c in color)


class ReelLayout:
    """Six story chapters over a floating scene, with one consistent visual system."""

    def __init__(self, story, model, width, height, duration, colors, legend, scene_height=None):
        self.story, self.model = story, model
        self.width, self.height, self.duration = width, height, duration
        self.s = width / 1080
        self.colors, self.legend = colors, legend
        self.accent = tuple(int(c) for c in story.accent)
        box = [round(v * self.s) for v in SCENE_BOX]
        if scene_height is not None:
            box[3] = box[1] + int(scene_height)
        self.scene_box = tuple(box)
        self.scene_size = (
            self.scene_box[2] - self.scene_box[0],
            self.scene_box[3] - self.scene_box[1],
        )
        # A shorter scene hands its spare height to the legend and card below.
        self.legend_y = self.scene_box[3] + round(18 * self.s)
        self.card_box = (
            round(CARD_BOX[0] * self.s),
            self.legend_y + round(112 * self.s),
            round(CARD_BOX[2] * self.s),
            round(CARD_BOX[3] * self.s),
        )
        self.chart_box = self.card_box
        self.background = self._background()
        self.static = [self._brand(), self._legend_row(), self._footnote()]
        self.card = self._card_panel()
        self.chart_frame, self.chart_curve, self.chart_head = (
            self._chart() if model.chart else (None, None, None)
        )
        self.note_card = None if model.chart or model.solid is not None else self._note_card()
        self.chapters = [self._chapter(index) for index in range(len(story.beats))]
        self.detail_label = Sprite(
            _text_image("CLOSE-UP · SAME COLOR SCALE", font(17 * self.s, "semibold"), MUTED, 2 * self.s)[0],
            self.card_box[0] + 22 * self.s,
            self.card_box[1] + 18 * self.s,
        )

    # ------------------------------------------------------------------ static art

    def _background(self):
        s, w, h = self.s, self.width, self.height
        yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
        v = yy / h
        top = np.array((8, 12, 24), np.float32)
        middle = np.array((11, 17, 33), np.float32)
        bottom = np.array((3, 5, 11), np.float32)
        ramp = np.where(v < 0.45, v / 0.45, 1 - (v - 0.45) / 0.55)[..., None]
        base = np.where(
            (v < 0.45)[..., None],
            top * (1 - ramp) + middle * ramp,
            bottom * (1 - ramp) + middle * ramp,
        )
        cx, cy = (self.scene_box[0] + self.scene_box[2]) / 2, (self.scene_box[1] + self.scene_box[3]) / 2
        glow = np.exp(-(((xx - cx) / (620 * s)) ** 2 + ((yy - cy) / (560 * s)) ** 2))
        accent = np.array(self.accent, np.float32)
        base += glow[..., None] * (0.15 * accent + np.array((4, 10, 22), np.float32))
        cool = np.exp(-(((xx - 160 * s) / (520 * s)) ** 2 + ((yy - 300 * s) / (420 * s)) ** 2))
        base += cool[..., None] * np.array((10, 22, 40), np.float32) * 0.6
        # A faint measurement grid that dissolves away from the scene.
        spacing = 54 * s
        lines = (
            (np.abs((xx + spacing / 2) % spacing - spacing / 2) < 0.6)
            | (np.abs((yy + spacing / 2) % spacing - spacing / 2) < 0.6)
        ).astype(np.float32)
        base += lines[..., None] * glow[..., None] * 9
        radius = ((xx / w - 0.5) ** 2 + (v - 0.5) ** 2 * 0.6) * 2
        base *= (1 - 0.42 * np.clip(radius, 0, 1))[..., None]
        grain = np.random.default_rng(3).normal(0, 1.4, (h, w, 1)).astype(np.float32)
        return np.clip(base + grain, 0, 255) / 255

    def _brand(self):
        s = self.s
        layer = Image.new("RGBA", (self.width, round(140 * s)), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        dot = 12 * s
        y = BRAND_Y * s
        draw.ellipse((MARGIN * s, y + 4 * s, MARGIN * s + dot, y + 4 * s + dot), fill=self.accent)
        face = font(20 * s, "semibold")
        image, pad = _text_image("VTK  FIELD NOTES", face, SOFT, 3.2 * s)
        layer.alpha_composite(image, (round(MARGIN * s + dot + 14 * s - pad), round(y - pad - 2 * s)))
        right = "SCIENCE IN 60s" if round(self.duration) == 60 else "SCIENCE IN MOTION"
        image, pad = _text_image(right, font(20 * s, "medium"), MUTED, 3.2 * s)
        layer.alpha_composite(image, (round(self.width - MARGIN * s - image.width + pad), round(y - pad - 2 * s)))
        return Sprite(layer, 0, 0)

    def _legend_row(self):
        s = self.s
        top = self.legend_y
        layer = Image.new("RGBA", (self.width, round(100 * s)), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        x0 = MARGIN * s
        label = self.model.field_name
        if self.legend is not None:
            label += f"  ·  {self.legend[1]}"
        image, pad = _text_image(label.upper(), font(17 * s, "semibold"), MUTED, 2.4 * s)
        layer.alpha_composite(image, (round(x0 - pad), 0))
        if self.legend is not None:
            (low, high), unit = self.legend
            bar_w, bar_h = round(430 * s), max(2, round(12 * s))
            ramp = Image.fromarray(np.broadcast_to(self.colors[None], (4, 256, 3)).copy()).resize(
                (bar_w, bar_h), Image.BILINEAR
            )
            mask = Image.new("L", (bar_w, bar_h), 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, bar_w - 1, bar_h - 1), radius=bar_h // 2, fill=255)
            layer.paste(ramp, (round(x0), round(34 * s)), mask)
            face = font(21 * s, "medium")
            # Round underflow such as 1e-66 to the zero it represents.
            if abs(low) < 1e-9 * max(1.0, abs(high)):
                low = 0.0
            draw.text((x0, 54 * s), f"{low:.3g}", font=face, fill=SOFT)
            text = f"{high:.3g}"
            draw.text((x0 + bar_w - face.getlength(text), 54 * s), text, font=face, fill=SOFT)
        else:
            draw.text((x0, 38 * s), "Original VTK scene", font=font(26 * s, "medium"), fill=SOFT)
        for index, (name, value) in enumerate(self.model.metrics[:2]):
            mx = (590 + index * 220) * s
            image, pad = _text_image(name.upper(), font(15 * s, "semibold"), MUTED, 2 * s)
            layer.alpha_composite(image, (round(mx - pad), 0))
            lines, face, _ = fit_lines(value, 44 * s, "bold", 205 * s, 1)
            draw.text((mx, 26 * s), lines[0], font=face, fill=_lighten(self.accent, 0.1))
        return Sprite(layer, 0, top)

    def _card_panel(self):
        s = self.s
        x0, y0, x1, y1 = self.card_box
        layer = Image.new("RGBA", (x1 - x0, y1 - y0), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        draw.rounded_rectangle(
            (0, 0, x1 - x0 - 1, y1 - y0 - 1),
            radius=round(26 * s),
            fill=(22, 32, 52, 150),
            outline=(255, 255, 255, 30),
            width=max(1, round(1.5 * s)),
        )
        return Sprite(layer, x0, y0)

    def _note_card(self):
        s = self.s
        x0, y0, x1, y1 = self.card_box
        layer = Image.new("RGBA", (x1 - x0, y1 - y0), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        image, pad = _text_image("WHAT YOU ARE SEEING", font(16 * s, "semibold"), self.accent, 2.4 * s)
        layer.alpha_composite(image, (round(28 * s - pad), round(24 * s - pad)))
        lines, face, size = fit_lines(self.story.model_note, 30 * s, "medium", (x1 - x0) - 56 * s, 4)
        y = 64 * s
        for line in lines:
            draw.text((28 * s, y), line, font=face, fill=SOFT)
            y += size * 1.3
        return Sprite(layer, x0, y0)

    def _footnote(self):
        s = self.s
        layer = Image.new("RGBA", (self.width, round(110 * s)), (0, 0, 0, 0))
        draw = ImageDraw.Draw(layer)
        lines, face, size = fit_lines(self.model.note, 22 * s, "regular", self.width - 2 * MARGIN * s, 2)
        y = 0
        for line in lines:
            draw.text((MARGIN * s, y), line, font=face, fill=MUTED)
            y += size * 1.35
        draw.text(
            (MARGIN * s, y + 6 * s),
            self.model.key.replace("_", " "),
            font=font(18 * s, "medium"),
            fill=(88, 104, 128),
        )
        return Sprite(layer, 0, NOTE_Y * s)

    def _chart(self):
        """Antialiased profile chart: a fixed frame layer plus a revealable curve layer."""
        k = 3  # draw supersampled, then downsample
        s = self.s * k
        x0, y0, x1, y1 = self.card_box
        width, height = (x1 - x0) * k, (y1 - y0) * k
        frame = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        curve = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        fd, cd = ImageDraw.Draw(frame), ImageDraw.Draw(curve)
        chart = self.model.chart
        x = np.asarray(chart["x"], dtype=float)
        series = [np.asarray(values, dtype=float) for values in chart["series"]]
        values = np.concatenate([values[np.isfinite(values)] for values in series])
        low, high = float(values.min()), float(values.max())
        if high == low:
            high = low + 1
        ticks = (low, (low + high) / 2, high)
        pad_y = 0.08 * (high - low)
        low, high = low - pad_y, high + pad_y
        left, right = 92 * s, width - 30 * s
        top, bottom = 62 * s, height - 50 * s
        label_face = font(16 * s, "semibold")
        tick_face = font(17 * s, "medium")
        title, pad = _text_image("COMPUTED PROFILE", label_face, MUTED, 2.4 * s)
        frame.alpha_composite(title, (round(26 * s - pad), round(20 * s - pad)))
        for number, value in enumerate(ticks):
            gy = bottom - (value - low) / (high - low) * (bottom - top)
            fd.line((left, gy, right, gy), fill=(255, 255, 255, 48 if number == 0 else 22), width=max(1, round(1.2 * s)))
            text = f"{value:.3g}"
            fd.text((left - 14 * s - tick_face.getlength(text), gy - 11 * s), text, font=tick_face, fill=MUTED)
        fd.text((left, bottom + 12 * s), f"{x.min():.3g}", font=tick_face, fill=MUTED)
        text = f"{x.max():.3g}"
        fd.text((right - tick_face.getlength(text), bottom + 12 * s), text, font=tick_face, fill=MUTED)
        label = f"{chart['xlabel']}  ·  {chart['ylabel']}"
        fd.text(((left + right) / 2 - tick_face.getlength(label) / 2, bottom + 12 * s), label, font=tick_face, fill=MUTED)
        palette = (_lighten(self.accent, 0.08), (226, 234, 248), (140, 192, 255))
        legend_x = right
        for index in reversed(range(len(series))):
            name = chart["labels"][index]
            face = font(17 * s, "semibold")
            legend_x -= face.getlength(name)
            fd.text((legend_x, 18 * s), name, font=face, fill=palette[index % 3])
            sw = 26 * s
            legend_x -= sw + 10 * s
            fd.line((legend_x, 29 * s, legend_x + sw, 29 * s), fill=palette[index % 3], width=round(3 * s))
            legend_x -= 26 * s

        def project(ys):
            return np.c_[
                left + (x - x.min()) / np.ptp(x) * (right - left),
                bottom - (ys - low) / (high - low) * (bottom - top),
            ]

        main = project(series[0])
        # Soft gradient fill under the primary series.
        fill = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        polygon = [tuple(p) for p in main] + [(main[-1, 0], bottom), (main[0, 0], bottom)]
        mask = Image.new("L", (width, height), 0)
        ImageDraw.Draw(mask).polygon(polygon, fill=255)
        fade = np.clip((bottom - np.arange(height)) / (bottom - top), 0, 1)[:, None] ** 0.8
        alpha = (np.asarray(mask, np.float32) * fade * 0.28).astype(np.uint8)
        fill.paste(Image.new("RGBA", (width, height), self.accent + (255,)), (0, 0), Image.fromarray(alpha))
        curve.alpha_composite(fill)
        glow = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        ImageDraw.Draw(glow).line([tuple(p) for p in main], fill=self.accent + (150,), width=round(9 * s), joint="curve")
        curve.alpha_composite(glow.filter(ImageFilter.GaussianBlur(7 * s)))
        for index, ys in enumerate(series):
            points = [tuple(p) for p in project(ys)]
            color = palette[index % 3]
            if index == 1:
                # Reference is dashed so coincident numerical/analytic curves are distinguishable.
                for start in range(0, len(points) - 1, 6):
                    cd.line(points[start : start + 4], fill=color + (235,), width=round(2.6 * s))
            else:
                cd.line(points, fill=color + (255,), width=round(3.4 * s), joint="curve")
        down = (x1 - x0, y1 - y0)
        frame = frame.resize(down, Image.LANCZOS)
        curve = curve.resize(down, Image.LANCZOS)
        head = main / k
        return Sprite(frame, x0, y0), curve, head

    def _chapter(self, index):
        """Per-chapter kinetic text: a chip, headline words, and body lines."""
        s = self.s
        story = self.story
        beat = story.beats[index]
        count = len(story.beats)
        sprites = {"chip": None, "words": [], "body": []}

        chip_text = f"{index + 1:02d} / {count:02d}    {beat.label.upper()}"
        face = font(19 * s, "bold")
        text, pad = _text_image(chip_text, face, _lighten(self.accent, 0.15), 2.6 * s)
        chip_h = round(46 * s)
        chip_w = text.width - 2 * pad + round(40 * s)
        chip = Image.new("RGBA", (chip_w + 2, chip_h + 2), (0, 0, 0, 0))
        ImageDraw.Draw(chip).rounded_rectangle(
            (0, 0, chip_w, chip_h),
            radius=chip_h // 2,
            fill=self.accent + (40,),
            outline=self.accent + (130,),
            width=max(1, round(1.5 * s)),
        )
        ascent, descent = face.getmetrics()
        chip.alpha_composite(text, (round(20 * s - pad), round((chip_h - ascent - descent) / 2 - pad + 2 * s)))
        sprites["chip"] = Sprite(chip, MARGIN * s, CHIP_Y * s)

        headline = story.title if index == 0 else beat.title
        column = self.width - 2 * MARGIN * s
        lines, face, size = fit_lines(headline, 84 * s, "black", column, 2, minimum=40 * s)
        line_height = size * 1.08
        y = HEADLINE_Y * s
        words = [(line_index, word) for line_index, line in enumerate(lines) for word in line.split()]
        space = face.getlength(" ")
        cursor = {}
        for order, (line_index, word) in enumerate(words):
            last = order == len(words) - 1
            if last:
                image, pad = _gradient_text(word, face, _lighten(self.accent, 0.35), self.accent)
            else:
                image, pad = _gradient_text(word, face, INK, (214, 223, 238))
            x = cursor.get(line_index, MARGIN * s)
            sprites["words"].append(
                (Sprite(image, x - pad, y + line_index * line_height - pad), order)
            )
            cursor[line_index] = x + face.getlength(word) + space
        body_top = y + len(lines) * line_height + 30 * s
        body = story.hook if index == 0 else story.takeaway if index == count - 1 else beat.body
        space_left = SCENE_BOX[1] * s - 26 * s - body_top
        body_lines, face, size = fit_lines(body, 38 * s, "regular", column, 3, minimum=22 * s)
        while len(body_lines) * size * 1.32 > space_left and size > 22 * s:
            body_lines, face, size = fit_lines(body, size - 1, "regular", column, 3, minimum=22 * s)
        for line_index, line in enumerate(body_lines):
            image, pad = _text_image(line, face, SOFT)
            sprites["body"].append(
                (Sprite(image, MARGIN * s - pad, body_top + line_index * size * 1.32 - pad), line_index)
            )
        return sprites

    # ------------------------------------------------------------------ per frame

    def _scene(self, frame, scene):
        s = self.s
        x0, y0, x1, y1 = self.scene_box
        pixels = np.asarray(scene.convert("RGBA").resize(self.scene_size), dtype=np.float32) / 255
        alpha = pixels[..., 3:4]
        # VTK renders over transparent black, so its colors are already premultiplied.
        rgb = np.minimum(pixels[..., :3], alpha)
        pad = round(140 * s)
        top, bottom = max(0, y0 - pad), min(self.height, y1 + pad)
        canvas = np.zeros((bottom - top, x1 - x0, 3), np.float32)
        canvas[y0 - top : y1 - top] = rgb
        source = Image.fromarray((canvas * 255).astype(np.uint8))
        halos = []
        for factor, radius in ((8, 9), (4, 5)):
            small = source.reduce(factor).filter(ImageFilter.GaussianBlur(radius * s))
            halos.append(np.asarray(small.resize(source.size, Image.BILINEAR), np.float32) / 255)
        region = frame[top:bottom, x0:x1]
        region += halos[0] * 0.55 + halos[1] * 0.3
        inner = region[y0 - top : y1 - top]
        inner *= 1 - alpha
        inner += rgb
        inner += halos[1][y0 - top : y1 - top] * 0.16

    def _progress(self, frame, position):
        s = self.s
        count = len(self.story.beats)
        gap = 10 * s
        total = self.width - 2 * MARGIN * s
        segment = (total - gap * (count - 1)) / count
        height = max(2, round(5 * s))
        band = Image.new("RGBA", (self.width, height + 2), (0, 0, 0, 0))
        draw = ImageDraw.Draw(band)
        for index in range(count):
            left = MARGIN * s + index * (segment + gap)
            draw.rounded_rectangle((left, 0, left + segment, height), radius=height / 2, fill=(255, 255, 255, 38))
            amount = float(np.clip(position - index, 0, 1))
            if amount > 0:
                draw.rounded_rectangle(
                    (left, 0, left + max(height, amount * segment), height),
                    radius=height / 2,
                    fill=self.accent + (255,),
                )
        Sprite(band, 0, PROGRESS_Y * s).draw(frame)

    def _status(self, frame, phase):
        s = self.s
        if self.model.steps is not None:
            step = float(np.interp(phase, np.linspace(0, 1, len(self.model.steps)), self.model.steps))
            text = f"STEP {step:,.0f} / {self.model.steps[-1]:,}"
            if "dt" in self.model.params:
                seconds = step * self.model.params["dt"]
                digits = 2 if self.model.steps[-1] * self.model.params["dt"] < 5 else 1
                text = f"t = {seconds:.{digits}f} s"
        elif self.model.kind == "generic":
            text = "SCENE STUDY"
        else:
            text = "STEADY FIELD"
        face = font(17 * s, "bold")
        image, pad = _text_image(text, face, INK, 2.2 * s)
        height = round(40 * s)
        width = image.width - 2 * pad + round(46 * s)
        pill = Image.new("RGBA", (width + 2, height + 2), (0, 0, 0, 0))
        draw = ImageDraw.Draw(pill)
        draw.rounded_rectangle((0, 0, width, height), radius=height // 2, fill=(8, 12, 22, 170), outline=(255, 255, 255, 40))
        dot = 9 * s
        live = self.model.steps is not None
        draw.ellipse((16 * s, height / 2 - dot / 2, 16 * s + dot, height / 2 + dot / 2), fill=(255, 92, 92) if live else self.accent)
        ascent, descent = face.getmetrics()
        pill.alpha_composite(image, (round(32 * s - pad), round((height - ascent - descent) / 2 - pad + 1 * s)))
        Sprite(pill, self.width - MARGIN * s - width, self.scene_box[3] - height - 6 * s).draw(frame)

    def _card(self, frame, time, detail):
        s = self.s
        self.card.draw(frame)
        x0, y0, x1, y1 = self.card_box
        if self.chart_frame is not None:
            self.chart_frame.draw(frame)
            reveal = _smootherstep((time - 0.5) / 4.5)
            width = self.chart_curve.width
            cut = max(0, min(width, round(reveal * width)))
            if cut:
                Sprite(self.chart_curve.crop((0, 0, cut, self.chart_curve.height)), x0, y0).draw(frame)
            if 0 < reveal < 1:
                hx = reveal * width
                hy = float(np.interp(hx, self.chart_head[:, 0], self.chart_head[:, 1]))
                radius = 9 * s
                dot = Image.new("RGBA", (round(radius * 6), round(radius * 6)), (0, 0, 0, 0))
                ImageDraw.Draw(dot).ellipse((radius * 2, radius * 2, radius * 4, radius * 4), fill=(255, 255, 255, 255))
                halo = dot.filter(ImageFilter.GaussianBlur(radius))
                Sprite(halo, x0 + hx - radius * 3, y0 + hy - radius * 3).draw(frame, 1.0)
                Sprite(dot, x0 + hx - radius * 3, y0 + hy - radius * 3).draw(frame, 0.9)
        elif detail is not None:
            inset = round(10 * s)
            label = round(44 * s)
            size = (x1 - x0 - 2 * inset, y1 - y0 - 2 * inset - label)
            crop = ImageOps.fit(detail.convert("RGBA"), size, Image.LANCZOS)
            mask = Image.new("L", size, 0)
            ImageDraw.Draw(mask).rounded_rectangle((0, 0, size[0] - 1, size[1] - 1), radius=round(18 * s), fill=255)
            alpha = np.asarray(crop, np.float32)[..., 3] / 255 * np.asarray(mask, np.float32) / 255
            rgb = np.minimum(np.asarray(crop, np.float32)[..., :3] / 255, alpha[..., None])
            region = frame[y0 + inset + label : y0 + inset + label + size[1], x0 + inset : x0 + inset + size[0]]
            region *= 1 - alpha[..., None]
            region += rgb
            self.detail_label.draw(frame)
        elif self.note_card is not None:
            self.note_card.draw(frame)

    def _text(self, frame, index, local, length):
        s = self.s
        sprites = self.chapters[index]
        exit_fade = 1.0 if index == len(self.chapters) - 1 else _ease_out((length - local) / 0.35)
        rise = 36 * s

        def animate(sprite, start, travel=rise):
            amount = _ease_out((local - start) / 0.65)
            sprite.draw(frame, amount * exit_fade, dy=travel * (1 - amount))

        animate(sprites["chip"], 0.05, 18 * s)
        for sprite, order in sprites["words"]:
            animate(sprite, 0.18 + 0.075 * order)
        words = len(sprites["words"])
        for sprite, order in sprites["body"]:
            animate(sprite, 0.35 + 0.075 * words + 0.12 * order, 20 * s)

    def compose(self, scene, time, phase, detail=None):
        count = len(self.story.beats)
        length = self.duration / count
        position = min(time / self.duration * count, count - 1e-8)
        index = int(position)
        local = time - index * length
        frame = self.background.copy()
        self._scene(frame, scene)
        for sprite in self.static:
            sprite.draw(frame)
        self._card(frame, time, detail)
        self._status(frame, phase)
        self._progress(frame, position)
        self._text(frame, index, local, length)
        # Soft entrance/exit, with a complete visible ending before fading.
        fade = min(1.0, time / 0.45, (self.duration - time) / 0.65)
        if fade < 1:
            frame *= max(0.0, fade)
        return Image.fromarray((np.clip(frame, 0, 1) * 255 + 0.5).astype(np.uint8))


def contact_sheet(images):
    """Six portrait frames in reading order, for quickly reviewing the story."""
    tile_width = 300
    tile_height = round(images[0].height / images[0].width * tile_width)
    sheet = Image.new("RGB", (tile_width * 3, tile_height * 2), (8, 13, 23))
    for index, image in enumerate(images):
        sheet.paste(
            image.resize((tile_width, tile_height), Image.LANCZOS),
            ((index % 3) * tile_width, (index // 3) * tile_height),
        )
    return sheet
