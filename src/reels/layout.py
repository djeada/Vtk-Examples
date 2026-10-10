"""Readable vertical film layout: cached typography, captions, and scientific plots."""

from functools import lru_cache
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageOps


INK = (239, 244, 250)
MUTED = (151, 168, 190)
PANEL = (12, 22, 36)


@lru_cache(maxsize=64)
def font(size, bold=False):
    # Prefer Inter when installed; Matplotlib bundles a portable DejaVu fallback.
    name = "InterDisplay-Bold.otf" if bold else "Inter-Regular.otf"
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


def text_block(draw, text, xy, width, size, fill=INK, bold=False, max_lines=None):
    face = font(size, bold)
    lines = wrap(text, face, width)
    if max_lines is not None:
        while len(lines) > max_lines and size > 10:
            size -= 1
            face = font(size, bold)
            lines = wrap(text, face, width)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=face, fill=fill, stroke_width=0)
        y += size * 1.22
    return y


class ReelLayout:
    """Six separate story beats share a consistent, mobile-readable visual system."""

    def __init__(self, story, model, width, height, duration, colors, legend):
        self.story, self.model = story, model
        self.width, self.height, self.duration = width, height, duration
        self.scale = width / 720
        self.colors, self.legend = colors, legend
        self.scene_box = tuple(
            round(value * self.scale) for value in (48, 325, 672, 760)
        )
        self.scene_size = (
            self.scene_box[2] - self.scene_box[0],
            self.scene_box[3] - self.scene_box[1],
        )
        self.chart_box = tuple(
            round(value * self.scale) for value in (48, 846, 672, 997)
        )
        self.cards = [self._card(index) for index in range(len(story.beats))]
        self.chart_image = self._chart() if model.chart else None

    def _background(self):
        yy, xx = np.mgrid[0 : self.height, 0 : self.width]
        glow = np.exp(
            -(
                (xx / self.width - 0.75) ** 2 / 0.2
                + (yy / self.height - 0.43) ** 2 / 0.13
            )
        )
        noise = (np.sin(xx * 0.017) * np.sin(yy * 0.013) + 1) * 0.5
        pixels = np.stack(
            [8 + 8 * glow + noise, 13 + 12 * glow + noise, 23 + 17 * glow + noise],
            axis=-1,
        ).astype(np.uint8)
        image = Image.fromarray(pixels)
        draw = ImageDraw.Draw(image)
        spacing = round(48 * self.scale)
        for px in range(0, self.width, spacing):
            draw.line((px, 0, px, self.height), fill=(15, 24, 37), width=1)
        return image

    def _card(self, index):
        s = self.scale
        story = self.story
        beat = story.beats[index]
        image = self._background()
        draw = ImageDraw.Draw(image)

        def line(text, x, y, size, color=INK, bold=False):
            draw.text(
                (round(x * s), round(y * s)),
                text,
                font=font(size * s, bold),
                fill=color,
            )

        accent = story.accent
        draw.rounded_rectangle(
            tuple(round(v * s) for v in (48, 48, 74, 74)),
            radius=round(7 * s),
            fill=accent,
        )
        line("VTK / FIELD NOTES", 86, 51, 15, MUTED, True)
        line(
            "SCIENCE IN 60s" if self.duration == 60 else "SCIENCE IN MOTION",
            483,
            51,
            13,
            MUTED,
        )
        line(f"{index + 1:02d} / {len(story.beats):02d}", 48, 109, 18, accent, True)
        line(beat.label.upper(), 142, 110, 15, accent, True)
        headline = story.title if index == 0 else beat.title
        text_block(
            draw, headline, (48 * s, 151 * s), 624 * s, 55 * s, bold=True, max_lines=2
        )
        draw.rounded_rectangle(
            tuple(round(v * s) for v in (46, 323, 674, 762)),
            radius=round(20 * s),
            fill=PANEL,
            outline=(34, 50, 70),
            width=max(1, round(s)),
        )
        line(self.model.field_name, 66, 343, 12, MUTED, True)
        if self.legend is not None:
            (low, high), unit = self.legend
            bar = Image.fromarray(
                np.broadcast_to(self.colors[None], (10, 256, 3)).copy()
            ).resize((round(256 * s), round(9 * s)))
            image.paste(bar, (round(48 * s), round(789 * s)))
            line(f"{low:.3g}", 48, 809, 16, MUTED)
            upper = f"{high:.3g} {unit}"
            line(upper, 304 - font(16 * s).getlength(upper) / s, 809, 16, MUTED)
        else:
            line("ORIGINAL VTK SCENE", 48, 796, 15, MUTED, True)
        for metric, (label, value) in enumerate(self.model.metrics[:2]):
            mx = 343 + metric * 173
            line(label.upper(), mx, 783, 10, MUTED, True)
            text_block(
                draw,
                value,
                (mx * s, 802 * s),
                164 * s,
                23 * s,
                fill=accent,
                bold=True,
                max_lines=1,
            )
        if not self.model.chart:
            draw.rounded_rectangle(
                tuple(round(v * s) for v in (48, 846, 672, 997)),
                radius=round(16 * s),
                fill=PANEL,
                outline=(34, 50, 70),
            )
        body = (
            story.hook
            if index == 0
            else story.takeaway
            if index == len(story.beats) - 1
            else beat.body
        )
        text_block(draw, body, (48 * s, 1030 * s), 624 * s, 27 * s, max_lines=3)
        draw.line(
            tuple(round(v * s) for v in (48, 1140, 672, 1140)),
            fill=(38, 53, 72),
            width=max(1, round(s)),
        )
        text_block(
            draw,
            self.model.note,
            (48 * s, 1157 * s),
            624 * s,
            15 * s,
            fill=MUTED,
            max_lines=2,
        )
        line(f"{self.model.key.replace('_', ' ')}", 48, 1233, 11, MUTED)
        return image

    def _chart(self):
        s = self.scale
        width = self.chart_box[2] - self.chart_box[0]
        height = self.chart_box[3] - self.chart_box[1]
        image = Image.new("RGB", (width, height), PANEL)
        draw = ImageDraw.Draw(image)
        chart = self.model.chart
        x = np.asarray(chart["x"])
        series = [np.asarray(values) for values in chart["series"]]
        values = np.concatenate([values[np.isfinite(values)] for values in series])
        low, high = float(values.min()), float(values.max())
        if high == low:
            high = low + 1
        left, right = 55 * s, width - 20 * s
        top, bottom = 34 * s, height - 33 * s
        draw.text(
            (12 * s, 8 * s), "COMPUTED PROFILE", font=font(10 * s, True), fill=MUTED
        )
        draw.text((12 * s, top - 7 * s), f"{high:.3g}", font=font(11 * s), fill=MUTED)
        draw.text((12 * s, bottom - 7 * s), f"{low:.3g}", font=font(11 * s), fill=MUTED)
        draw.line(
            (left, bottom, right, bottom), fill=(69, 87, 109), width=max(1, round(s))
        )
        draw.line((left, top, left, bottom), fill=(69, 87, 109), width=max(1, round(s)))
        draw.line(
            (left, (top + bottom) / 2, right, (top + bottom) / 2),
            fill=(29, 46, 67),
            width=1,
        )
        palette = (self.story.accent, (197, 213, 238), (140, 192, 255))
        for index, ys in enumerate(series):
            points = list(
                zip(
                    left + (x - x.min()) / np.ptp(x) * (right - left),
                    bottom - (ys - low) / (high - low) * (bottom - top),
                )
            )
            # Reference is dashed so coincident numerical/analytic curves are distinguishable.
            if index == 1:
                for start in range(0, len(points) - 1, 5):
                    draw.line(
                        points[start : start + 3],
                        fill=palette[index],
                        width=max(1, round(2 * s)),
                    )
            else:
                draw.line(points, fill=palette[index % 3], width=max(1, round(2 * s)))
            draw.text(
                ((330 + index * 125) * s, 8 * s),
                chart["labels"][index],
                font=font(11 * s),
                fill=palette[index % 3],
            )
        draw.text(
            (left, bottom + 9 * s), f"{x.min():.3g}", font=font(11 * s), fill=MUTED
        )
        draw.text(
            (right - 20 * s, bottom + 9 * s),
            f"{x.max():.3g}",
            font=font(11 * s),
            fill=MUTED,
        )
        label = chart["xlabel"] + " / " + chart["ylabel"]
        draw.text(
            ((left + right) / 2 - font(11 * s).getlength(label) / 2, bottom + 9 * s),
            label,
            font=font(11 * s),
            fill=MUTED,
        )
        return image

    def compose(self, scene, time, phase, detail=None):
        s = self.scale
        count = len(self.story.beats)
        position = min(time / self.duration * count, count - 1e-8)
        index = int(position)
        image = self.cards[index].copy()
        # Short typography crossfades provide transitions without inventing data.
        local = position - index
        if index > 0 and local < 0.035:
            alpha = local / 0.035
            image = Image.blend(
                self.cards[index - 1], image, alpha * alpha * (3 - 2 * alpha)
            )
        scene = scene.resize(self.scene_size)
        image.paste(scene, (self.scene_box[0], self.scene_box[1]))
        draw = ImageDraw.Draw(image)
        draw.text(
            (66 * s, 343 * s),
            self.model.field_name,
            font=font(12 * s, True),
            fill=MUTED,
        )
        if self.chart_image is not None:
            image.paste(self.chart_image, self.chart_box[:2])
        elif detail is not None:
            detail = ImageOps.fit(
                detail,
                (
                    self.chart_box[2] - self.chart_box[0],
                    self.chart_box[3] - self.chart_box[1],
                ),
            )
            image.paste(detail, self.chart_box[:2])
            draw = ImageDraw.Draw(image)
            draw.rectangle(
                (
                    self.chart_box[0],
                    self.chart_box[1],
                    self.chart_box[2],
                    self.chart_box[1] + 24 * s,
                ),
                fill=PANEL,
            )
            draw.text(
                (self.chart_box[0] + 12 * s, self.chart_box[1] + 6 * s),
                "DETAIL / SAME COLOR SCALE",
                font=font(10 * s, True),
                fill=MUTED,
            )
        else:
            draw = ImageDraw.Draw(image)
            text_block(
                draw,
                self.story.model_note,
                (65 * s, 876 * s),
                580 * s,
                18 * s,
                fill=MUTED,
                max_lines=3,
            )
        draw = ImageDraw.Draw(image)
        if self.model.steps is not None:
            step = float(
                np.interp(
                    phase, np.linspace(0, 1, len(self.model.steps)), self.model.steps
                )
            )
            draw.text(
                (48 * s, 1192 * s),
                f"RECORDED STEP {step:,.0f} / {self.model.steps[-1]:,}",
                font=font(13 * s, True),
                fill=self.story.accent,
            )
        else:
            draw.text(
                (48 * s, 1192 * s),
                "STEADY FIELD" if self.model.kind != "generic" else "SCENE STUDY",
                font=font(13 * s, True),
                fill=self.story.accent,
            )
        for segment in range(count):
            left = (48 + segment * 105) * s
            right = left + 95 * s
            draw.rounded_rectangle(
                (left, 1217 * s, right, 1221 * s), radius=2 * s, fill=(43, 58, 78)
            )
            amount = float(np.clip(position - segment, 0, 1))
            if amount > 0:
                draw.rounded_rectangle(
                    (left, 1217 * s, left + amount * (right - left), 1221 * s),
                    radius=2 * s,
                    fill=self.story.accent,
                )
        # Soft entrance/exit, with a complete visible ending before fading.
        fade = min(1, time / 0.45, (self.duration - time) / 0.65)
        if fade < 1:
            image = Image.blend(
                Image.new("RGB", image.size, (8, 13, 23)), image, max(0, fade)
            )
        return image

    def detail(self, scene, renderer):
        """A real crop around the obstacle, located with VTK's world/display map."""
        if self.model.solid is None:
            return None
        rows, cols = np.where(self.model.solid)
        halfx = (cols.max() - cols.min() + 1) * 3.2
        halfy = (rows.max() - rows.min() + 1) * 1.65
        cx, cy = (cols.max() + cols.min()) / 2, (rows.max() + rows.min()) / 2
        pixels = []
        for x, y in ((cx - halfx, cy - halfy), (cx + halfx, cy + halfy)):
            renderer.SetWorldPoint(float(x), float(y), 0, 1)
            renderer.WorldToDisplay()
            px, py, _ = renderer.GetDisplayPoint()
            pixels.append((px, scene.height - py))
        left = max(0, min(p[0] for p in pixels))
        right = min(scene.width, max(p[0] for p in pixels))
        top = max(0, min(p[1] for p in pixels))
        bottom = min(scene.height, max(p[1] for p in pixels))
        return (
            scene.crop((round(left), round(top), round(right), round(bottom)))
            if right > left and bottom > top
            else None
        )


def contact_sheet(images):
    """Six portrait frames in reading order, for quickly reviewing the story."""
    tile_width = 300
    tile_height = round(images[0].height / images[0].width * tile_width)
    sheet = Image.new("RGB", (tile_width * 3, tile_height * 2), (8, 13, 23))
    for index, image in enumerate(images):
        sheet.paste(
            image.resize((tile_width, tile_height)),
            ((index % 3) * tile_width, (index // 3) * tile_height),
        )
    return sheet
