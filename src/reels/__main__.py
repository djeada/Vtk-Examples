"""One example, one educational reel: ``python -m src.reels``."""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from src.__main__ import discover_examples, select_example
from src.reels.stories import Beat, Story, STORIES

QUALITIES = {"studio": (1080, 1920, 30), "fast": (720, 1280, 24)}


def load_story(key, path=None):
    title = key.replace("_", " ").title()
    base = STORIES.get(key) or Story(
        key,
        title,
        "Explore the geometry behind this visualization.",
        "",
        (70, 220, 235),
        tuple(
            Beat(label, heading, body)
            for label, heading, body in (
                (
                    "Explore",
                    title,
                    "Start with the scene produced by the original example.",
                ),
                (
                    "Geometry",
                    "Look at the shapes",
                    "Observe the geometry and how the objects occupy space.",
                ),
                (
                    "View",
                    "Change your perspective",
                    "Camera motion reveals depth and the relationships between objects.",
                ),
                (
                    "Read",
                    "Follow the visual structure",
                    "Look for boundaries, connections, and repeated patterns.",
                ),
                (
                    "Inspect",
                    "Study the details",
                    "This is the example's original VTK scene. Camera motion is presentation.",
                ),
                (
                    "Try it",
                    "Explore the example yourself",
                    "Run the interactive example to investigate its geometry further.",
                ),
            )
        ),
        "Open the original example to explore its scene interactively.",
        "Geometry, lighting, and viewpoint turn the dataset into a visible scene.",
    )
    if path is None:
        return base
    data = json.loads(Path(path).read_text())
    allowed = {"title", "hook", "takeaway", "model_note", "units", "accent", "beats"}
    if not isinstance(data, dict) or set(data) - allowed:
        raise ValueError(
            "Story must be an object containing only documented story fields"
        )
    limits = {"title": 60, "hook": 160, "takeaway": 180, "model_note": 220, "units": 60}
    for name, limit in limits.items():
        if name in data and (
            not isinstance(data[name], str)
            or not data[name].strip()
            or len(data[name]) > limit
        ):
            raise ValueError(f"Story {name} must contain 1–{limit} characters")
    if "accent" in data:
        accent = data["accent"]
        if (
            not isinstance(accent, list)
            or len(accent) != 3
            or any(type(n) is not int or not 0 <= n <= 255 for n in accent)
        ):
            raise ValueError("Story accent must be three RGB integers from 0 to 255")
        data["accent"] = tuple(accent)
    if "beats" in data:
        beats = data["beats"]
        if not isinstance(beats, list) or len(beats) != 6:
            raise ValueError("Story must contain exactly six beats")
        for beat in beats:
            if not isinstance(beat, dict) or set(beat) != {"label", "title", "body"}:
                raise ValueError("Each beat needs label, title, and body")
            for name, limit in (("label", 40), ("title", 60), ("body", 180)):
                if (
                    not isinstance(beat[name], str)
                    or not beat[name].strip()
                    or len(beat[name]) > limit
                ):
                    raise ValueError(f"Beat {name} must contain 1–{limit} characters")
        data["beats"] = tuple(Beat(**beat) for beat in beats)
    return replace(base, **data)


def subtitle_text(story, duration):
    def stamp(seconds):
        milliseconds = round(seconds * 1000)
        hours, remainder = divmod(milliseconds, 3600000)
        minutes, remainder = divmod(remainder, 60000)
        seconds, milliseconds = divmod(remainder, 1000)
        return f"{hours:02}:{minutes:02}:{seconds:02},{milliseconds:03}"

    captions = []
    for index, beat in enumerate(story.beats):
        heading = story.title if index == 0 else beat.title
        body = story.hook if index == 0 else story.takeaway if index == 5 else beat.body
        captions.append(
            f"{index + 1}\n{stamp(index * duration / 6)} --> {stamp((index + 1) * duration / 6)}\n{heading}\n{body}\n"
        )
    return "\n".join(captions)


def render_reel(path, output, args):
    import numpy as np
    from src.reels.encoding import VideoEncoder, probe_video, synth_soundtrack
    from src.reels.layout import ReelLayout, contact_sheet
    from src.reels.models import prepare_model, source_fingerprint
    from src.reels.scene import SceneRenderer, capture_example

    started = time.monotonic()
    width, height, fps = QUALITIES[args.quality]
    frame_count = round(args.duration * fps)
    duration = frame_count / fps
    story = load_story(path.stem, args.story)
    output.parent.mkdir(parents=True, exist_ok=True)
    sidecars = [
        output.with_suffix(suffix)
        for suffix in (".jpg", ".preview.jpg", ".json", ".srt")
    ]
    protected = sidecars + ([] if args.preview_only else [output])
    if not args.overwrite:
        for destination in protected:
            if destination.exists():
                raise FileExistsError(
                    f"Output already exists: {destination}. Use --overwrite to replace it."
                )
    print(f"\n{path.stem}: solving / capturing ({args.quality})", flush=True)
    last_progress = [-1]

    def progress(fraction):
        percent = int(fraction * 10) * 10
        if percent > last_progress[0]:
            print(f"  solver {percent}%", flush=True)
            last_progress[0] = percent

    scene_width, scene_height = round(624 * width / 720), round(435 * width / 720)
    generic_renderer = None
    if path.stem in STORIES:
        model = prepare_model(
            path.stem,
            args.quality,
            None if args.no_cache else Path(".cache/reels"),
            progress,
        )
        fingerprint = source_fingerprint(path.stem)
    else:
        model, generic_renderer = capture_example(path, scene_width, scene_height)
        fingerprint = hashlib.sha256(path.read_bytes()).hexdigest()
    solving_seconds = time.monotonic() - started
    scene = SceneRenderer(model, scene_width, scene_height, generic_renderer)
    try:
        layout = ReelLayout(
            story,
            model,
            width,
            height,
            duration,
            scene.colors,
            None if generic_renderer is not None else scene.legend,
        )
        with tempfile.TemporaryDirectory(prefix=".reel-", dir=output.parent) as folder:
            temporary = Path(folder)
            previews = []
            for index in range(6):
                phase = (index + 0.5) / 6
                image = scene.render(phase)
                previews.append(
                    layout.compose(
                        image,
                        phase * duration,
                        phase,
                        layout.detail(image, scene.renderer),
                    )
                )
            contact_sheet(previews).save(temporary / "preview.jpg", quality=95)
            previews[0].save(temporary / "thumbnail.jpg", quality=95)
            (temporary / "captions.srt").write_text(subtitle_text(story, duration))
            audio = args.audio
            if not args.silent and audio is None:
                audio = temporary / "soundtrack.wav"
                synth_soundtrack(audio, duration, seed=int(fingerprint[:8], 16))
            if not args.preview_only:
                print(
                    f"  encoding {frame_count:,} frames; {duration:g}s video",
                    flush=True,
                )
                with VideoEncoder(
                    output,
                    width,
                    height,
                    fps,
                    frame_count,
                    audio=None if args.silent else audio,
                    overwrite=args.overwrite,
                ) as encoder:
                    for index in range(frame_count):
                        phase = index / max(1, frame_count - 1)
                        image = scene.render(phase)
                        frame = layout.compose(
                            image,
                            index / fps,
                            phase,
                            layout.detail(image, scene.renderer),
                        )
                        encoder.write(np.asarray(frame))
                        if (index + 1) % max(1, frame_count // 10) == 0:
                            print(
                                f"  video {(index + 1) / frame_count:.0%} · elapsed {time.monotonic() - started:.1f}s",
                                flush=True,
                            )
            manifest = {
                "example": str(path.relative_to(Path(__file__).resolve().parents[1])),
                "source_sha256": fingerprint,
                "story": asdict(story),
                "mode": "preview" if args.preview_only else "video",
                "video": {
                    "width": width,
                    "height": height,
                    "fps": fps,
                    "duration_seconds": duration,
                    "frames": frame_count,
                    "quality": args.quality,
                },
                "model": {
                    "parameters": model.params,
                    "note": model.note,
                    "field": model.field_name,
                    "units": model.unit,
                    "range": list(scene.legend[0])
                    if generic_renderer is None
                    else None,
                    "metrics": model.metrics,
                    "cache_hit": model.cache_hit,
                },
                "audio": "silent"
                if args.silent
                else str(args.audio)
                if args.audio
                else "original synthesized ambient soundtrack",
                "generation_seconds": round(time.monotonic() - started, 3),
                "solving_seconds": round(solving_seconds, 3),
                "encoding": None if args.preview_only else probe_video(output),
            }
            (temporary / "manifest.json").write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)
                + "\n"
            )
            for source, destination in zip(
                ("thumbnail.jpg", "preview.jpg", "manifest.json", "captions.srt"),
                sidecars,
            ):
                (temporary / source).replace(destination)
        print(
            f"  saved {sidecars[1] if args.preview_only else output} ({time.monotonic() - started:.1f}s generation)",
            flush=True,
        )
    finally:
        scene.close()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Create one educational vertical reel per example."
    )
    parser.add_argument(
        "example", nargs="?", help="example name; defaults to obstacle_flow_simulation"
    )
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--all", choices=["cfd"], help="create nine separate CFD reels")
    parser.add_argument("--duration", type=float, default=60)
    parser.add_argument("--quality", choices=QUALITIES, default="studio")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--story", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--preview-only", action="store_true")
    audio = parser.add_mutually_exclusive_group()
    audio.add_argument("--silent", action="store_true")
    audio.add_argument("--audio", type=Path)
    args = parser.parse_args(argv)
    if args.list:
        for key, story in STORIES.items():
            if not args.example or args.example.lower() in key:
                print(f"{key:42} {story.title}")
        print(
            "Other native VTK examples accept their name and an optional --story JSON."
        )
        return 0
    if not math.isfinite(args.duration) or args.duration < 1:
        parser.error("--duration must be a finite number of seconds, at least 1")
    if args.all and (args.example or args.story):
        parser.error(
            "--all cannot be combined with an example or --story; each story belongs to one example"
        )
    if args.audio and not args.audio.is_file():
        parser.error(f"Audio file does not exist: {args.audio}")
    if not args.preview_only and any(
        shutil.which(binary) is None for binary in ("ffmpeg", "ffprobe")
    ):
        parser.error("Install FFmpeg (including ffprobe) and put it on PATH")
    if (
        sys.platform.startswith("linux")
        and not os.environ.get("DISPLAY")
        and not os.environ.get("VTK_REELS_HEADLESS")
    ):
        if not shutil.which("xvfb-run"):
            parser.error(
                "Headless Linux needs Xvfb: install xvfb, then run this command again"
            )
        environment = dict(os.environ, VTK_REELS_HEADLESS="1")
        try:
            return subprocess.call(
                [
                    "xvfb-run",
                    "-a",
                    sys.executable,
                    "-m",
                    "src.reels",
                    *(sys.argv[1:] if argv is None else argv),
                ],
                env=environment,
            )
        except KeyboardInterrupt:
            return 130
    try:
        examples = discover_examples()
        paths = (
            [select_example(key, examples) for key in STORIES]
            if args.all
            else [select_example(args.example or "obstacle_flow_simulation", examples)]
        )
        if args.all and args.output and args.output.suffix.lower() == ".mp4":
            raise ValueError(
                "Batch --output must be a directory; every example gets its own MP4"
            )
        if not args.all and args.output and args.output.suffix.lower() != ".mp4":
            raise ValueError("Single-example --output must be an .mp4 path")
        for path in paths:
            output = (
                (args.output or Path("reels")) / f"{path.stem}.mp4"
                if args.all
                else args.output or Path("reels") / f"{path.stem}.mp4"
            )
            render_reel(path, output, args)
        return 0
    except KeyboardInterrupt:
        print("Rendering cancelled.", file=sys.stderr)
        return 130
    except ImportError as error:
        print(
            f"Missing dependency: {error}. Install requirements-reels.txt (and optional dependencies for the selected example).",
            file=sys.stderr,
        )
        return 2
    except (ValueError, OSError, RuntimeError) as error:
        print(f"Reel export failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
