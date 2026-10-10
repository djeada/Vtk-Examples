"""Streaming FFmpeg encoding helpers for short educational videos."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
import wave

import numpy as np


_VALID_PRESETS = {
    "ultrafast",
    "superfast",
    "veryfast",
    "faster",
    "fast",
    "medium",
    "slow",
    "slower",
    "veryslow",
}
_STDERR_LIMIT = 64 * 1024
_AUDIO_SAMPLE_RATE = 22_050


def _positive_integer(name: str, value: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be a positive integer")
    if value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _stderr_text(data: bytearray) -> str:
    return data.decode("utf-8", errors="replace").strip()


def _write_all(stream, payload: memoryview) -> None:
    """Write an entire frame payload even when an unbuffered pipe is partial."""
    remaining = payload
    while remaining:
        written = stream.write(remaining)
        if written is None or written <= 0:
            raise BrokenPipeError("FFmpeg pipe accepted no frame data")
        remaining = remaining[written:]


class VideoEncoder:
    """Write RGB frames to an atomically published MP4 through FFmpeg.

    Use as a context manager and call :meth:`write` once per frame. If rendering
    raises, the FFmpeg process is stopped and its temporary output is removed.
    The requested output is only published after all frames encode successfully.
    """

    def __init__(
        self,
        output: Path,
        width: int,
        height: int,
        fps: int,
        frame_count: int,
        audio: Path | None = None,
        preset: str = "fast",
        overwrite: bool = False,
    ) -> None:
        self.output = Path(output)
        self.width = _positive_integer("width", width)
        self.height = _positive_integer("height", height)
        self.fps = _positive_integer("fps", fps)
        self.frame_count = _positive_integer("frame_count", frame_count)
        if self.width % 2 or self.height % 2:
            raise ValueError("width and height must be even for yuv420p output")
        if self.output.suffix.lower() != ".mp4":
            raise ValueError("output filename must have an .mp4 extension")
        if not isinstance(preset, str) or preset not in _VALID_PRESETS:
            raise ValueError(
                f"preset must be one of: {', '.join(sorted(_VALID_PRESETS))}"
            )
        if not isinstance(overwrite, bool):
            raise ValueError("overwrite must be a boolean")
        self.preset = preset
        self.overwrite = overwrite
        self.audio = Path(audio) if audio is not None else None
        self._process: subprocess.Popen | None = None
        self._temp_path: Path | None = None
        self._stderr = bytearray()
        self._stderr_thread: threading.Thread | None = None
        self._written_frames = 0
        self._closed = False

        if not self.output.parent.exists():
            raise FileNotFoundError(
                f"Output directory does not exist: {self.output.parent}"
            )
        if not self.output.parent.is_dir():
            raise NotADirectoryError(
                f"Output parent is not a directory: {self.output.parent}"
            )
        if self.audio is not None and not self.audio.is_file():
            raise FileNotFoundError(f"Audio file does not exist: {self.audio}")
        if self.output.exists() and not self.overwrite:
            raise FileExistsError(f"Output already exists: {self.output}")

    @property
    def duration(self) -> float:
        """Requested video duration in seconds."""
        return self.frame_count / self.fps

    def _command(self, temp_path: Path) -> list[str]:
        ffmpeg = shutil.which("ffmpeg")
        if ffmpeg is None:
            raise FileNotFoundError(
                "FFmpeg executable was not found; install FFmpeg and make `ffmpeg` available on PATH."
            )
        command = [
            ffmpeg,
            "-hide_banner",
            "-loglevel",
            "error",
            "-nostdin",
            "-y",
            "-f",
            "rawvideo",
            "-vcodec",
            "rawvideo",
            "-pix_fmt",
            "rgb24",
            "-s:v",
            f"{self.width}x{self.height}",
            "-r",
            str(self.fps),
            "-i",
            "pipe:0",
        ]
        if self.audio is not None:
            command.extend(["-stream_loop", "-1", "-i", str(self.audio)])
        command.extend(
            [
                "-frames:v",
                str(self.frame_count),
                "-t",
                f"{self.duration:.12f}",
                "-c:v",
                "libx264",
                "-preset",
                self.preset,
                "-crf",
                "18",
                "-pix_fmt",
                "yuv420p",
                "-movflags",
                "+faststart",
            ]
        )
        if self.audio is not None:
            command.extend(["-c:a", "aac", "-b:a", "128k"])
        else:
            command.append("-an")
        command.append(str(temp_path))
        return command

    def __enter__(self) -> VideoEncoder:
        if self._process is not None or self._closed:
            raise RuntimeError("VideoEncoder context cannot be entered more than once")
        if self.output.exists() and not self.overwrite:
            raise FileExistsError(f"Output already exists: {self.output}")

        fd, temp_name = tempfile.mkstemp(
            prefix=f".{self.output.stem}.", suffix=".tmp.mp4", dir=self.output.parent
        )
        os.close(fd)
        self._temp_path = Path(temp_name)
        try:
            command = self._command(self._temp_path)
            self._process = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                bufsize=0,
            )
        except FileNotFoundError as error:
            self._cleanup_temp()
            if "FFmpeg executable" in str(error):
                raise
            raise FileNotFoundError(f"Could not start FFmpeg: {error}") from error
        except Exception:
            self._cleanup_temp()
            raise

        self._stderr_thread = threading.Thread(
            target=self._drain_stderr,
            args=(self._process.stderr,),
            name="ffmpeg-stderr-reader",
            daemon=True,
        )
        self._stderr_thread.start()
        return self

    def _drain_stderr(self, stream) -> None:
        while True:
            chunk = stream.read(4096)
            if not chunk:
                break
            self._stderr.extend(chunk)
            if len(self._stderr) > _STDERR_LIMIT:
                del self._stderr[:-_STDERR_LIMIT]

    def write(self, frame: np.ndarray) -> None:
        """Encode one uint8 RGB frame with shape ``(height, width, 3)``."""
        if self._process is None or self._closed:
            raise RuntimeError("write() requires an active VideoEncoder context")
        if self._written_frames >= self.frame_count:
            raise ValueError(
                f"received more than the expected {self.frame_count} frames"
            )
        if not isinstance(frame, np.ndarray):
            raise TypeError("frame must be a numpy.ndarray")
        if frame.dtype != np.uint8:
            raise ValueError("frame dtype must be uint8")
        if frame.shape != (self.height, self.width, 3):
            raise ValueError(
                f"frame shape must be {(self.height, self.width, 3)}, got {frame.shape}"
            )
        try:
            payload = memoryview(np.ascontiguousarray(frame)).cast("B")
            _write_all(self._process.stdin, payload)
        except (BrokenPipeError, OSError) as error:
            message = _stderr_text(self._stderr)
            self._abort()
            detail = f": {message}" if message else ""
            raise RuntimeError(
                f"FFmpeg stopped while accepting frames{detail}"
            ) from error
        self._written_frames += 1

    def _cleanup_temp(self) -> None:
        if self._temp_path is not None:
            try:
                self._temp_path.unlink(missing_ok=True)
            finally:
                self._temp_path = None

    def _abort(self) -> None:
        process = self._process
        if process is not None:
            if process.stdin is not None and not process.stdin.closed:
                try:
                    process.stdin.close()
                except OSError:
                    pass
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
        if self._stderr_thread is not None:
            self._stderr_thread.join(timeout=2)
        if process is not None and process.stderr is not None:
            process.stderr.close()
        self._process = None
        self._closed = True
        self._cleanup_temp()

    def _publish(self) -> None:
        assert self._temp_path is not None
        if self.overwrite:
            os.replace(self._temp_path, self.output)
        else:
            # link() provides atomic no-clobber publication if another process
            # created the destination after the initial existence check.
            os.link(self._temp_path, self.output)
            self._temp_path.unlink()
        self._temp_path = None

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        process = self._process
        if process is None:
            self._closed = True
            return False
        if exc_type is not None:
            self._abort()
            return False

        if self._written_frames != self.frame_count:
            written = self._written_frames
            self._abort()
            raise ValueError(
                f"expected {self.frame_count} frames but received {written}"
            )

        try:
            try:
                process.stdin.close()
            except BrokenPipeError:
                pass
            return_code = process.wait()
            if self._stderr_thread is not None:
                self._stderr_thread.join()
            if process.stderr is not None:
                process.stderr.close()
            self._process = None
            self._closed = True
            if return_code != 0:
                message = _stderr_text(self._stderr)
                detail = f": {message}" if message else ""
                raise RuntimeError(
                    f"FFmpeg encoding failed with exit code {return_code}{detail}"
                )
            if (
                self._temp_path is None
                or not self._temp_path.is_file()
                or self._temp_path.stat().st_size == 0
            ):
                raise RuntimeError("FFmpeg completed without producing a valid MP4")
            self._publish()
        except BaseException:
            if self._process is not None:
                self._abort()
            else:
                self._cleanup_temp()
            raise
        return False


def probe_video(path: Path) -> dict:
    """Return FFprobe's JSON metadata for a video or audio/video file."""
    video_path = Path(path)
    if not video_path.is_file():
        raise FileNotFoundError(f"Video file does not exist: {video_path}")
    ffprobe = shutil.which("ffprobe")
    if ffprobe is None:
        raise FileNotFoundError(
            "FFprobe executable was not found; install FFmpeg and make `ffprobe` available on PATH."
        )
    result = subprocess.run(
        [
            ffprobe,
            "-v",
            "error",
            "-show_format",
            "-show_streams",
            "-of",
            "json",
            str(video_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        detail = result.stderr.strip()
        raise RuntimeError(
            f"FFprobe failed with exit code {result.returncode}: {detail}"
        )
    try:
        metadata = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError("FFprobe returned invalid JSON metadata") from error
    if not isinstance(metadata, dict):
        raise RuntimeError("FFprobe metadata must be a JSON object")
    return metadata


def synth_soundtrack(path: Path, duration: float, seed: int = 0) -> Path:
    """Write a deterministic, soft ambient pulse bed as a stereo PCM WAV.

    The sound is built from quiet pitched oscillators with slow amplitude swells,
    gentle stereo variation, and fade-in/fade-out envelopes. No sampled material
    or noise source is used.
    """
    output = Path(path)
    try:
        duration_value = float(duration)
    except (TypeError, ValueError) as error:
        raise ValueError("duration must be finite and positive") from error
    if not math.isfinite(duration_value) or duration_value <= 0:
        raise ValueError("duration must be finite and positive")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    if output.suffix.lower() != ".wav":
        raise ValueError("soundtrack filename must have a .wav extension")
    if not output.parent.exists():
        raise FileNotFoundError(f"Output directory does not exist: {output.parent}")
    if not output.parent.is_dir():
        raise NotADirectoryError(f"Output parent is not a directory: {output.parent}")

    frames = max(1, round(duration_value * _AUDIO_SAMPLE_RATE))
    times = np.arange(frames, dtype=np.float64) / _AUDIO_SAMPLE_RATE
    rng = np.random.default_rng(int(seed))
    phase = rng.uniform(0, 2 * np.pi, size=4)
    pulse_phase = rng.uniform(0, 2 * np.pi)
    slow_phase = rng.uniform(0, 2 * np.pi)
    pulse = (
        0.82 + 0.18 * (0.5 + 0.5 * np.cos(2 * np.pi * times / 4.8 + pulse_phase)) ** 2
    )
    slow = 0.91 + 0.09 * np.sin(2 * np.pi * times / 17.0 + slow_phase)
    fade_seconds = min(2.0, duration_value / 4.0)
    envelope = np.ones(frames, dtype=np.float64)
    if fade_seconds:
        fade_frames = max(1, int(fade_seconds * _AUDIO_SAMPLE_RATE))
        fade_in = np.linspace(0.0, 1.0, fade_frames, endpoint=True)
        fade = 0.5 - 0.5 * np.cos(np.pi * fade_in)
        envelope[:fade_frames] *= fade[: min(fade_frames, frames)]
        envelope[-fade_frames:] *= fade[-min(fade_frames, frames) :]

    left = (
        (
            0.30 * np.sin(2 * np.pi * 55.0 * times + phase[0])
            + 0.19 * np.sin(2 * np.pi * 82.41 * times + phase[1])
            + 0.12 * np.sin(2 * np.pi * 110.0 * times + phase[2])
            + 0.045 * np.sin(2 * np.pi * 220.0 * times + phase[3])
        )
        * pulse
        * slow
    )
    right = (
        (
            0.30 * np.sin(2 * np.pi * 55.0 * times + phase[0] + 0.08)
            + 0.19 * np.sin(2 * np.pi * 82.41 * times + phase[1] - 0.06)
            + 0.12 * np.sin(2 * np.pi * 110.0 * times + phase[2] + 0.05)
            + 0.045 * np.sin(2 * np.pi * 220.0 * times + phase[3] - 0.04)
        )
        * pulse
        * slow
    )
    stereo = np.column_stack((left, right)) * envelope[:, None] * 0.28
    pcm = np.rint(np.clip(stereo, -1.0, 1.0) * 32767).astype("<i2")

    fd, temp_name = tempfile.mkstemp(
        prefix=f".{output.stem}.", suffix=".tmp.wav", dir=output.parent
    )
    os.close(fd)
    temp_path = Path(temp_name)
    try:
        with wave.open(str(temp_path), "wb") as wav_file:
            wav_file.setnchannels(2)
            wav_file.setsampwidth(2)
            wav_file.setframerate(_AUDIO_SAMPLE_RATE)
            wav_file.writeframes(pcm.tobytes())
        os.replace(temp_path, output)
    finally:
        temp_path.unlink(missing_ok=True)
    return output
