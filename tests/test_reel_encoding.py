"""Tests for the FFmpeg reel encoder and deterministic audio generator."""

from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import wave

import numpy as np

from src.reels.encoding import VideoEncoder, _write_all, probe_video, synth_soundtrack


HAS_FFMPEG = shutil.which("ffmpeg") is not None
HAS_FFPROBE = shutil.which("ffprobe") is not None


class VideoEncoderInputTests(unittest.TestCase):
    def test_frame_pipe_writer_retries_partial_writes(self):
        class ShortWriter:
            def __init__(self):
                self.data = bytearray()
                self.calls = 0

            def write(self, remaining):
                self.calls += 1
                chunk = remaining[:3]
                self.data.extend(chunk)
                return len(chunk)

        writer = ShortWriter()
        payload = bytes(range(17))
        _write_all(writer, memoryview(payload))
        self.assertEqual(writer.data, payload)
        self.assertGreater(writer.calls, 1)

        class StalledWriter:
            def write(self, remaining):
                return 0

        with self.assertRaises(BrokenPipeError):
            _write_all(StalledWriter(), memoryview(b"frame"))

    def test_ffmpeg_command_sets_high_quality_crf(self):
        with tempfile.TemporaryDirectory() as directory:
            encoder = VideoEncoder(Path(directory) / "movie.mp4", 4, 4, 2, 2)
            with patch("src.reels.encoding.shutil.which", return_value="/usr/bin/ffmpeg"):
                command = encoder._command(Path(directory) / ".movie.tmp.mp4")
            crf_index = command.index("-crf")
            self.assertEqual(command[crf_index + 1], "14")
            self.assertIn("+faststart", command)
            # BT.709 conversion and tags keep phone playback colors faithful.
            self.assertEqual(command[command.index("-colorspace") + 1], "bt709")
            self.assertIn("out_color_matrix=bt709", command[command.index("-vf") + 1])

    def test_rejects_invalid_video_configuration(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = (
                ({"width": 5, "height": 4}, "even"),
                ({"width": 4, "height": 0}, "height"),
                ({"fps": True}, "fps"),
                ({"frame_count": 0}, "frame_count"),
                ({"preset": "invented"}, "preset"),
            )
            for kwargs, message in cases:
                with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, message):
                    options = dict(width=4, height=4, fps=2, frame_count=2)
                    options.update(kwargs)
                    VideoEncoder(root / "movie.mp4", **options)
            with self.assertRaisesRegex(ValueError, r"\.mp4"):
                VideoEncoder(root / "movie.mov", 4, 4, 2, 2)

    def test_existing_output_rejected_without_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "movie.mp4"
            output.write_bytes(b"existing")
            with self.assertRaises(FileExistsError):
                VideoEncoder(output, 4, 4, 2, 2)
            self.assertEqual(output.read_bytes(), b"existing")

    def test_frame_validation_and_underwrite_clean_temporary_file(self):
        if not HAS_FFMPEG:
            self.skipTest("FFmpeg is not installed")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "movie.mp4"
            with self.assertRaisesRegex(ValueError, "expected 2 frames"):
                with VideoEncoder(output, 4, 4, 2, 2) as writer:
                    writer.write(np.zeros((4, 4, 3), dtype=np.uint8))
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".*.tmp.mp4")), [])

            with self.assertRaisesRegex(ValueError, "shape"):
                with VideoEncoder(output, 4, 4, 2, 1) as writer:
                    writer.write(np.zeros((4, 6, 3), dtype=np.uint8))
            self.assertEqual(list(root.glob(".*.tmp.mp4")), [])

    def test_render_exception_cleans_up_partial_output(self):
        if not HAS_FFMPEG:
            self.skipTest("FFmpeg is not installed")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "cancelled.mp4"
            with self.assertRaisesRegex(RuntimeError, "render stopped"):
                with VideoEncoder(output, 4, 4, 2, 3) as writer:
                    writer.write(np.zeros((4, 4, 3), dtype=np.uint8))
                    raise RuntimeError("render stopped")
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".*.tmp.mp4")), [])

    def test_missing_ffmpeg_has_install_hint_and_removes_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            encoder = VideoEncoder(root / "movie.mp4", 4, 4, 2, 1)
            with patch("src.reels.encoding.shutil.which", return_value=None):
                with self.assertRaisesRegex(FileNotFoundError, "install FFmpeg"):
                    encoder.__enter__()
            self.assertEqual(list(root.glob(".*.tmp.mp4")), [])

    def test_ffmpeg_error_is_reported_and_existing_output_survives(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake_ffmpeg = root / "ffmpeg-fail"
            fake_ffmpeg.write_text("#!/bin/sh\necho deliberate encoder failure >&2\nexit 17\n")
            fake_ffmpeg.chmod(0o755)
            output = root / "existing.mp4"
            output.write_bytes(b"preserve this")
            encoder = VideoEncoder(output, 4, 4, 1, 1, overwrite=True)
            with patch("src.reels.encoding.shutil.which", return_value=str(fake_ffmpeg)):
                with self.assertRaisesRegex(RuntimeError, "deliberate encoder failure"):
                    with encoder as writer:
                        writer.write(np.zeros((4, 4, 3), dtype=np.uint8))
            self.assertEqual(output.read_bytes(), b"preserve this")
            self.assertEqual(list(root.glob(".*.tmp.mp4")), [])


@unittest.skipUnless(HAS_FFMPEG and HAS_FFPROBE, "FFmpeg/FFprobe are not installed")
class FFmpegIntegrationTests(unittest.TestCase):
    def test_tiny_rgb_video_encodes_atomically_and_probes(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "tiny.mp4"
            frames = [
                np.full((4, 6, 3), (255, 0, 0), dtype=np.uint8),
                np.full((4, 6, 3), (0, 0, 255), dtype=np.uint8),
            ]
            with VideoEncoder(output, width=6, height=4, fps=2, frame_count=2) as writer:
                for frame in frames:
                    writer.write(frame[:, ::-1, :])  # exercise contiguous conversion
            self.assertTrue(output.is_file())
            self.assertEqual(list(Path(directory).glob(".*.tmp.mp4")), [])
            metadata = probe_video(output)
            video = next(stream for stream in metadata["streams"] if stream["codec_type"] == "video")
            self.assertEqual((video["width"], video["height"]), (6, 4))
            self.assertEqual(video["codec_name"], "h264")
            self.assertEqual(video["pix_fmt"], "yuv420p")
            self.assertEqual(int(video.get("nb_frames", video.get("nb_read_frames", 0))), 2)

    def test_soundtrack_is_deterministic_mastered_and_can_be_looped_into_video(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = synth_soundtrack(root / "first.wav", 0.3, seed=17)
            second = synth_soundtrack(root / "second.wav", 0.3, seed=17)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with wave.open(str(first), "rb") as wav_file:
                self.assertEqual(wav_file.getnchannels(), 2)
                self.assertEqual(wav_file.getframerate(), 48_000)
                self.assertEqual(wav_file.getnframes(), round(0.3 * 48_000))
                samples = np.frombuffer(wav_file.readframes(wav_file.getnframes()), dtype="<i2")
            self.assertGreater(np.max(np.abs(samples)), 0)
            # Mastered with 1 dB of headroom for the AAC encoder.
            self.assertLessEqual(np.max(np.abs(samples)), round(32767 * 10 ** (-1 / 20)) + 1)

            output = root / "with-audio.mp4"
            with VideoEncoder(output, 6, 4, 2, 2, audio=first) as writer:
                writer.write(np.zeros((4, 6, 3), dtype=np.uint8))
                writer.write(np.full((4, 6, 3), 128, dtype=np.uint8))
            metadata = probe_video(output)
            self.assertTrue(any(stream["codec_type"] == "audio" for stream in metadata["streams"]))
            self.assertTrue(any(stream["codec_type"] == "video" for stream in metadata["streams"]))

    def test_probe_rejects_nonvideo_files_with_ffprobe_error(self):
        with tempfile.TemporaryDirectory() as directory:
            invalid = Path(directory) / "not-video.mp4"
            invalid.write_text("not an mp4")
            with self.assertRaisesRegex(RuntimeError, "FFprobe failed"):
                probe_video(invalid)


if __name__ == "__main__":
    unittest.main()
