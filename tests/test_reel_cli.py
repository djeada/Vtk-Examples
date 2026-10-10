"""Pure command-line contract tests for educational reels."""

from contextlib import redirect_stderr, redirect_stdout
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from src.reels import __main__ as cli
from src.reels.stories import Beat


class ReelStoryTests(unittest.TestCase):
    def test_curated_and_generic_default_stories_have_six_beats(self):
        curated = cli.load_story("fixed_end_heat_transfer_1d")
        generic = cli.load_story("circle")

        self.assertEqual(curated.key, "fixed_end_heat_transfer_1d")
        self.assertEqual(len(curated.beats), 6)
        self.assertEqual(generic.key, "circle")
        self.assertEqual(len(generic.beats), 6)
        self.assertTrue(generic.model_note.strip())
        self.assertEqual(generic.units, "")

    def test_custom_story_validates_beats_and_inherits_omitted_fields(self):
        beats = [
            {"label": f"Part {i}", "title": f"Heading {i}", "body": f"Caption {i}"}
            for i in range(6)
        ]
        with tempfile.TemporaryDirectory() as directory:
            story_path = Path(directory) / "story.json"
            story_path.write_text(
                json.dumps(
                    {
                        "title": "Custom title",
                        "hook": "Custom hook",
                        "accent": [10, 20, 30],
                        "beats": beats,
                    }
                )
            )
            story = cli.load_story("fixed_end_heat_transfer_1d", story_path)

        curated = cli.load_story("fixed_end_heat_transfer_1d")
        self.assertEqual(story.title, "Custom title")
        self.assertEqual(story.hook, "Custom hook")
        self.assertEqual(story.accent, (10, 20, 30))
        self.assertEqual(story.model_note, curated.model_note)
        self.assertEqual(story.beats, tuple(Beat(**beat) for beat in beats))

    def test_custom_story_rejects_invalid_shapes_and_values(self):
        invalid_stories = (
            {"unknown": "field"},
            {"title": ""},
            {"accent": [0, 0, 256]},
            {"beats": [{"label": "only one"}]},
            {
                "beats": [
                    {"label": "Label", "title": "Title", "body": "Body"}
                    for _ in range(5)
                ]
            },
        )
        with tempfile.TemporaryDirectory() as directory:
            story_path = Path(directory) / "story.json"
            for data in invalid_stories:
                with self.subTest(data=data):
                    story_path.write_text(json.dumps(data))
                    with self.assertRaises(ValueError):
                        cli.load_story("circle", story_path)

    def test_subtitles_use_hook_first_and_takeaway_last_with_exact_timing(self):
        story = cli.load_story("fixed_end_heat_transfer_1d")
        subtitles = cli.subtitle_text(story, 60).split("\n\n")

        self.assertEqual(len(subtitles), 6)
        self.assertEqual(
            subtitles[0],
            "1\n00:00:00,000 --> 00:00:10,000\n"
            f"{story.title}\n{story.hook}",
        )
        self.assertEqual(
            subtitles[-1],
            "6\n00:00:50,000 --> 00:01:00,000\n"
            f"{story.beats[-1].title}\n{story.takeaway}\n",
        )


class ReelCommandTests(unittest.TestCase):
    def setUp(self):
        self.display = patch.dict(os.environ, {"DISPLAY": ":test-display"})
        self.display.start()
        self.addCleanup(self.display.stop)
        self.which = patch.object(cli.shutil, "which", return_value="/usr/bin/ffmpeg")
        self.which.start()
        self.addCleanup(self.which.stop)

    def test_list_needs_no_encoder_or_model_solver(self):
        output = io.StringIO()
        with (
            redirect_stdout(output),
            patch.object(cli.shutil, "which", side_effect=AssertionError("encoder lookup")),
            patch.object(cli, "discover_examples", side_effect=AssertionError("solver lookup")),
            patch.object(cli.subprocess, "call", side_effect=AssertionError("subprocess")),
        ):
            result = cli.main(["--list", "heat"])

        self.assertEqual(result, 0)
        listed = output.getvalue()
        self.assertIn("fixed_end_heat_transfer_1d", listed)
        self.assertIn("heat_convection_solver_2d", listed)
        self.assertNotIn("fluid_flow_simulator", listed)
        self.assertIn("Other native VTK examples", listed)

    def test_all_and_custom_story_are_rejected_before_work(self):
        stderr = io.StringIO()
        with (
            redirect_stderr(stderr),
            patch.object(cli, "discover_examples", side_effect=AssertionError("discovery")),
        ):
            with self.assertRaises(SystemExit) as error:
                cli.main(["--all", "cfd", "--story", "story.json"])

        self.assertEqual(error.exception.code, 2)
        self.assertIn("--all cannot be combined", stderr.getvalue())

    def test_invalid_output_selectors_and_example_return_clear_errors(self):
        output = io.StringIO()
        with redirect_stderr(output), patch.object(cli, "render_reel") as render:
            single_result = cli.main(["circle", "--output", "out-dir"])
            batch_result = cli.main(["--all", "cfd", "--output", "all.mp4"])
            example_result = cli.main(["does-not-exist"])

        self.assertEqual((single_result, batch_result, example_result), (2, 2, 2))
        self.assertEqual(render.call_count, 0)
        self.assertIn("Single-example --output must be an .mp4 path", output.getvalue())
        self.assertIn("Batch --output must be a directory", output.getvalue())
        self.assertIn("No example matches", output.getvalue())

    def test_single_and_batch_select_one_distinct_mp4_per_example(self):
        with tempfile.TemporaryDirectory() as directory:
            batch_dir = Path(directory) / "batch"
            with patch.object(cli, "render_reel") as render:
                single_result = cli.main(["circle"])
                single_path, single_output, single_args = render.call_args.args
                self.assertEqual(single_path.stem, "circle")
                self.assertEqual(single_output, Path("reels/circle.mp4"))
                self.assertEqual(single_args.quality, "studio")
                self.assertEqual(single_args.duration, 60)

                render.reset_mock()
                batch_result = cli.main(["--all", "cfd", "--output", str(batch_dir)])
                calls = render.call_args_list

        self.assertEqual((single_result, batch_result), (0, 0))
        self.assertEqual(len(calls), 9)
        paths = [call.args[0] for call in calls]
        outputs = [call.args[1] for call in calls]
        self.assertEqual({path.stem for path in paths}, set(cli.STORIES))
        self.assertEqual(len(set(outputs)), 9)
        self.assertTrue(all(path.parent == batch_dir for path in outputs))
        self.assertEqual({path.suffix for path in outputs}, {".mp4"})


if __name__ == "__main__":
    unittest.main()
