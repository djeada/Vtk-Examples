import contextlib
import io
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

from src import __main__ as launcher


class LauncherTests(unittest.TestCase):
    def test_discovery_excludes_helpers_and_finds_all_categories(self):
        examples = launcher.discover_examples()
        self.assertTrue(examples)
        self.assertEqual(len({path.parent for path in examples}), 9)
        self.assertNotIn("converter_interface", {path.stem for path in examples})
        self.assertEqual(examples, sorted(examples))

    def test_selection_accepts_names_paths_and_unique_fragments(self):
        examples = launcher.discover_examples()
        expected = launcher.SOURCE_DIR / "01_basic_shapes" / "circle.py"
        for selector in (
            "circle",
            "CIRCLE",
            "01_basic_shapes/circle",
            "src/01_basic_shapes/circle.py",
            "circl",
        ):
            with self.subTest(selector=selector):
                self.assertEqual(launcher.select_example(selector, examples), expected)

    def test_exact_name_wins_over_substring(self):
        examples = launcher.discover_examples()
        selected = launcher.select_example("structured_grid", examples)
        self.assertEqual(selected.stem, "structured_grid")

    def test_unknown_or_ambiguous_selection_is_actionable(self):
        examples = launcher.discover_examples()
        with self.assertRaisesRegex(ValueError, "No example matches.*--list"):
            launcher.select_example("unknown-example", examples)
        with self.assertRaisesRegex(ValueError, "Multiple examples.*full name"):
            launcher.select_example("heat", examples)

    def test_list_needs_no_vtk_and_starts_no_process(self):
        # -S omits installed packages, so this also verifies discovery is usable
        # before installing the optional GUI/scientific dependencies.
        completed = subprocess.run(
            [sys.executable, "-S", "-m", "src", "--list", "circle"],
            cwd=launcher.SOURCE_DIR.parent,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertEqual(completed.stdout.strip(), "01_basic_shapes/circle")

    def test_forwarding_preserves_interpreter_arguments_and_exit_code(self):
        with patch.object(launcher.subprocess, "call", return_value=7) as call:
            result = launcher.main(
                ["naca_airfoil_generator", "--", "0012", "--export", "my airfoil.vtp"]
            )
        self.assertEqual(result, 7)
        command = call.call_args.args[0]
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(Path(command[1]).stem, "naca_airfoil_generator")
        self.assertEqual(command[2:], ["0012", "--export", "my airfoil.vtp"])

    def test_no_selection_shows_help(self):
        with contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(launcher.main([]), 0)
        self.assertIn("--list", output.getvalue())


if __name__ == "__main__":
    unittest.main()
