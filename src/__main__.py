"""Discover and launch examples with ``python -m src``."""

import argparse
from pathlib import Path
import subprocess
import sys


SOURCE_DIR = Path(__file__).resolve().parent


def discover_examples(source_dir=SOURCE_DIR):
    """Return example paths without importing code or starting a GUI."""
    return sorted(
        path
        for category in source_dir.iterdir()
        if category.is_dir() and category.name[:2].isdigit()
        for path in category.glob("*.py")
        if not path.name.startswith("_") and path.stem != "converter_interface"
    )


def select_example(selector, examples, source_dir=SOURCE_DIR):
    """Prefer exact names/paths; otherwise require one substring match."""
    selector = selector.removeprefix("src/").removesuffix(".py").lower()
    exact = [
        path
        for path in examples
        if selector
        in {
            path.stem.lower(),
            path.relative_to(source_dir).with_suffix("").as_posix().lower(),
        }
    ]
    matches = exact or [
        path
        for path in examples
        if selector in path.relative_to(source_dir).as_posix().lower()
    ]
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError(
            f"No example matches {selector!r}. Use --list to browse examples."
        )
    choices = "\n  ".join(
        path.relative_to(source_dir).with_suffix("").as_posix() for path in matches
    )
    raise ValueError(
        f"Multiple examples match {selector!r}; use a full name:\n  {choices}"
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Browse and run VTK examples from the repository root.",
        epilog="Example: python -m src naca_airfoil_generator -- 0012 --export airfoil.vtp",
    )
    parser.add_argument(
        "example", nargs="?", help="name, category/path, or unique substring"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list examples (optionally filtered by example)",
    )
    arguments = list(sys.argv[1:] if argv is None else argv)
    script_arguments = []
    if "--" in arguments:
        separator = arguments.index("--")
        script_arguments = arguments[separator + 1 :]
        arguments = arguments[:separator]
    args = parser.parse_args(arguments)
    if script_arguments and (args.list or not args.example):
        parser.error("arguments after -- require an example to run")

    examples = discover_examples()
    if args.list:
        for path in examples:
            name = path.relative_to(SOURCE_DIR).with_suffix("").as_posix()
            if not args.example or args.example.lower() in name.lower():
                print(name)
        return 0
    if not args.example:
        parser.print_help()
        return 0
    try:
        path = select_example(args.example, examples)
    except ValueError as error:
        parser.error(str(error))
    # A separate process preserves each script's direct-run and sibling imports.
    # Use the active interpreter and preserve the caller's output directory.
    try:
        return subprocess.call([sys.executable, str(path), *script_arguments])
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
