# Contributing

Bug reports, fixes, new examples, and documentation updates are welcome. For larger changes, open an issue first to discuss the approach.

## Set up the project

Use Python 3.10 or newer. From the repository root, create a virtual environment and install either the common dependencies or the full set:

```bash
# Linux or macOS
python3 -m venv .venv
source .venv/bin/activate
# Windows PowerShell:
# py -m venv .venv
# .venv\Scripts\Activate.ps1
python -m pip install -r requirements-core.txt
```

Install `requirements.txt` instead when working on examples that use the optional Qt, Spark, or other integrations.

## Work on an example

Python examples live under `src/`, grouped by topic. Run them from the repository root through the launcher so shared files resolve consistently:

```bash
python -m src --list
python -m src circle
python -m src 01_basic_shapes/circle
```

Examples that use input files should use the shared files under `data/` where appropriate. Avoid committing generated output or machine-specific paths. Keep each example focused, and follow the naming and organization of nearby files.

## Check your changes

Run the relevant example and, when changing shared behavior, run the unit tests:

```bash
python -m pip install -r requirements-dev.txt
ruff check src tests
python -m compileall -q src tests
python -m unittest discover -s tests
```

Some VTK examples open a rendering window and need a graphical display. If an example requires a special dataset, UI toolkit, or runtime setup, document that requirement in its source or the README.

On Linux without a display, run the tests with `xvfb-run -a python -m unittest discover -s tests`. The checks workflow runs lint, Python compilation, these tests, and the VTK.js production build. For browser changes, also run `npm ci` and `npm run build` inside `vtk_js/basic`, and check the controls in a browser.

## Submit a change

1. Fork the repository and create a branch for your change.
2. Make a focused change and update documentation when behavior or setup changes.
3. Commit with a short message describing the change.
4. Push the branch and open a pull request with a brief description and the checks you ran.
