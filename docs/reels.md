# Educational science reels

The reel command turns one selected example into one vertical video. For the
eleven CFD examples, numerical adapters run the corresponding model and render
its computed field with six concise science-story captions curated in
[`src/reels/stories.py`](../src/reels/stories.py). Other examples use their
original VTK scene with a generic story, or a custom JSON story when supplied.

## Install and run

Start from the repository root. Install the reel dependencies and make sure the
`ffmpeg` and `ffprobe` executables are available on your system `PATH`:

```bash
python -m pip install -r requirements-reels.txt
ffmpeg -version
ffprobe -version
python -m src.reels --list
```

With no arguments, the command renders the obstacle-flow story using the
default studio quality and 60-second duration. Choose a fast render explicitly
when turnaround matters more than resolution:

```bash
python -m src.reels obstacle_flow_simulation --quality fast \
  --output reels/obstacle.mp4
```

Render one separate video per curated CFD example:

```bash
python -m src.reels --all cfd
```

Single-scene reels default to `reels/<example-stem>.mp4`. A batch writes one
MP4 per curated example under `reels/` or the selected output directory. Each reel also
gets `<stem>.jpg` (thumbnail), `<stem>.preview.jpg` (contact sheet),
`<stem>.json` (manifest), and `<stem>.srt` (captions) alongside the video.
Use `--output PATH` to choose the MP4 path for one scene or the destination
directory for a batch.

The command reports measured progress and elapsed generation time. A requested
60-second reel is a video-duration target, not a promise that rendering takes
60 seconds; generation time depends on the scene, machine, and quality setting.
Use `--overwrite` to replace existing outputs. `--no-cache` forces a fresh
scene capture; `--preview-only` creates a contact sheet without encoding video.

## Options

| Option | Purpose |
| --- | --- |
| `--list` | List scenes and their available curated stories. |
| `--all cfd` | Generate reels for the CFD story catalog. |
| `--duration SECONDS` | Set video duration (60 seconds by default). |
| `--quality studio` | Default: render at 1080 × 1920 and 30 frames per second. |
| `--quality fast` | Render at 720 × 1280 and 24 frames per second. |
| `--output PATH` | Choose the output MP4 path, or the output location for a batch. |
| `--silent` | Render without an audio track. |
| `--audio PATH` | Use an audio file in place of the generated soundtrack. |
| `--overwrite` | Replace an existing output. |
| `--story PATH.json` | Supply a custom story for the selected example. |
| `--no-cache` | Capture a fresh scene instead of reusing cached output. |
| `--preview-only` | Generate a contact sheet without an MP4. |

For example, render a single story at studio quality with an audio file:

```bash
python -m src.reels fixed_end_heat_transfer_1d --duration 45 \
  --quality studio --audio music.wav --output reels/heat.mp4
```

Use `--silent` for no audio:

```bash
python -m src.reels heat_convection_solver_2d --duration 45 \
  --quality fast --silent --output reels/scalar-transport.mp4
```

If the destination exists, add `--overwrite`. `--silent` and `--audio PATH`
are mutually exclusive.

### Custom story JSON

`--story` accepts an object with optional `title`, `hook`, `takeaway`,
`model_note`, `units`, and `accent` fields. `accent` is an RGB array with
integer components from 0 to 255. Supply exactly six beats, each containing a
`label`, `title`, and `body`. Missing top-level fields inherit their value
from the selected example's story when one exists. The selected example stem
remains the story key.

```json
{
  "title": "Follow the Heat",
  "hook": "Where does the energy leave?",
  "units": "°C",
  "accent": [70, 220, 235],
  "beats": [
    {"label": "Set the field", "title": "A warm boundary", "body": "The boundary temperature is prescribed."},
    {"label": "Add energy", "title": "A source inside", "body": "Heat generation warms the interior."},
    {"label": "Solve", "title": "Balance each cell", "body": "The numerical field settles to a steady balance."},
    {"label": "Read color", "title": "Temperature along x", "body": "Color shows the computed temperature field."},
    {"label": "Follow flux", "title": "Energy moves outward", "body": "Heat crosses the boundaries that can carry it away."},
    {"label": "Check", "title": "Close the budget", "body": "At steady state, input and output heat balance."}
  ]
}
```

For example: `python -m src.reels some_example --story story.json`. A custom
story gives a generic scene precise educational captions. For non-curated
examples, the fallback captures the selected example's original VTK scene; it
needs that example's optional dependencies and does not apply to converters or
projects without a scene. It retains all original viewports and scalar legends;
camera motion supplies the presentation animation.

## Look and sound

Scenes render at 3× (studio) or 2× (fast) supersampling on a transparent
background, then float over a designed backdrop with a soft glow. A perspective
camera is fitted tightly to each model and eases to a new pose in every
chapter. Headlines animate in word by word, and the profile chart draws itself.

The default soundtrack is generated per example (deterministic from the source
hash): a pad, pluck arpeggio, bass, lead, and drums at 48 kHz. Its tempo is
chosen so each chapter spans whole bars, and the arrangement builds chapter by
chapter, with risers and impacts on caption changes. Every soundtrack, including
`--audio` files, is loudness-normalized to −14 LUFS and faded out at the end.
Video is H.264 High profile with BT.709 color tags.

## Headless Linux

VTK needs a display to render its scene. On a Linux machine without a desktop
session, install Xvfb. The command detects the missing display and uses Xvfb
automatically. The listing command does not render a scene and does not need a
display. On macOS or Windows, run in a desktop session with VTK's rendering
dependencies installed. For non-CFD examples that import UI or integration
packages, install the full `requirements.txt` set as well.

```bash
xvfb-run -a python -m src.reels obstacle_flow_simulation \
  --duration 60 --quality fast --output reels/obstacle.mp4
```
This explicit wrapper is useful if you want to manage the virtual display
yourself; it is not required when `xvfb-run` is on the system `PATH`.

## What the reel represents

Curated reels use numerical adapters that call the example solvers, not a
recording of each script's original layout. The generic fallback instead
captures the original VTK scene. Adapter parameters and grid sizes are recorded
in the manifest. Quality changes both video resolution/frame rate and the
curated adapter grid, so fast and studio versions can have slightly different
computed fields. The curated obstacle reel uses `Re=50`; the standalone
obstacle example defaults to `Re=100`.

Reels retain each model's units. Heat and scalar transport adapters show steady
fields while captions advance; illustrative particles do not make them
transient. The obstacle adapter shows simulated snapshots in lattice units,
not physical seconds. Channel velocity is prescribed and its adapter computes
temperature; it does not solve momentum or pressure. Quality changes video
resolution, frame rate, and curated adapter grid resolution, but not the model
equations or boundary conditions. The manifest records the selected parameters.

Color legends retain the scalar ranges used by the scene, and captions avoid
implying more than the underlying model computes.
