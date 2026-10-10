# VTK.js basic example

This small browser demo connects a `vtkConeSource` to a mapper and actor, then displays the cone in a VTK.js render window. Use the **Base detail** slider to change the cone's polygon count. Drag the model to rotate it, scroll or pinch to zoom, and select **Reset view** to restore the camera.

## Requirements

- Node.js 22.15 or newer
- npm

## Run locally

From this directory, install the dependencies and start the development server:

```sh
npm ci
npm start
```

Open the local URL printed by webpack-dev-server (usually <http://localhost:8080>). The server recompiles the example when source files change.

To create a production bundle in `dist/`, run:

```sh
npm run build
```

The generated bundle is `dist/bundle.js`. Serve `vtk_js/basic/`, which contains `index.html`, over HTTP to view the build. For example, from this directory:

```sh
python -m http.server 8080 --bind 127.0.0.1
```

Open <http://localhost:8080>. Keep `index.html` and `dist/bundle.js` together with this folder structure when sharing the build.

## How it works

`src/main.js` creates a cone source, connects it to a mapper and actor, adds the actor to a renderer, and renders the scene. The slider changes the source resolution and triggers a render. The VTK.js fullscreen render helper tracks viewport changes so the canvas follows the responsive layout.
