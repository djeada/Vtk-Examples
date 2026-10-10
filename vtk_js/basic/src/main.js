import '@kitware/vtk.js/Rendering/Profiles/Geometry';
import vtkFullScreenRenderWindow from '@kitware/vtk.js/Rendering/Misc/FullScreenRenderWindow';
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkConeSource from '@kitware/vtk.js/Filters/Sources/ConeSource';

const viewport = document.querySelector('#viewport');
const resolutionInput = document.querySelector('#resolution');
const resolutionValue = document.querySelector('#resolution-value');
const resetButton = document.querySelector('#reset-view');

const fullScreenRenderWindow = vtkFullScreenRenderWindow.newInstance({
  rootContainer: viewport,
  background: [0.055, 0.09, 0.15],
  containerStyle: {
    width: '100%',
    height: '100%',
    position: 'absolute',
    top: '0',
    left: '0',
    overflow: 'hidden',
  },
});

const renderer = fullScreenRenderWindow.getRenderer();
const renderWindow = fullScreenRenderWindow.getRenderWindow();

const coneSource = vtkConeSource.newInstance({ resolution: Number(resolutionInput.value) });
const mapper = vtkMapper.newInstance();
mapper.setInputConnection(coneSource.getOutputPort());

const actor = vtkActor.newInstance();
actor.setMapper(mapper);
actor.getProperty().setColor(0.96, 0.42, 0.22);

renderer.addActor(actor);
renderer.resetCamera();
renderWindow.render();

const camera = renderer.getActiveCamera();
const initialCamera = {
  position: camera.getPosition(),
  focalPoint: camera.getFocalPoint(),
  viewUp: camera.getViewUp(),
};

fullScreenRenderWindow.setResizeCallback(({ width, height }) => {
  camera.setUseHorizontalViewAngle(width < height);
});

resolutionInput.addEventListener('input', () => {
  const resolution = Number(resolutionInput.value);
  coneSource.setResolution(resolution);
  resolutionValue.value = resolution;
  resolutionValue.textContent = resolution;
  renderWindow.render();
});

resetButton.addEventListener('click', () => {
  camera.setPosition(...initialCamera.position);
  camera.setFocalPoint(...initialCamera.focalPoint);
  camera.setViewUp(...initialCamera.viewUp);
  renderer.resetCameraClippingRange();
  renderWindow.render();
});
