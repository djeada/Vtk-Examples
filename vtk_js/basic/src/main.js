import '@kitware/vtk.js/Rendering/Profiles/Geometry';  
import vtkFullScreenRenderWindow from '@kitware/vtk.js/Rendering/Misc/FullScreenRenderWindow';  
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';  
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';  
import vtkCalculator from '@kitware/vtk.js/Filters/General/Calculator';  
import vtkConeSource from '@kitware/vtk.js/Filters/Sources/ConeSource';  

const fullScreenRenderWindow = vtkFullScreenRenderWindow.newInstance();
const renderer = fullScreenRenderWindow.getRenderer();
const renderWindow = fullScreenRenderWindow.getRenderWindow();

const coneSource = vtkConeSource.newInstance();
const mapper = vtkMapper.newInstance();
mapper.setInputConnection(coneSource.getOutputPort());

const actor = vtkActor.newInstance();
actor.setMapper(mapper);

renderer.addActor(actor);
renderer.resetCamera();
renderWindow.render();

