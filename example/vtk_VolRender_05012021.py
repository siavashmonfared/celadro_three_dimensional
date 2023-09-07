#!/usr/bin/env python

import vtk
import numpy as np

def VolRendering(fin,fout):

    colors = vtk.vtkNamedColors()

    # This is a simple volume rendering example that
    # uses a vtkFixedPointVolumeRayCastMapper

    # Create the standard renderer, render window
    # and interactor.
    ren1 = vtk.vtkRenderer()

    #renWin = vtk.vtkRenderWindow()
    #renWin.AddRenderer(ren1)

    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)

    # Create the reader for the data.
    reader = vtk.vtkStructuredPointsReader()
    #reader.SetFileName("tissue.vtk")
    reader.SetFileName(fin)

    # Create transfer mapping scalar value to opacity.
    opacityTransferFunction = vtk.vtkPiecewiseFunction()
    opacityTransferFunction.AddPoint(0.5, 0.)
    opacityTransferFunction.AddPoint(1.0, 1.1)

    # Create transfer mapping scalar value to color.
    colorTransferFunction = vtk.vtkColorTransferFunction()
    colorTransferFunction.AddRGBPoint(0., 0.82421875, 0.82421875, 0.82421875);


    # The property describes how the data will look.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(colorTransferFunction)
    volumeProperty.SetScalarOpacity(opacityTransferFunction)
    volumeProperty.ShadeOn()
    volumeProperty.SetInterpolationTypeToLinear()

    # The mapper / ray cast function know how to render the data.
    #volumeMapper = vtk.vtkFixedPointVolumeRayCastMapper()
    #volumeMapper.SetInputConnection(reader.GetOutputPort())
    
    volumeMapper = vtk.vtkOpenGLGPUVolumeRayCastMapper();
    volumeMapper.SetInputConnection(reader.GetOutputPort())

    # The volume holds the mapper and the property and
    # can be used to position/orient the volume.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    ren1.AddVolume(volume)
    ren1.SetBackground(colors.GetColor3d('White'))
    #ren1.GetActiveCamera().Azimuth(45)
    #ren1.GetActiveCamera().Elevation(30)
    #ren1.ResetCameraClippingRange()
    #ren1.ResetCamera()

    #renWin.SetSize(750,750);
    #renWin.SetWindowName('SimpleRayCast')
    #renWin.Render()

    #iren.Start()
    
    camera = vtk.vtkCamera()
    ren1.SetActiveCamera(camera) 
    camera.SetViewUp(0,0,1);
    camera.SetPosition(32,300,200);
    camera.SetFocalPoint(32,32,0);
    camera.Zoom(1)
    camera.SetClippingRange(.1, 12000); 
    # self.camera.SetClippingRange(0.0, 100000)

    
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(ren1)
    renderWindow.SetOffScreenRendering(1)
    renderWindow.SetPosition(+1,+1);
    renderWindow.SetSize(750,750);
    renderWindow.Render();
    
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow);
    windowToImageFilter.Update();
    
    writer = vtk.vtkPNGWriter()
    #writer.SetFileName("config.png")
    writer.SetFileName(fout)
    writer.SetInputData(windowToImageFilter.GetOutput())
    writer.Write() 



print('running vtk_VolRender_05012021.py')
rng = np.arange(1,101,1)
#rng = [165,195,215,300]
for fr in (rng):
    print('time step: ',fr,flush=True)
    fin = 'frame_' + str(fr) + '.vtk'
    fout = 'config_' + str(fr) + '.png'
    VolRendering(fin,fout)

    

    
  
