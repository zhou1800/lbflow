#### import the simple module from the paraview
from paraview.simple import *
from os import listdir
import gc
#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()

gc.enable()

inPath = './'
outFile = 'G.tif'
fNames = listdir(inPath)
fNames  = [x for x in fNames if "vti" in x]

for i in range(0, len(fNames)):
    fName = inPath + fNames[i]

    # create a new 'XML Image Data Reader'
    G = XMLImageDataReader(FileName = fName);
    G.PointArrayStatus = ['r0', 'r1']
    
    # get animation scene
    animationScene1 = GetAnimationScene()
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    renderView1.ViewSize = [1000, 1000]
    
    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [1133.0, 1333.5, 10000.0]
    renderView1.CameraFocalPoint = [1133.0, 1333.5, 0.0]
    renderView1.CameraParallelScale = 1446.1419649825814
    renderView1.OrientationAxesVisibility = 0
    
    # create a new 'Calculator'
    calculator1 = Calculator(Input=G)
    calculator1.Function = '(r0-r1)/(r0+r1)'
    calculator1.ReplacementValue = 3.0
    
    # get color transfer function/color map for 'Result'
    resultLUT = GetColorTransferFunction('Result')
    resultLUT.RGBPoints = [-1.0, 0.231373, 0.298039, 0.752941, 1.0, 0.7058823529411765, 0.01568627450980392, 0.14901960784313725, 3.0, 0.0, 0.0, 0.0]
    resultLUT.ScalarRangeInitialized = 1.0
    
    # show data in view
    calculator1Display = Show(calculator1, renderView1)
    # trace defaults for the display properties.
    calculator1Display.Representation = 'Slice'
    calculator1Display.ColorArrayName = ['POINTS', 'Result']
    calculator1Display.LookupTable = resultLUT
    calculator1Display.OSPRayScaleArray = 'Result'
    calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    calculator1Display.SelectOrientationVectors = 'None'
    calculator1Display.ScaleFactor = 266.7
    calculator1Display.SelectScaleArray = 'Result'
    calculator1Display.GlyphType = 'Arrow'
    calculator1Display.PolarAxes = 'PolarAxesRepresentation'
    calculator1Display.ScalarOpacityUnitDistance = 19.21315358834187
    calculator1Display.SliceMode = 'XZ Plane'
    calculator1Display.Slice = 0
    
    # hide color bar/color legend
    calculator1Display.SetScalarBarVisibility(renderView1, False)
    
    # save screenshot
    #for i in range(start, end):
    SaveScreenshot(fName[:-4] + '.tif', magnification=3, quality=100, view=renderView1)
    #    animationScene1.GoToNext()

    # save animation images/movie
    #WriteAnimation(fNames[i][:-14] + str(start) + '.tif', Magnification=3, FrameRate=15.0, Compression=False)

    Delete(renderView1)
    del renderView1
    gc.collect()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
