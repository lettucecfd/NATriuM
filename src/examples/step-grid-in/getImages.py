# trace generated using paraview version 5.13.3
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview import simple as pvs
#### disable automatic camera reset on 'Show'
pvs._DisableFirstRenderCameraReset()

iT=20000
foldername="9682930_standard"
# create a new 'XML Partitioned Unstructured Grid Reader'
t_0200000pvtu = pvs.XMLPartitionedUnstructuredGridReader(registrationName=f't_0.{iT}.pvtu', FileName=['C:\\Users\\phili\\Desktop\\{foldername}\\output\\t_0.{iT}.pvtu'])

# Properties modified on t_0200000pvtu
t_0200000pvtu.TimeArray = 'None'

# get active view
renderView1 = pvs.GetActiveViewOrCreate('RenderView')

# show data in view
t_0200000pvtuDisplay = pvs.Show(t_0200000pvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
t_0200000pvtuDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.5, 0.0, 107.2]
renderView1.CameraFocalPoint = [1.5, 0.0, 0.0]

# get the material library
materialLibrary1 = pvs.GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
pvs.ColorBy(t_0200000pvtuDisplay, ('POINTS', 'T'))

# rescale color and/or opacity maps used to include current data range
t_0200000pvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
t_0200000pvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'T'
tLUT = pvs.GetColorTransferFunction('T')

# get opacity transfer function/opacity map for 'T'
tPWF = pvs.GetOpacityTransferFunction('T')

# get 2D transfer function for 'T'
tTF2D = pvs.GetTransferFunction2D('T')

# get layout
layout1 = pvs.GetLayout()

# layout/tab size in pixels
layout1.SetSize(2082, 976)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.305108721842908, -0.017864360196445533, 107.2]
renderView1.CameraFocalPoint = [5.305108721842908, -0.017864360196445533, 0.0]
renderView1.CameraParallelScale = 4.376768248129167

# save screenshot
pvs.SaveScreenshot(filename=f'C:/Users/phili/Desktop/{foldername}/images/rhoScreenshot_iT{iT}.png', viewOrLayout=renderView1, location=16, ImageResolution=[2082, 976])

# create a new 'Calculator'
calculator1 = pvs.Calculator(registrationName='Calculator1', Input=t_0200000pvtu)

# Properties modified on calculator1
calculator1.ResultArrayName = 'U'
calculator1.Function = 'iHat*ux+jHat*uy'

# show data in view
calculator1Display = pvs.Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'

# hide data in view
pvs.Hide(t_0200000pvtu, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Calculator'
calculator2 = pvs.Calculator(registrationName='Calculator2', Input=calculator1)

# Properties modified on calculator2
calculator2.ResultArrayName = 'Ma'
calculator2.Function = 'sqrt(ux*ux+uy*uy)*1.5*1.5*1.5*1.4*1/3'

# show data in view
calculator2Display = pvs.Show(calculator2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator2Display.Representation = 'Surface'

# hide data in view
pvs.Hide(calculator1, renderView1)

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Ma'
maLUT = pvs.GetColorTransferFunction('Ma')

# get opacity transfer function/opacity map for 'Ma'
maPWF = pvs.GetOpacityTransferFunction('Ma')

# get 2D transfer function for 'Ma'
maTF2D = pvs.GetTransferFunction2D('Ma')

# create a new 'Contour'
contour1 = pvs.Contour(registrationName='Contour1', Input=calculator2)

# Properties modified on contour1
contour1.Isosurfaces = [0.0, 0.03677437641404573, 0.07354875282809147, 0.11032312924213719, 0.14709750565618293, 0.18387188207022867, 0.22064625848427438, 0.2574206348983201, 0.29419501131236586, 0.3309693877264116, 0.36774376414045734, 0.4045181405545031, 0.44129251696854876, 0.4780668933825945, 0.5148412697966402, 0.551615646210686, 0.5883900226247317, 0.6251643990387774, 0.6619387754528232, 0.6987131518668689, 0.7354875282809147, 0.7722619046949604, 0.8090362811090062, 0.8458106575230518, 0.8825850339370975, 0.9193594103511433, 0.956133786765189, 0.9929081631792348, 1.0296825395932805, 1.0664569160073263, 1.103231292421372, 1.1400056688354177, 1.1767800452494634, 1.2135544216635092, 1.2503287980775548, 1.2871031744916006, 1.3238775509056464, 1.3606519273196922, 1.3974263037337378, 1.4342006801477836, 1.4709750565618294, 1.507749432975875, 1.5445238093899207, 1.5812981858039665, 1.6180725622180123, 1.654846938632058, 1.6916213150461037, 1.7283956914601495, 1.765170067874195, 1.8019444442882409]

# show data in view
contour1Display = pvs.Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'

# hide data in view
pvs.Hide(calculator2, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# layout/tab size in pixels
layout1.SetSize(2086, 980)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.305108721842908, -0.017864360196445533, 107.2]
renderView1.CameraFocalPoint = [5.305108721842908, -0.017864360196445533, 0.0]
renderView1.CameraParallelScale = 4.376768248129167

# save screenshot
pvs.SaveScreenshot(filename=f'C:/Users/phili/Desktop/{foldername}/images/MaStreamlinesScreenshot_iT{iT}.png', viewOrLayout=renderView1, location=16, ImageResolution=[2086, 980])

# set active source
pvs.SetActiveSource(calculator2)

# show data in view
calculator2Display = pvs.Show(calculator2, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
calculator2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
pvs.Hide(contour1, renderView1)

# layout/tab size in pixels
layout1.SetSize(2086, 980)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.305108721842908, -0.017864360196445533, 107.2]
renderView1.CameraFocalPoint = [5.305108721842908, -0.017864360196445533, 0.0]
renderView1.CameraParallelScale = 4.376768248129167

# save screenshot
pvs.SaveScreenshot(filename=f'C:/Users/phili/Desktop/{foldername}/images/MaScreenshot_iT{iT}.png', viewOrLayout=renderView1, location=16, ImageResolution=[2086, 980])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2090, 984)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [5.305108721842908, -0.017864360196445533, 107.2]
renderView1.CameraFocalPoint = [5.305108721842908, -0.017864360196445533, 0.0]
renderView1.CameraParallelScale = 4.376768248129167


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------