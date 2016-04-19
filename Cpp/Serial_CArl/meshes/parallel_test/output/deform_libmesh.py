#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


nb_of_blocks = 30

carl_multi_crystal_test_macro_filename = '/Users/breubreubreu/Programming/CArl/Cpp/Mesh_intersection/meshes/3D/output/carl_multi_crystal_test_macro_w_0_2__p_20.exo'
carl_multi_crystal_test_micro_filename = '/Users/breubreubreu/Programming/CArl/Cpp/Mesh_intersection/meshes/3D/output/carl_multi_crystal_test_micro_w_0_2__p_20.exo'

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]

# create the two 'ExodusIIReader'
carl_multi_crystal_test_macro = ExodusIIReader(FileName=[carl_multi_crystal_test_macro_filename])
carl_multi_crystal_test_macro.GenerateObjectIdCellArray = 1
carl_multi_crystal_test_macro.GenerateGlobalElementIdArray = 1
carl_multi_crystal_test_macro.ElementVariables = []
carl_multi_crystal_test_macro.FaceVariables = []
carl_multi_crystal_test_macro.EdgeVariables = []
carl_multi_crystal_test_macro.SideSetResultArrayStatus = []
carl_multi_crystal_test_macro.NodeSetResultArrayStatus = []
carl_multi_crystal_test_macro.FaceSetResultArrayStatus = []
carl_multi_crystal_test_macro.EdgeSetResultArrayStatus = []
carl_multi_crystal_test_macro.GenerateGlobalNodeIdArray = 1
carl_multi_crystal_test_macro.ElementSetResultArrayStatus = []
carl_multi_crystal_test_macro.PointVariables = []
carl_multi_crystal_test_macro.GlobalVariables = []
carl_multi_crystal_test_macro.ApplyDisplacements = 1
carl_multi_crystal_test_macro.DisplacementMagnitude = 1.0
carl_multi_crystal_test_macro.EdgeBlocks = []
carl_multi_crystal_test_macro.NodeSetArrayStatus = []
carl_multi_crystal_test_macro.SideSetArrayStatus = []
carl_multi_crystal_test_macro.FaceSetArrayStatus = []
carl_multi_crystal_test_macro.EdgeSetArrayStatus = []
carl_multi_crystal_test_macro.ElementSetArrayStatus = []
carl_multi_crystal_test_macro.NodeMapArrayStatus = []
carl_multi_crystal_test_macro.EdgeMapArrayStatus = []
carl_multi_crystal_test_macro.FaceMapArrayStatus = []
carl_multi_crystal_test_macro.ElementMapArrayStatus = []
carl_multi_crystal_test_macro.ElementBlocks = []
carl_multi_crystal_test_macro.FaceBlocks = []
carl_multi_crystal_test_macro.HasModeShapes = 0
carl_multi_crystal_test_macro.ModeShape = 1
carl_multi_crystal_test_macro.AnimateVibrations = 1
carl_multi_crystal_test_macro.GenerateFileIdArray = 0

carl_multi_crystal_test_micro = ExodusIIReader(FileName=[carl_multi_crystal_test_micro_filename])
carl_multi_crystal_test_micro.GenerateObjectIdCellArray = 1
carl_multi_crystal_test_micro.GenerateGlobalElementIdArray = 1
carl_multi_crystal_test_micro.ElementVariables = []
carl_multi_crystal_test_micro.FaceVariables = []
carl_multi_crystal_test_micro.EdgeVariables = []
carl_multi_crystal_test_micro.SideSetResultArrayStatus = []
carl_multi_crystal_test_micro.NodeSetResultArrayStatus = []
carl_multi_crystal_test_micro.FaceSetResultArrayStatus = []
carl_multi_crystal_test_micro.EdgeSetResultArrayStatus = []
carl_multi_crystal_test_micro.GenerateGlobalNodeIdArray = 1
carl_multi_crystal_test_micro.ElementSetResultArrayStatus = []
carl_multi_crystal_test_micro.PointVariables = []
carl_multi_crystal_test_micro.GlobalVariables = []
carl_multi_crystal_test_micro.ApplyDisplacements = 1
carl_multi_crystal_test_micro.DisplacementMagnitude = 1.0
carl_multi_crystal_test_micro.EdgeBlocks = []
carl_multi_crystal_test_micro.NodeSetArrayStatus = []
carl_multi_crystal_test_micro.SideSetArrayStatus = []
carl_multi_crystal_test_micro.FaceSetArrayStatus = []
carl_multi_crystal_test_micro.EdgeSetArrayStatus = []
carl_multi_crystal_test_micro.ElementSetArrayStatus = []
carl_multi_crystal_test_micro.NodeMapArrayStatus = []
carl_multi_crystal_test_micro.EdgeMapArrayStatus = []
carl_multi_crystal_test_micro.FaceMapArrayStatus = []
carl_multi_crystal_test_micro.ElementMapArrayStatus = []
carl_multi_crystal_test_micro.ElementBlocks = []
carl_multi_crystal_test_micro.FaceBlocks = []
carl_multi_crystal_test_micro.HasModeShapes = 0
carl_multi_crystal_test_micro.ModeShape = 1
carl_multi_crystal_test_micro.AnimateVibrations = 1
carl_multi_crystal_test_micro.GenerateFileIdArray = 0

# Properties modified on carl_multi_crystal_test_micro
carl_multi_crystal_test_micro.ElementVariables = ['E', 'mu', 'sigma_00', 'sigma_01', 'sigma_02', 'sigma_10', 'sigma_11', 'sigma_12', 'sigma_20', 'sigma_21', 'sigma_22', 'vonMises']
carl_multi_crystal_test_micro.PointVariables = ['u', 'v', 'w']

element_block_list = range(nb_of_blocks)
for iii in xrange(nb_of_blocks):
    element_block_list[iii] = 'Unnamed block ID: ' + str(iii + 1) + ' Type: TETRA4'

carl_multi_crystal_test_micro.ElementBlocks = element_block_list

# Properties modified on carl_multi_crystal_test_macro
carl_multi_crystal_test_macro.ElementVariables = ['E', 'mu', 'sigma_00', 'sigma_01', 'sigma_02', 'sigma_10', 'sigma_11', 'sigma_12', 'sigma_20', 'sigma_21', 'sigma_22', 'vonMises']
carl_multi_crystal_test_macro.PointVariables = ['u', 'v', 'w']
carl_multi_crystal_test_macro.ElementBlocks = ['Unnamed block ID: 33 Type: TETRA4']

# create a new 'Calculator'
calculator1 = Calculator(Input=carl_multi_crystal_test_macro)
calculator1.AttributeMode = 'Point Data'
calculator1.CoordinateResults = 0
calculator1.ResultNormals = 0
calculator1.ResultTCoords = 0
calculator1.ResultArrayName = 'deformation'
calculator1.Function = 'u*iHat+v*jHat+w*kHat'
calculator1.ReplaceInvalidResults = 1
calculator1.ReplacementValue = 0.0

# create a new 'Calculator'
calculator2 = Calculator(Input=carl_multi_crystal_test_micro)
calculator2.AttributeMode = 'Point Data'
calculator2.CoordinateResults = 0
calculator2.ResultNormals = 0
calculator2.ResultTCoords = 0
calculator2.ResultArrayName = 'deformation'
calculator2.Function = 'u*iHat+v*jHat+w*kHat'
calculator2.ReplaceInvalidResults = 1
calculator2.ReplacementValue = 0.0

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(Input=calculator1)
warpByVector1.Vectors = ['POINTS', 'deformation']
warpByVector1.ScaleFactor = 1.0

# create a new 'Warp By Vector'
warpByVector2 = WarpByVector(Input=calculator2)
warpByVector2.Vectors = ['POINTS', 'deformation']
warpByVector2.ScaleFactor = 1.0

# show data in view
warpByVector1Display = Show(warpByVector1)
# trace defaults for the display properties.
warpByVector1Display.CubeAxesVisibility = 0
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.AmbientColor = [1.0, 1.0, 1.0]
warpByVector1Display.ColorArrayName = [None, '']
warpByVector1Display.DiffuseColor = [1.0, 1.0, 1.0]
warpByVector1Display.LookupTable = None
warpByVector1Display.MapScalars = 1
warpByVector1Display.InterpolateScalarsBeforeMapping = 1
warpByVector1Display.Opacity = 1.0
warpByVector1Display.PointSize = 2.0
warpByVector1Display.LineWidth = 1.0
warpByVector1Display.Interpolation = 'Gouraud'
warpByVector1Display.Specular = 1.0
warpByVector1Display.SpecularColor = [1.0, 1.0, 1.0]
warpByVector1Display.SpecularPower = 100.0
warpByVector1Display.Ambient = 0.0
warpByVector1Display.Diffuse = 1.0
warpByVector1Display.EdgeColor = [0.0, 0.0, 0.0]
warpByVector1Display.BackfaceRepresentation = 'Follow Frontface'
warpByVector1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
warpByVector1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
warpByVector1Display.BackfaceOpacity = 1.0
warpByVector1Display.Position = [0.0, 0.0, 0.0]
warpByVector1Display.Scale = [1.0, 1.0, 1.0]
warpByVector1Display.Orientation = [0.0, 0.0, 0.0]
warpByVector1Display.Origin = [0.0, 0.0, 0.0]
warpByVector1Display.Pickable = 1
warpByVector1Display.Texture = None
warpByVector1Display.Triangulate = 0
warpByVector1Display.NonlinearSubdivisionLevel = 1
warpByVector1Display.CubeAxesColor = [1.0, 1.0, 1.0]
warpByVector1Display.CubeAxesCornerOffset = 0.0
warpByVector1Display.CubeAxesFlyMode = 'Closest Triad'
warpByVector1Display.CubeAxesInertia = 1
warpByVector1Display.CubeAxesTickLocation = 'Inside'
warpByVector1Display.CubeAxesXAxisMinorTickVisibility = 1
warpByVector1Display.CubeAxesXAxisTickVisibility = 1
warpByVector1Display.CubeAxesXAxisVisibility = 1
warpByVector1Display.CubeAxesXGridLines = 0
warpByVector1Display.CubeAxesXTitle = 'X-Axis'
warpByVector1Display.CubeAxesUseDefaultXTitle = 1
warpByVector1Display.CubeAxesYAxisMinorTickVisibility = 1
warpByVector1Display.CubeAxesYAxisTickVisibility = 1
warpByVector1Display.CubeAxesYAxisVisibility = 1
warpByVector1Display.CubeAxesYGridLines = 0
warpByVector1Display.CubeAxesYTitle = 'Y-Axis'
warpByVector1Display.CubeAxesUseDefaultYTitle = 1
warpByVector1Display.CubeAxesZAxisMinorTickVisibility = 1
warpByVector1Display.CubeAxesZAxisTickVisibility = 1
warpByVector1Display.CubeAxesZAxisVisibility = 1
warpByVector1Display.CubeAxesZGridLines = 0
warpByVector1Display.CubeAxesZTitle = 'Z-Axis'
warpByVector1Display.CubeAxesUseDefaultZTitle = 1
warpByVector1Display.CubeAxesGridLineLocation = 'All Faces'
warpByVector1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
warpByVector1Display.CustomBoundsActive = [0, 0, 0]
warpByVector1Display.OriginalBoundsRangeActive = [0, 0, 0]
warpByVector1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
warpByVector1Display.CustomRangeActive = [0, 0, 0]
warpByVector1Display.UseAxesOrigin = 0
warpByVector1Display.AxesOrigin = [0.0, 0.0, 0.0]
warpByVector1Display.CubeAxesXLabelFormat = '%-#6.3g'
warpByVector1Display.CubeAxesYLabelFormat = '%-#6.3g'
warpByVector1Display.CubeAxesZLabelFormat = '%-#6.3g'
warpByVector1Display.StickyAxes = 0
warpByVector1Display.CenterStickyAxes = 0
warpByVector1Display.SelectionCellLabelBold = 0
warpByVector1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
warpByVector1Display.SelectionCellLabelFontFamily = 'Arial'
warpByVector1Display.SelectionCellLabelFontSize = 18
warpByVector1Display.SelectionCellLabelItalic = 0
warpByVector1Display.SelectionCellLabelJustification = 'Left'
warpByVector1Display.SelectionCellLabelOpacity = 1.0
warpByVector1Display.SelectionCellLabelShadow = 0
warpByVector1Display.SelectionPointLabelBold = 0
warpByVector1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
warpByVector1Display.SelectionPointLabelFontFamily = 'Arial'
warpByVector1Display.SelectionPointLabelFontSize = 18
warpByVector1Display.SelectionPointLabelItalic = 0
warpByVector1Display.SelectionPointLabelJustification = 'Left'
warpByVector1Display.SelectionPointLabelOpacity = 1.0
warpByVector1Display.SelectionPointLabelShadow = 0
warpByVector1Display.ScalarOpacityUnitDistance = 0.46260204443706504
warpByVector1Display.SelectMapper = 'Projected tetra'

# change representation type
warpByVector1Display.SetRepresentationType('Wireframe')

# set scalar coloring
ColorBy(warpByVector1Display, ('POINTS', 'u'))

# rescale color and/or opacity maps used to include current data range
warpByVector1Display.RescaleTransferFunctionToDataRange(True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')

# show data in view
warpByVector2Display = Show(warpByVector2,renderView1)
# trace defaults for the display properties.
warpByVector2Display.CubeAxesVisibility = 0
warpByVector2Display.Representation = 'Surface'
warpByVector2Display.AmbientColor = [1.0, 1.0, 1.0]
warpByVector2Display.ColorArrayName = [None, '']
warpByVector2Display.DiffuseColor = [1.0, 1.0, 1.0]
warpByVector2Display.LookupTable = None
warpByVector2Display.MapScalars = 1
warpByVector2Display.InterpolateScalarsBeforeMapping = 1
warpByVector2Display.Opacity = 1.0
warpByVector2Display.PointSize = 2.0
warpByVector2Display.LineWidth = 1.0
warpByVector2Display.Interpolation = 'Gouraud'
warpByVector2Display.Specular = 1.0
warpByVector2Display.SpecularColor = [1.0, 1.0, 1.0]
warpByVector2Display.SpecularPower = 100.0
warpByVector2Display.Ambient = 0.0
warpByVector2Display.Diffuse = 1.0
warpByVector2Display.EdgeColor = [0.0, 0.0, 0.0]
warpByVector2Display.BackfaceRepresentation = 'Follow Frontface'
warpByVector2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
warpByVector2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
warpByVector2Display.BackfaceOpacity = 1.0
warpByVector2Display.Position = [0.0, 0.0, 0.0]
warpByVector2Display.Scale = [1.0, 1.0, 1.0]
warpByVector2Display.Orientation = [0.0, 0.0, 0.0]
warpByVector2Display.Origin = [0.0, 0.0, 0.0]
warpByVector2Display.Pickable = 1
warpByVector2Display.Texture = None
warpByVector2Display.Triangulate = 0
warpByVector2Display.NonlinearSubdivisionLevel = 1
warpByVector2Display.CubeAxesColor = [1.0, 1.0, 1.0]
warpByVector2Display.CubeAxesCornerOffset = 0.0
warpByVector2Display.CubeAxesFlyMode = 'Closest Triad'
warpByVector2Display.CubeAxesInertia = 1
warpByVector2Display.CubeAxesTickLocation = 'Inside'
warpByVector2Display.CubeAxesXAxisMinorTickVisibility = 1
warpByVector2Display.CubeAxesXAxisTickVisibility = 1
warpByVector2Display.CubeAxesXAxisVisibility = 1
warpByVector2Display.CubeAxesXGridLines = 0
warpByVector2Display.CubeAxesXTitle = 'X-Axis'
warpByVector2Display.CubeAxesUseDefaultXTitle = 1
warpByVector2Display.CubeAxesYAxisMinorTickVisibility = 1
warpByVector2Display.CubeAxesYAxisTickVisibility = 1
warpByVector2Display.CubeAxesYAxisVisibility = 1
warpByVector2Display.CubeAxesYGridLines = 0
warpByVector2Display.CubeAxesYTitle = 'Y-Axis'
warpByVector2Display.CubeAxesUseDefaultYTitle = 1
warpByVector2Display.CubeAxesZAxisMinorTickVisibility = 1
warpByVector2Display.CubeAxesZAxisTickVisibility = 1
warpByVector2Display.CubeAxesZAxisVisibility = 1
warpByVector2Display.CubeAxesZGridLines = 0
warpByVector2Display.CubeAxesZTitle = 'Z-Axis'
warpByVector2Display.CubeAxesUseDefaultZTitle = 1
warpByVector2Display.CubeAxesGridLineLocation = 'All Faces'
warpByVector2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
warpByVector2Display.CustomBoundsActive = [0, 0, 0]
warpByVector2Display.OriginalBoundsRangeActive = [0, 0, 0]
warpByVector2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
warpByVector2Display.CustomRangeActive = [0, 0, 0]
warpByVector2Display.UseAxesOrigin = 0
warpByVector2Display.AxesOrigin = [0.0, 0.0, 0.0]
warpByVector2Display.CubeAxesXLabelFormat = '%-#6.3g'
warpByVector2Display.CubeAxesYLabelFormat = '%-#6.3g'
warpByVector2Display.CubeAxesZLabelFormat = '%-#6.3g'
warpByVector2Display.StickyAxes = 0
warpByVector2Display.CenterStickyAxes = 0
warpByVector2Display.SelectionCellLabelBold = 0
warpByVector2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
warpByVector2Display.SelectionCellLabelFontFamily = 'Arial'
warpByVector2Display.SelectionCellLabelFontSize = 18
warpByVector2Display.SelectionCellLabelItalic = 0
warpByVector2Display.SelectionCellLabelJustification = 'Left'
warpByVector2Display.SelectionCellLabelOpacity = 1.0
warpByVector2Display.SelectionCellLabelShadow = 0
warpByVector2Display.SelectionPointLabelBold = 0
warpByVector2Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
warpByVector2Display.SelectionPointLabelFontFamily = 'Arial'
warpByVector2Display.SelectionPointLabelFontSize = 18
warpByVector2Display.SelectionPointLabelItalic = 0
warpByVector2Display.SelectionPointLabelJustification = 'Left'
warpByVector2Display.SelectionPointLabelOpacity = 1.0
warpByVector2Display.SelectionPointLabelShadow = 0
warpByVector2Display.ScalarOpacityUnitDistance = 0.3346560459821775
warpByVector2Display.SelectMapper = 'Projected tetra'

# set scalar coloring
ColorBy(warpByVector2Display, ('FIELD', 'vtkBlockColors'))

# change representation type
warpByVector2Display.SetRepresentationType('Surface With Edges')

# show color bar/color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# show color bar/color legend
warpByVector2Display.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).