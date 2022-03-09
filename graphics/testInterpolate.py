#!/usr/bin/env python3

from basic_plot_functions import scatterMapFields
from copy import deepcopy
import Interpolate
import modelsp_utils as modelUtils
import numpy as np

initDate = '2018041500'
mapCentralLongitude = 0.

imageFileExtension = 'png'

# choose the test variables/level
testVariables = ['theta', 'qv', 'uReconstructZonal']
testLevel = 0

#referenceCase = 'barycentricScalar'
referenceCase = 'barycentric'
weightCases = {
#  'barycentricScalar': {
#    'weightMethod': 'barycentricScalar',
#    'nInterpPoints': 5,
#  },
  'barycentric': {
    'weightMethod': 'barycentric',
    'nInterpPoints': 6,
  },
  'pbarycentric': {
    'weightMethod': 'pbarycentric',
    'nInterpPoints': 4,
  },
#  'unstinterpScalar3': {
#    'weightMethod': 'unstinterpScalar',
#    'nInterpPoints': 3,
#  },
  'unstinterp3': {
    'weightMethod': 'unstinterp',
    'nInterpPoints': 3,
  },
#  'inverseD3': {
#    'weightMethod': 'inverseD',
#    'nInterpPoints': 3,
#  },
}

meshConfigs = {
  'coarse':{
    'directory': '/glade/p/mmm/parc/liuz/pandac_common/120km_GFSANA',
    'nCells': 40962,
  },
  'fine':{
    'directory': '/glade/p/mmm/parc/liuz/pandac_common/30km_GFSANA',
    'nCells': 655362,
  },
}

gridFiles = {}
meshLatDeg = {}
meshLonDeg = {}
meshLat = {}
meshLon = {}
for meshName, conf in meshConfigs.items():
  gridFiles[meshName] = modelUtils.getGridFile(initDate,
    conf['directory'], conf['nCells'],
  )
  grid = modelUtils.readGrid(gridFile=gridFiles[meshName])
  meshLatDeg[meshName] = deepcopy(grid['latitude'])
  meshLonDeg[meshName] = deepcopy(grid['longitude'])
  del grid
  meshLat[meshName] = meshLatDeg[meshName] * np.pi / 180.0
  meshLon[meshName] = meshLonDeg[meshName] * np.pi / 180.0

lons = {
  'coarse': meshLonDeg['coarse'],
  'fine': meshLonDeg['fine'],
}
lats = {
  'coarse': meshLatDeg['coarse'],
  'fine': meshLatDeg['fine'],
}
markers = {
  'coarse': '.',
  'fine': '.',
}

sizes = {
  'coarse': 1.5,
  'fine': 0.5,
}

minAbsNormError = {}
maxAbsNormError = {}
meanAbsNormError = {}

minAbsError = {}
maxAbsError = {}
meanAbsError = {}

minBias = {}
maxBias = {}
meanBias = {}

meshField = {}
for testVariable in testVariables:
  meshField[testVariable] = {}
  for meshName, conf in meshConfigs.items():
    meshField[testVariable][meshName] = \
      modelUtils.varRead(testVariable, gridFiles[meshName])[:,testLevel]

  fields = {
    'coarse': meshField[testVariable]['coarse']
  }

  print('Plotting coarse field: '+testVariable)
  scatterMapFields(lons, lats, fields, 'original_coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                   cLon = mapCentralLongitude,
                   markers = markers, sizes = sizes,
                   projection = 'moll',
  )

  minAbsNormError[testVariable] = np.NaN
  maxAbsNormError[testVariable] = np.NaN

  minAbsError[testVariable] = np.NaN
  maxAbsError[testVariable] = np.NaN

  minBias[testVariable] = np.NaN
  maxBias[testVariable] = np.NaN

## generate the double-interpolated field and calculate errors
absNormErrorFields = {}
absErrorFields = {}
biasFields = {}
for caseName, case in weightCases.items():
  weightMethod = case['weightMethod']
  print('Testing caseName = '+caseName)
  print('generating weights for interpolating from coarse to fine mesh')
  coarse2fine = Interpolate.InterpolateLonLat(meshLon['coarse'], meshLat['coarse'],
    nInterpPoints = case['nInterpPoints'],
    weightMethod = weightMethod,
    calculateDiagnostics = True)
  coarse2fine.initWeights(meshLon['fine'], meshLat['fine'])

  print('generating weights for interpolating from fine to coarse mesh')
  fine2coarse = Interpolate.InterpolateLonLat(meshLon['fine'], meshLat['fine'],
    nInterpPoints = case['nInterpPoints'],
    weightMethod = weightMethod,
    calculateDiagnostics = True)
  fine2coarse.initWeights(meshLon['coarse'], meshLat['coarse'])

  if 'barycentric' in caseName:
    print('plotting barycentric combination used, triangle area, and determinant')
    fields['coarse'] = fine2coarse.combinationUsed
    scatterMapFields(lons, lats, fields, caseName+'_combinationUsedCoarse.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     cmap = 'tab10',
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['coarse'] = fine2coarse.triangleArea
    scatterMapFields(lons, lats, fields, caseName+'_triangleAreaCoarse.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     cmap = 'gist_rainbow',
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['coarse'] = fine2coarse.determinant
    scatterMapFields(lons, lats, fields, caseName+'_determinantCoarse.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['coarse'] = fine2coarse.ycoord
    scatterMapFields(lons, lats, fields, caseName+'_ycoordCoarse.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes, cbarType = 'Log',
                     projection = 'moll',
    )
    fields['coarse'] = fine2coarse.aspectRatio
    scatterMapFields(lons, lats, fields, caseName+'_aspectRatioCoarse.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    del fields['coarse']

    fields['fine'] = coarse2fine.combinationUsed
    scatterMapFields(lons, lats, fields, caseName+'_combinationUsedFine.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     cmap = 'tab10',
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['fine'] = coarse2fine.triangleArea
    scatterMapFields(lons, lats, fields, caseName+'_triangleAreaFine.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     cmap = 'gist_rainbow',
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['fine'] = coarse2fine.determinant
    scatterMapFields(lons, lats, fields, caseName+'_determinantFine.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    fields['fine'] = coarse2fine.ycoord
    scatterMapFields(lons, lats, fields, caseName+'_ycoordFine.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes, cbarType = 'Log',
                     projection = 'moll',
    )
    fields['fine'] = coarse2fine.aspectRatio
    scatterMapFields(lons, lats, fields, caseName+'_aspectRatioFine.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     markers = markers, sizes = sizes,
                     projection = 'moll',
    )
    del fields['fine']

  absNormErrorFields[caseName] = {}
  absErrorFields[caseName] = {}
  biasFields[caseName] = {}

  meanAbsNormError[caseName] = {}
  meanAbsError[caseName] = {}
  meanBias[caseName] = {}


  for testVariable in testVariables:
    print('Interpolating and calculating errors for testVariable = '+testVariable)
    coarse2fineField = coarse2fine.apply(meshField[testVariable]['coarse'])
    coarse2fine2coarseField = fine2coarse.apply(coarse2fineField)

#    print('plotting double-interpolated field')
#    fields['coarse'] = coarse2fine2coarseField
#    scatterMapFields(lons, lats, fields, caseName+'_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
#                     cLon = mapCentralLongitude,
#                     markers = markers, sizes = sizes,
#                     projection = 'moll',
#    )

    absNormErrorFields[caseName][testVariable] = \
      np.divide(np.abs(coarse2fine2coarseField - meshField[testVariable]['coarse']), np.abs(meshField[testVariable]['coarse']))
    minAbsNormError_ = np.min(absNormErrorFields[caseName][testVariable])
    maxAbsNormError_ = np.max(absNormErrorFields[caseName][testVariable])
    print(minAbsNormError_, maxAbsNormError_)
    minAbsNormError[testVariable] = np.nanmin([minAbsNormError_, minAbsNormError[testVariable]])
    maxAbsNormError[testVariable] = np.nanmax([maxAbsNormError_, maxAbsNormError[testVariable]])

    meanAbsNormError[caseName][testVariable] = np.nanmean(absNormErrorFields[caseName][testVariable])

    absErrorFields[caseName][testVariable] = \
      np.abs(coarse2fine2coarseField - meshField[testVariable]['coarse'])
    minAbsError_ = np.min(absErrorFields[caseName][testVariable])
    maxAbsError_ = np.max(absErrorFields[caseName][testVariable])
    print(minAbsError_, maxAbsError_)
    minAbsError[testVariable] = np.nanmin([minAbsError_, minAbsError[testVariable]])
    maxAbsError[testVariable] = np.nanmax([maxAbsError_, maxAbsError[testVariable]])

    meanAbsError[caseName][testVariable] = np.nanmean(absErrorFields[caseName][testVariable])

    biasFields[caseName][testVariable] = \
      coarse2fine2coarseField - meshField[testVariable]['coarse']
    minBias_ = np.min(biasFields[caseName][testVariable])
    maxBias_ = np.max(biasFields[caseName][testVariable])
    print(minBias_, maxBias_)
    minBias[testVariable] = np.nanmin([minBias_, minBias[testVariable]])
    maxBias[testVariable] = np.nanmax([maxBias_, maxBias[testVariable]])

    meanBias[caseName][testVariable] = np.nanmean(biasFields[caseName][testVariable])

print('meanAbsNormError = ', meanAbsNormError)
print('meanAbsError = ', meanAbsError)
print('meanBias = ', meanBias)

## plot the errors
for caseName, case in weightCases.items():
  for testVariable in testVariables:
    print('Plotting error metrics for caseName = '+caseName+' and testVariable = '+testVariable)

    weightMethod = case['weightMethod']
    fields['coarse'] = absNormErrorFields[caseName][testVariable]
    scatterMapFields(lons, lats, fields, caseName+'_NormalizedAbsError_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                     cLon = mapCentralLongitude,
  #                   dmin = minAbsNormError[testVariable], dmax = maxAbsNormError[testVariable],
                     dmin = 1.e-8, dmax = 1.0,
                     cmap = 'gist_rainbow',
                     markers = markers, sizes = sizes, cbarType = 'Log',
                     projection = 'moll',
    )

    fields['coarse'] = absErrorFields[caseName][testVariable]
    scatterMapFields(lons, lats, fields, caseName+'_AbsError_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     dmin = minAbsError[testVariable], dmax = maxAbsError[testVariable],
                     cmap = 'gist_rainbow',
                     markers = markers, sizes = sizes, cbarType = 'Log',
                     projection = 'moll',
    )

    dd = np.max([minBias[testVariable], maxBias[testVariable]])
    fields['coarse'] = biasFields[caseName][testVariable]
    scatterMapFields(lons, lats, fields, caseName+'_Bias_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                     cLon = mapCentralLongitude,
                     dmin = -dd, dmax = dd,
                     cmap = 'BrBG',
                     markers = markers, sizes = sizes, cbarType = 'SymLog',
                     projection = 'moll',
    )

    if caseName != referenceCase:
      fields['coarse'] = absNormErrorFields[caseName][testVariable] - absNormErrorFields[referenceCase][testVariable]
      scatterMapFields(lons, lats, fields, caseName+'-'+referenceCase+'_NormalizedAbsError_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                       cLon = mapCentralLongitude,
                       dmin = -0.1, dmax = 0.1,
                       cmap = 'BrBG',
                       markers = markers, sizes = sizes, cbarType = 'SymLog',
                       projection = 'moll',
      )

      referenceAbsError = deepcopy(absErrorFields[referenceCase][testVariable])
      referenceAbsError[referenceAbsError == 0.0] = 1.e-16*meshField[testVariable]['coarse'][referenceAbsError == 0.0]

      caseAbsError = deepcopy(absErrorFields[caseName][testVariable])
      caseAbsError[caseAbsError == 0.0] = 1.e-16*meshField[testVariable]['coarse'][caseAbsError == 0.0]

      fields['coarse'] = np.divide(caseAbsError, referenceAbsError)
      scatterMapFields(lons, lats, fields, caseName+'DIV'+referenceCase+'_AbsError_coarse2fine2coarse_'+testVariable+'_l'+str(testLevel)+'.'+imageFileExtension,
                       cLon = mapCentralLongitude,
                       dmin = 1.e-3, dmax = 1.e3,
                       cmap = 'BrBG',
                       markers = markers, sizes = sizes, cbarType = 'Log',
                       projection = 'moll',
      )
