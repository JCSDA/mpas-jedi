import argparse
import predefined_configs as pconf
import binning_utils as bu
import config as conf
import datetime as dt
import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from copy import deepcopy
from datetime import datetime, timedelta
import modelsp_utils as mu
import stat_utils as su
import var_utils as vu

def write_diag_stats():
  # Parse command line for date
  ap = argparse.ArgumentParser()
  ap.add_argument('date', default=os.getcwd().split('/')[-3], type=str, nargs = '?',
                    help='Valid date (%Y%m%d%H)')
  ap.add_argument('-r', '--referenceAnalysis', default = mu.GFSANA_DIR+'/x1.40962.init', type = str,
                        help='Path/prefix of reference analysis state')
  ap.add_argument('-m', '--mpasState', default = '../restart', type = str,
                        help='Path/prefix of arbitrary MPAS state')

  args = ap.parse_args()
  date = str(args.date)
  initDate = datetime.strptime(date,'%Y%m%d%H')
  fileDate= initDate.strftime('%Y-%m-%d_%H.%M.%S')

  #TODO: allow for entire ReferenceAnalysis and MPAS as args instead of date
  ReferenceAnalysis = str(args.referenceAnalysis)+'.'+fileDate+'.nc'
  MPAS = str(args.mpasState)+'.'+fileDate+'.nc'

  lats, lons = mu.readGrid(gridFile=ReferenceAnalysis)

  # Initialize a dictionary to contain all statistical info for this osKey
  statsDict = {}
  for attribName in su.fileStatAttributes:
    statsDict[attribName] = []
  for statName in su.allFileStats:
    statsDict[statName] = []

  dsKey = 'mpas'
  DiagSpaceGrp = conf.DiagSpaceConfig[dsKey]['DiagSpaceGrp']

  for varName in vu.modVarNames2d+vu.modVarNames3d:
    print('Working on '+varName)
    varShort, varUnits = vu.modelVarAttributes(varName)

    field = np.asarray(mu.varDiff(varName, ReferenceAnalysis, MPAS))
    diagName = 'mmgfsan'

    dims = field.shape
    nLevels = 1
    nDims = len(dims)
    if nDims == 2: nLevels = dims[1]

# TODO(JJG): generalize binning and diagnostics similar to DiagnoseObsStatistics
# TODO(JJG): extend to ACC diagnostic, 500mb geopotential height variable

    # no binVar/binMethod (all data)
    binMethod = bu.identityBinMethod
    binVar = vu.noBinVar
    binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
    binVal = vu.noBinVar
    binnedDiag = deepcopy(field).flatten()
    statsVal = su.calcStats(binnedDiag)
    for statName in su.allFileStats:
        statsDict[statName].append(statsVal[statName])
    statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
    statsDict['varName'].append(varShort)
    statsDict['varUnits'].append(varUnits)
    statsDict['diagName'].append(diagName)
    statsDict['binMethod'].append(binMethod)
    statsDict['binVar'].append(binVarShort)
    statsDict['binUnits'].append(binVarUnits)
    statsDict['binVal'].append(binVal)

    # binVar == vu.modVarLat (all levels)
    ## binMethods:
    ## (1) named latitude bands (str)
    ## (2) equidistant latitude ranges (float)
    binVar = vu.modVarLat
    binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
    for (binMethod, latitudeBins) in list(zip(
        [bu.latbandsMethod, bu.identityBinMethod],
        [pconf.namedLatBands, pconf.binLims[vu.modVarLat]],
        )):

      for (binVal, minBound, maxBound) in list(zip(
          latitudeBins['values'],
          latitudeBins['minBounds'],
          latitudeBins['maxBounds'],
          )):

        binnedDiagField = deepcopy(field)
        if nDims == 1:
          binnedDiagField[lats < minBound] = np.NaN
          binnedDiagField[lats > maxBound] = np.NaN
        elif nDims == 2:
          binnedDiagField[lats < minBound,:] = np.NaN
          binnedDiagField[lats > maxBound,:] = np.NaN

        binnedDiag = binnedDiagField.flatten()
        statsVal = su.calcStats(binnedDiag)
        for statName in su.allFileStats:
            statsDict[statName].append(statsVal[statName])
        statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
        statsDict['varName'].append(varShort)
        statsDict['varUnits'].append(varUnits)
        statsDict['diagName'].append(diagName)
        statsDict['binMethod'].append(binMethod)
        statsDict['binVar'].append(binVarShort)
        statsDict['binUnits'].append(binVarUnits)
        statsDict['binVal'].append(binVal)

    # GEO IR instrument lat/lon box binMethod, Region binVar
    binMethod = bu.geoirlatlonboxMethod
    binVar = vu.modelRegionBinVar
    binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
    for (binVal, minLon, maxLon, minLat, maxLat) in list(zip(
        pconf.geoirLonBands['values'],
        pconf.geoirLonBands['minBounds'],
        pconf.geoirLonBands['maxBounds'],
        pconf.geoirLatBands['minBounds'],
        pconf.geoirLatBands['maxBounds'],
        )):

      binnedDiagField = deepcopy(field)
      if nDims == 1:
        binnedDiagField[lons < minLon] = np.NaN
        binnedDiagField[lons > maxLon] = np.NaN
        binnedDiagField[lats < minLat] = np.NaN
        binnedDiagField[lats > maxLat] = np.NaN
      elif nDims == 2:
        binnedDiagField[lons < minLon,:] = np.NaN
        binnedDiagField[lons > maxLon,:] = np.NaN
        binnedDiagField[lats < minLat,:] = np.NaN
        binnedDiagField[lats > maxLat,:] = np.NaN

      binnedDiag = binnedDiagField.flatten()
      statsVal = su.calcStats(binnedDiag)
      for statName in su.allFileStats:
          statsDict[statName].append(statsVal[statName])
      statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
      statsDict['varName'].append(varShort)
      statsDict['varUnits'].append(varUnits)
      statsDict['diagName'].append(diagName)
      statsDict['binMethod'].append(binMethod)
      statsDict['binVar'].append(binVarShort)
      statsDict['binUnits'].append(binVarUnits)
      statsDict['binVal'].append(binVal)

    # binVar == vu.modVarLev
    if nDims > 1:
      binVar = vu.modVarLev
      binVarShort, binVarUnits = vu.modelVarAttributes(binVar)

      ## identity binMethod
      binMethod = bu.identityBinMethod
      binnedDiagField = deepcopy(field)

      for ilev in list(range(0,nLevels)):
        binnedDiag = binnedDiagField[:,ilev]
        statsVal = su.calcStats(binnedDiag)
        for statName in su.allFileStats:
            statsDict[statName].append(statsVal[statName])
        statsDict['binVal'].append(str(ilev))
      statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevels
      statsDict['varName'] += [varShort]*nLevels
      statsDict['varUnits'] += [varUnits]*nLevels
      statsDict['diagName'] += [diagName]*nLevels
      statsDict['binMethod'] += [binMethod]*nLevels
      statsDict['binVar'] += [binVarShort]*nLevels
      statsDict['binUnits'] += [binVarUnits]*nLevels

      ## named latitude band binMethods
      latitudeBins = pconf.namedLatBands
      for (latBand, minBound, maxBound) in list(zip(
          latitudeBins['values'],
          latitudeBins['minBounds'],
          latitudeBins['maxBounds'],
          )):

        binMethod = latBand

        binnedDiagField = deepcopy(field)
        binnedDiagField[lats < minBound,:] = np.NaN
        binnedDiagField[lats > maxBound,:] = np.NaN

        for ilev in list(range(0,nLevels)):
          binnedDiag = binnedDiagField[:,ilev]
          statsVal = su.calcStats(binnedDiag)
          for statName in su.allFileStats:
              statsDict[statName].append(statsVal[statName])
          statsDict['binVal'].append(str(ilev))
        statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevels
        statsDict['varName'] += [varShort]*nLevels
        statsDict['varUnits'] += [varUnits]*nLevels
        statsDict['diagName'] += [diagName]*nLevels
        statsDict['binMethod'] += [binMethod]*nLevels
        statsDict['binVar'] += [binVarShort]*nLevels
        statsDict['binUnits'] += [binVarUnits]*nLevels

      ## GEO IR instrument lat/lon box binMethod
      for (geoirInst, minLon, maxLon, minLat, maxLat) in list(zip(
          pconf.geoirLonBands['values'],
          pconf.geoirLonBands['minBounds'],
          pconf.geoirLonBands['maxBounds'],
          pconf.geoirLatBands['minBounds'],
          pconf.geoirLatBands['maxBounds'],
          )):

        binMethod = geoirInst

        binnedDiagField = deepcopy(field)
        binnedDiagField[lons < minLon,:] = np.NaN
        binnedDiagField[lons > maxLon,:] = np.NaN
        binnedDiagField[lats < minLat,:] = np.NaN
        binnedDiagField[lats > maxLat,:] = np.NaN

        for ilev in list(range(0,nLevels)):
          binnedDiag = binnedDiagField[:,ilev]
          statsVal = su.calcStats(binnedDiag)
          for statName in su.allFileStats:
              statsDict[statName].append(statsVal[statName])
          statsDict['binVal'].append(str(ilev))
        statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevels
        statsDict['varName'] += [varShort]*nLevels
        statsDict['varUnits'] += [varUnits]*nLevels
        statsDict['diagName'] += [diagName]*nLevels
        statsDict['binMethod'] += [binMethod]*nLevels
        statsDict['binVar'] += [binVarShort]*nLevels
        statsDict['binUnits'] += [binVarUnits]*nLevels

  su.write_stats_nc(dsKey, statsDict)

def main():
  write_diag_stats()

if __name__ == '__main__': main()
