import argparse
import binning_utils as bu
import config as conf
from copy import deepcopy
from datetime import datetime
import datetime as dt
import logging
import logsetup
import multiprocessing as mp
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import modelsp_utils as mu
import os
import predefined_configs as pconf
import stat_utils as su
import var_utils as vu

_logger = logging.getLogger(__name__)

def write_diag_stats():
  # Parse command line for date
  ap = argparse.ArgumentParser()
  ap.add_argument('date', default=os.getcwd().split('/')[-3], type=str, nargs = '?',
                  help='Valid date (YYYYMMDDHH)')
  ap.add_argument('-r', '--referenceAnalysis', default = mu.GFSANA_DIR+'/x1.'+str(mu.ncells)+'.init',
                  type = str,
                  help='Path/prefix of reference analysis state')
  ap.add_argument('-m', '--mpasState', default = '../restart', type = str,
                  help='Path/prefix of arbitrary MPAS state')
  ap.add_argument('-n', '--nprocs', default = 1, type = int,
                  help='Number of tasks/processors for multiprocessing')

  args = ap.parse_args()
  date = str(args.date)
  initDate = datetime.strptime(date,'%Y%m%d%H')
  fileDate= initDate.strftime('%Y-%m-%d_%H.%M.%S')

  #TODO: allow for entire ReferenceAnalysis and MPAS as args instead of date
  ReferenceFile = str(args.referenceAnalysis)+'.'+fileDate+'.nc'
  MPASFile = str(args.mpasState)+'.'+fileDate+'.nc'

  _logger.info('ReferenceAnalysis: '+ReferenceFile)
  _logger.info('MPAS: '+MPASFile)

  ReferenceAnalysis = mu.getNCData(ReferenceFile)
  MPAS = mu.getNCData(MPASFile)

  nprocs = min(mp.cpu_count(), args.nprocs)

  _logger.info('Reading mesh description')
  grid = mu.readGrid(gridFile=ReferenceFile)
  latitude = grid['latitude']
  longitude = grid['longitude']

  spaceKey = 'mpas'
  DiagSpaceGrp = conf.DiagSpaceConfig[spaceKey]['DiagSpaceGrp']

  if nprocs > 1:
    workers = mp.Pool(processes = nprocs)
  else:
    workers = None


  subStats = []
  for varName in vu.modVarNames2d+vu.modVarNames3d:

    for diagName in mu.variableSpecificDiagnostics(varName):

      diagFunction = mu.diagnosticFunctions[diagName]
      modelVarName = mu.aggModelVariable(varName)
      diagnostic = np.asarray(diagFunction(modelVarName, ReferenceAnalysis, MPAS))

      if workers is None:
        subStats.append(singleVarDiagnostic(
          varName, diagName, DiagSpaceGrp,
          latitude, longitude,
          diagnostic))
      else:
        subStats.append(workers.apply_async(singleVarDiagnostic,
          args=(
            varName, diagName, DiagSpaceGrp,
            latitude, longitude,
            diagnostic,
          )
        ))

  # Fill in the global statsDict for all variables and diagnostics
  statsDict = {}
  for aa in su.fileStatAttributes:
    statsDict[aa] = []
  for ss in su.allFileStats:
    statsDict[ss] = []

  if workers is None:
    for stats in subStats:
      goodValues = True
      for name, values in stats.items():
        if len(values) < 0: goodValues = False
      if goodValues:
        for name, values in stats.items():
          statsDict[name] += values
  else:
    workers.close()
    workers.join()
    for stats in subStats:
      goodValues = True
      for name, values in stats.get().items():
        if len(values) < 0: goodValues = False
      if goodValues:
        for name, values in stats.get().items():
          statsDict[name] += values
      else:
        self.logger.info('badValues')

#  for name, values in statsDict.items():
#    _logger.info(name+': '+str(len(values)))

  _logger.info('Writing statistics file')

  stats = su.BinnedStatisticsFile(statSpace=spaceKey)
  stats.write(statsDict)

  _logger.info('Finished')


def singleVarDiagnostic(
  varName, diagName, DiagSpaceGrp,
  latitude, longitude,
  diagnostic_):

# TODO(JJG): generalize binning and diagnostics similar to DiagnoseObsStatistics
# TODO(JJG): extend to ACC diagnostic, 500mb geopotential height variable

  logSuffix = ' for ('+varName+', '+diagName+'):'
  _logger.info('Starting singleVarDiagnostic'+logSuffix)

  statsDict = {}
  for aa in su.fileStatAttributes:
    statsDict[aa] = []
  for ss in su.allFileStats:
    statsDict[ss] = []

  varShort, varUnits = vu.modelVarAttributes(varName)

  diagnostic = deepcopy(diagnostic_)

  dShape = diagnostic.shape
  nDims = len(dShape)
  nLevels = 1
  minLevel = 1
  maxLevel = 1
  if nDims == 2:
    nLevels = dShape[1]
    minLevel = mu.aggMinLevel(varName)
    maxLevel = mu.aggMaxLevel(varName, nLevels)

    # mask levels > maxLevel, where index of first level is 1
    if maxLevel < nLevels: diagnostic[:, maxLevel:] = np.NaN

    # mask levels < minLevel, where index of first level is 1
    if minLevel > 1: diagnostic[:, :(minLevel-1)] = np.NaN


  # no binVar/binMethod (all data)
  _logger.info(' all-cell aggregation'+logSuffix)
  binMethod = bu.identityBinMethod
  binVar = vu.noBinVar
  binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
  binVal = vu.noBinVar
  binnedDiag = deepcopy(diagnostic).flatten()
  statsVal = su.calcStats(binnedDiag)
  for ss in su.allFileStats:
      statsDict[ss].append(statsVal[ss])
  statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
  statsDict['varName'].append(varShort)
  statsDict['varUnits'].append(varUnits)
  statsDict['diagName'].append(diagName)
  statsDict['binMethod'].append(binMethod)
  statsDict['binVar'].append(binVarShort)
  statsDict['binUnits'].append(binVarUnits)
  statsDict['binVal'].append(binVal)


  # binVar == vu.modVarLat (all levels)
  _logger.info(' binned by latitude'+logSuffix)
  ## binMethods:
  ## (1) named tropical latitude bands (str)
  ## (2) equidistant latitude ranges (float)
  binVar = vu.modVarLat

  latLims1D = {}
  latLims1D['starts'] = pconf.binAxes1D[binVar].starts()
  latLims1D['stops'] = pconf.binAxes1D[binVar].stops()
  latLims1D['values'] = pconf.binAxes1D[binVar].values()

  binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
  for (binMethod, latitudeBins) in list(zip(
      [bu.troplatbandsMethod, bu.identityBinMethod],
      [pconf.namedTropLatBands, latLims1D],
      )):

    for (binVal, start, stop) in list(zip(
        latitudeBins['values'],
        latitudeBins['starts'],
        latitudeBins['stops'],
        )):

      binnedDiagField = deepcopy(diagnostic)

      if nDims == 1:
        binnedDiagField[latitude < start] = np.NaN
        binnedDiagField[latitude >= stop] = np.NaN
      elif nDims == 2:
        binnedDiagField[latitude < start, :] = np.NaN
        binnedDiagField[latitude >= stop, :] = np.NaN

      binnedDiag = binnedDiagField.flatten()
      statsVal = su.calcStats(binnedDiag)
      for ss in su.allFileStats:
          statsDict[ss].append(statsVal[ss])
      statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
      statsDict['varName'].append(varShort)
      statsDict['varUnits'].append(varUnits)
      statsDict['diagName'].append(diagName)
      statsDict['binMethod'].append(binMethod)
      statsDict['binVar'].append(binVarShort)
      statsDict['binUnits'].append(binVarUnits)
      statsDict['binVal'].append(binVal)


  # GEO IR instrument lat/lon box binMethod, Region binVar
  _logger.info(' binned by GEOIR footprints'+logSuffix)
  binMethod = bu.geoirlatlonboxMethod
  binVar = vu.modelRegionBinVar
  binVarShort, binVarUnits = vu.modelVarAttributes(binVar)
  for (binVal, minLon, maxLon, minLat, maxLat) in list(zip(
      pconf.geoirLonBands['values'],
      pconf.geoirLonBands['starts'],
      pconf.geoirLonBands['stops'],
      pconf.geoirLatBands['starts'],
      pconf.geoirLatBands['stops'],
      )):

    binnedDiagField = deepcopy(diagnostic)
    if nDims == 1:
      binnedDiagField[longitude < minLon] = np.NaN
      binnedDiagField[longitude > maxLon] = np.NaN
      binnedDiagField[latitude < minLat] = np.NaN
      binnedDiagField[latitude > maxLat] = np.NaN
    elif nDims == 2:
      binnedDiagField[longitude < minLon, :] = np.NaN
      binnedDiagField[longitude > maxLon, :] = np.NaN
      binnedDiagField[latitude < minLat, :] = np.NaN
      binnedDiagField[latitude > maxLat, :] = np.NaN

    binnedDiag = binnedDiagField.flatten()
    statsVal = su.calcStats(binnedDiag)
    for ss in su.allFileStats:
        statsDict[ss].append(statsVal[ss])
    statsDict['DiagSpaceGrp'].append(DiagSpaceGrp)
    statsDict['varName'].append(varShort)
    statsDict['varUnits'].append(varUnits)
    statsDict['diagName'].append(diagName)
    statsDict['binMethod'].append(binMethod)
    statsDict['binVar'].append(binVarShort)
    statsDict['binUnits'].append(binVarUnits)
    statsDict['binVal'].append(binVal)


  # level-resolved bins for level-dependent fields only
  if nDims == 2 and minLevel==1 and maxLevel==nLevels:
    _logger.info(' binned by model level'+logSuffix)

    allLevels = np.arange(1, nLevels+1)
    nLevelBins = len(allLevels)

    binVar = vu.modVarLev
    binVarShort, binVarUnits = vu.modelVarAttributes(binVar)


    ## identity binMethod
    _logger.info('   + all-cells'+logSuffix)

    binMethod = bu.identityBinMethod
    binnedDiagField = deepcopy(diagnostic)

    for ilev in allLevels:
      binnedDiag = binnedDiagField[:, ilev-1].flatten()
      statsVal = su.calcStats(binnedDiag)
      for ss in su.allFileStats:
          statsDict[ss].append(statsVal[ss])
      statsDict['binVal'].append(str(ilev))
    statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevelBins
    statsDict['varName'] += [varShort]*nLevelBins
    statsDict['varUnits'] += [varUnits]*nLevelBins
    statsDict['diagName'] += [diagName]*nLevelBins
    statsDict['binMethod'] += [binMethod]*nLevelBins
    statsDict['binVar'] += [binVarShort]*nLevelBins
    statsDict['binUnits'] += [binVarUnits]*nLevelBins


    ## named latitude band binMethods
    _logger.info('   + latitude'+logSuffix)
    latitudeBins = pconf.namedTropLatBands
    for (latBand, start, stop) in list(zip(
        latitudeBins['values'],
        latitudeBins['starts'],
        latitudeBins['stops'],
        )):

      binMethod = latBand

      binnedDiagField = deepcopy(diagnostic)
      binnedDiagField[latitude < start, :] = np.NaN
      binnedDiagField[latitude >= stop, :] = np.NaN

      for ilev in allLevels:
        binnedDiag = binnedDiagField[:, ilev-1].flatten()
        statsVal = su.calcStats(binnedDiag)
        for ss in su.allFileStats:
            statsDict[ss].append(statsVal[ss])
        statsDict['binVal'].append(str(ilev))
      statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevelBins
      statsDict['varName'] += [varShort]*nLevelBins
      statsDict['varUnits'] += [varUnits]*nLevelBins
      statsDict['diagName'] += [diagName]*nLevelBins
      statsDict['binMethod'] += [binMethod]*nLevelBins
      statsDict['binVar'] += [binVarShort]*nLevelBins
      statsDict['binUnits'] += [binVarUnits]*nLevelBins


    ## GEO IR instrument lat/lon box binMethod
    _logger.info('   + GEOIR footprints'+logSuffix)
    for (geoirInst, minLon, maxLon, minLat, maxLat) in list(zip(
        pconf.geoirLonBands['values'],
        pconf.geoirLonBands['starts'],
        pconf.geoirLonBands['stops'],
        pconf.geoirLatBands['starts'],
        pconf.geoirLatBands['stops'],
        )):

      binMethod = geoirInst

      binnedDiagField = deepcopy(diagnostic)
      binnedDiagField[longitude < minLon, :] = np.NaN
      binnedDiagField[longitude > maxLon, :] = np.NaN
      binnedDiagField[latitude < minLat, :] = np.NaN
      binnedDiagField[latitude > maxLat, :] = np.NaN

      for ilev in allLevels:
        binnedDiag = binnedDiagField[:, ilev-1].flatten()
        statsVal = su.calcStats(binnedDiag)
        for ss in su.allFileStats:
            statsDict[ss].append(statsVal[ss])
        statsDict['binVal'].append(str(ilev))
      statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nLevelBins
      statsDict['varName'] += [varShort]*nLevelBins
      statsDict['varUnits'] += [varUnits]*nLevelBins
      statsDict['diagName'] += [diagName]*nLevelBins
      statsDict['binMethod'] += [binMethod]*nLevelBins
      statsDict['binVar'] += [binVarShort]*nLevelBins
      statsDict['binUnits'] += [binVarUnits]*nLevelBins

  return statsDict

def main():
  write_diag_stats()
  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
