import DiagnoseModelStatisticsArgs

import binning_utils as bu
import config as conf
from copy import deepcopy
from datetime import datetime
import datetime as dt
import logging
import logsetup
import multiprocessing as mp
import numpy as np
import modelsp_utils as mu
import os
import predefined_configs as pconf
import stat_utils as su
import var_utils as vu

_logger = logging.getLogger(__name__)

class DiagnoseModelStatistics():
  '''
  Diagnose model-space statistics
  Driven by
    - static selections in conf
    - command-line arguments in DiagnoseModelStatisticsArgs
  '''
  def __init__(self):
    self.name = 'DiagnoseModelStatistics'
    self.args = DiagnoseModelStatisticsArgs.args
    self.logger = logging.getLogger(self.name)
    self.nprocs = min(mp.cpu_count(), self.args.nprocs)

    # construct mean DB into 0th member slot
    date = str(self.args.date)
    initDate = datetime.strptime(date,'%Y%m%d%H')
    fileDate= initDate.strftime('%Y-%m-%d_%H.%M.%S')

    self.files = {'state': {}, 'diagnostics':{}}

    # mean/deterministic files['state']
    #TODO: allow for entire referenceState and meanState as args instead of date
    self.files['state']['reference'] = str(self.args.referenceState)+'.'+fileDate+'.nc'
    self.files['state']['mean'] = str(self.args.meanState)+'.'+fileDate+'.nc'

    self.logger.info('ReferenceState: '+self.files['state']['reference'])
    self.logger.info('MeanState: '+self.files['state']['mean'])

    # mean/deterministic files['diagnostics'] (optional)
    file = str(self.args.referenceDiagnostics)+'.'+fileDate+'.nc'
    if os.path.exists(file):
      self.files['diagnostics']['reference'] = file
    else:
      self.files['diagnostics']['reference'] = None

    file = str(self.args.meanDiagnostics)+'.'+fileDate+'.nc'
    if os.path.exists(file):
      self.files['diagnostics']['mean'] = file
    else:
      self.files['diagnostics']['mean'] = None

    self.logger.info('ReferenceDiagnostics: '+str(self.files['diagnostics']['reference']))
    self.logger.info('MeanDiagnostics: '+str(self.files['diagnostics']['mean']))


    # ensemble member files['state']
    # background
    self.files['state']['bgEns'] = []
    ensemblePath = (str(self.args.bgEnsemblePath)+'.'+fileDate+'.nc').replace('YYYYMMDDHH', date)
    for member in list(range(1, self.args.nMembers+1)):
      file = ensemblePath.format(member)
      if os.path.exists(file):
        self.logger.info('adding background member file: '+file)
        self.files['state']['bgEns'].append(file)

    # analysis
    self.files['state']['anEns'] = []
    ensemblePath = (str(self.args.anEnsemblePath)+'.'+fileDate+'.nc').replace('YYYYMMDDHH', date)
    for member in list(range(1, self.args.nMembers+1)):
      file = ensemblePath.format(member)
      if os.path.exists(file):
        self.logger.info('adding analysis member file: '+file)
        self.files['state']['anEns'].append(file)

    # inflated
    self.files['state']['infEns'] = []
    ensemblePath = (str(self.args.inflatedEnsemblePath)+'.'+fileDate+'.nc').replace('YYYYMMDDHH', date)
    for member in list(range(1, self.args.nMembers+1)):
      file = ensemblePath.format(member)
      if os.path.exists(file):
        self.logger.info('adding inflated member file: '+file)
        self.files['state']['infEns'].append(file)

    # account for when there is no inflation
    if len(self.files['state']['anEns']) == 0 and self.args.nMembers > 1:
      self.logger.info('no valid analysis member files; assuming inflated members are identical to analysis members')
      self.files['state']['anEns'] = deepcopy(self.files['state']['infEns'])
      self.files['state']['infEns'] = []

    # forecast
#    self.files['state']['fcEns'] = []
#    ensemblePath = (str(self.args.fcEnsemblePath)+'.'+fileDate+'.nc').replace('YYYYMMDDHH', date)
#    for member in list(range(1, self.args.nMembers+1)):
#      file = ensemblePath.format(member)
#      self.logger.info('adding member file: '+file)
#      self.files['state']['fcEns'].append(file)

  def diagnose(self):

    self.logger.info('Populating binning metadata')

    # initialize netcdf4 Dataset's
    dataSets = {}
    for ds in ['reference', 'mean']:
      dataSets[ds] = {}
      dataSets[ds]['state'] = mu.getNCData(self.files['state'][ds])
      if self.files['diagnostics'][ds] is None:
        dataSets[ds]['diagnostics'] = None
      else:
        dataSets[ds]['diagnostics'] = mu.getNCData(self.files['diagnostics'][ds])

    for ds in ['bgEns', 'anEns', 'infEns']: # , 'fcEns']
      dataSets[ds] = []
      #for stateFile, diagnosticFile in zip(self.files['state'][ds], self.files['diagnostics'][ds]):
      for stateFile in self.files['state'][ds]:
        dataSets[ds].append({
          'state': mu.getNCData(stateFile),
          'diagnostics': None,
          #'diagnostics': mu.getNCData(diagnosticFile),
        })

    # read "metadata" fields that are used for binning and
    #  create reshaped copies to fit 1D and 2D arrays
    pressure = mu.varRead(vu.modVarPrs, dataSets['reference'])
    nCells = pressure.shape[0]
    nVertLevels = pressure.shape[1]
    nVertLevelsP1 = nVertLevels+1
    nDiagLevs = len(mu.diagnosticPressures)

    dbVals = {}
    allLevelSizes = [1, nVertLevels, nVertLevelsP1, nDiagLevs]
    for nn in allLevelSizes:
      dbVals[(nn,)] = {}

    # modVarLev
    # unique for variables with different numbers of levels, not defined for nDiagLevs
    for nn in allLevelSizes:
    #for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)][vu.modVarLev] = np.arange(1, nn+1)

    # modVarLat, modVarLon (same for all level-distributions)
    grid = mu.readGrid(gridFile=self.files['state']['reference'])
    for nn in allLevelSizes:
      dbVals[(nn,)][vu.modVarLat] = grid['latitude']
      dbVals[(nn,)][vu.modVarLon] = grid['longitude']

    # add pressure as a binning variable
    # TODO: enable binning by altitude
    dbVals[(nDiagLevs,)][vu.modVarDiagPrs] = mu.diagnosticPressures


    ###############################################
    ## Extract constructor info about the DiagSpace
    ###############################################

    DiagSpaceName = 'mpas'
    DiagSpaceInfo = conf.DiagSpaceConfig[DiagSpaceName]
    DiagSpaceGrp = DiagSpaceInfo['DiagSpaceGrp']
    binVarConfigs = DiagSpaceInfo.get('binVarConfigs',{})
    selectDiagNames = DiagSpaceInfo.get('diagNames',{})

    modelVars = vu.modVarNames2d+vu.modVarNames3d


    ########################################################
    ## Construct dictionary of binMethods for ModelSpace
    ########################################################

    self.logger.info('Initializing binMethods')

    binMethods = {}

    for binVarKey, binMethodNames in binVarConfigs.items():
        binVarConfig = pconf.binVarConfigs.get(binVarKey,pconf.nullBinVarConfig)
        for binMethodName in binMethodNames:
            config = binVarConfig.get(binMethodName,pconf.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['dsName'] = DiagSpaceName
            config['fileFormat'] = 'model'

            binMethods[(binVarKey, binMethodName)] = bu.BinMethod(config)

#    ######################################
#    ## Construct diagnostic configurations
#    ######################################
#
#    logger.info('Initializing diagnosticConfigs')
#
#    nMembers = 1
#
#    diagnosticConfigs = du.diagnosticConfigs(
#        selectDiagNames, ObsSpaceName,
#        includeEnsembleDiagnostics = (nMembers > 1),
#        fileFormat = 'model')

    ######################################
    ## Collect statistics for all modelVars
    ######################################

    self.logger.info('Calculating diagnostic statistics')
    self.logger.info("with "+str(self.nprocs)+" out of "+str(mp.cpu_count())+" processors")

    # Fill in the global statsDict for all variables and diagnostics
    statsDict = {}
    for aa in su.fileStatAttributes:
      statsDict[aa] = []
    for ss in su.allFileStats:
      statsDict[ss] = []

    if self.nprocs > 1:
      workers = mp.Pool(processes = self.nprocs)
    else:
      workers = None

    subStats = []

    for varName in modelVars: 

      self.logger.info('varName: '+varName)

      # TODO: first, get varName arrays from all dataSets that are needed by all diagnostics
      # then expensive quanities can be reused as needed

      for diagName in mu.variableSpecificDiagnostics(varName, len(dataSets['bgEns'])):
        # TODO(JJG): extend to ACC diagnostic, e.g., for 500mb geopotential height

        modelVarName = mu.aggModelVariable(varName)

        diagFunction = mu.diagnosticFunctions[diagName]
        diagnostic = diagFunction.evaluate(modelVarName, dataSets, self.nprocs)

        self.logger.info('diagnostic calculated: '+diagName)

        dShape = diagnostic.shape
        nDims = len(dShape)
        assert (nDims==1 or nDims==2), 'Number of dimensions must be 1 or 2'

        nn = 1
        if nDims==2:
          nn = dShape[1]
        dbValsNN = dbVals[(nn,)]

        # TODO: move this binning to a generic binMethod instead of creating
        #       unique level-binned variables for, e.g., qv

        if len(dShape)==2 and varName not in vu.modDiagnosticVarNames:
          minLevel = mu.aggMinLevel(varName)
          maxLevel = mu.aggMaxLevel(varName, nn)

          # mask levels < minLevel
          mask = bu.lessBound(dbValsNN[vu.modVarLev], minLevel)
          diagnostic[:,mask] = np.NaN

          # mask levels > maxLevel
          mask = bu.greatBound(dbValsNN[vu.modVarLev], maxLevel)
          diagnostic[:,mask] = np.NaN

        if np.isfinite(diagnostic).sum() == 0:
          self.logger.warning('All missing values for (varName, diagnostic): '+varName+', '+diagName)
          continue

        # parallelize across binMethods
        for (binVarKey, binMethodName), binMethod in binMethods.items():
          if binMethod.excludeDiag(diagName): continue
          if binMethod.excludeVariable(varName): continue

          if workers is None:
            subStats.append(self._processBinMethod(
              dbValsNN,
              DiagSpaceGrp,
              varName,
              diagName,
              binVarKey, binMethodName, binMethod,
              diagnostic,
            ))
          else:
            subStats.append(workers.apply_async(self._processBinMethod,
              args=(
                dbValsNN,
                DiagSpaceGrp,
                varName,
                diagName,
                binVarKey, binMethodName, binMethod,
                diagnostic,
              )
            ))

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


    for ds in (
      [dataSets['reference'], dataSets['mean']] +
      dataSets['bgEns'] +
      dataSets['anEns'] +
      dataSets['infEns']
    ):
      for typ in ['state', 'diagnostics']:
        if ds[typ] is not None:
          ds[typ].close()

    self.logger.info('Writing statistics file')

    stats = su.BinnedStatisticsFile(statSpace=DiagSpaceName)
    stats.write(statsDict)

    self.logger.info('Finished')

  def _processBinMethod(self,
    dbVals,
    DiagSpaceGrp,
    varName,
    diagName,
    binVarKey, binMethodName, binMethod,
    diagValues,
  ):

    varShort, varUnits = vu.varAttributes(varName)

    #self.logger.info('  starting '+varShort+', '+diagName+', '+binVarKey+', '+binMethodName)

    outerIter = None

    statsDict = {}
    for aa in su.fileStatAttributes:
      statsDict[aa] = []
    for ss in su.allFileStats:
      statsDict[ss] = []

    # initialize binMethod filter function result
    # NOTE: binning can be performed using either mean
    #       or ensemble values, but only limited quantities
    #       are available from ensemble members (e.g., HofX)
    binMethod.evaluate(dbVals, varName, outerIter)

    # initialize binVarShort and binVarUnits for generating the varStatsDict entries
    binVarName, binGrpName = vu.splitObsVarGrp(binVarKey)
    binVarShort, binVarUnits = vu.varAttributes(binVarName)

    binVals = binMethod.getvalues()
    nBins = len(binVals)
    for binVal in binVals:
      # apply binMethod filters for binVal
      binnedDiagnostic = binMethod.apply(diagValues, diagName, binVal)

      # store value and statistics associated with this bin
      statsDict['binVal'].append(binVal)
      statsVal = su.calcStats(binnedDiagnostic)
      for ss in su.allFileStats:
        statsDict[ss].append(statsVal[ss])

    #END binVals LOOP

    # store metadata common to all bins
    statsDict['DiagSpaceGrp'] += [DiagSpaceGrp]*nBins
    statsDict['varName'] += [varShort]*nBins
    statsDict['varUnits'] += [varUnits]*nBins
    statsDict['diagName'] += [diagName]*nBins
    statsDict['binMethod'] += [binMethodName]*nBins
    statsDict['binVar'] += [binVarShort]*nBins
    statsDict['binUnits'] += [binVarUnits]*nBins

    self.logger.info('  completed '+varShort+', '+diagName+', '+binVarKey+', '+binMethodName)

    return statsDict


#=========================================================================
# main program
def main():

  statistics = DiagnoseModelStatistics()
  statistics.diagnose()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
