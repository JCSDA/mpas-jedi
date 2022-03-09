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

    #mean/deterministic states
    #TODO: allow for entire ReferenceStateFile and MPASStateFile as args instead of date
    self.ReferenceStateFile = str(self.args.referenceState)+'.'+fileDate+'.nc'
    self.MPASStateFile = str(self.args.mpasState)+'.'+fileDate+'.nc'

    self.logger.info('ReferenceState: '+self.ReferenceStateFile)
    self.logger.info('MPASState: '+self.MPASStateFile)

    #TODO: enable ensemble statistics similar to DiagnoseObsStatistics
    # for member in list(range(1, self.args.nMembers+1)):
    #   ensemblePath = str(self.args.ensemblePath).format(member)
    #   self.logger.info('adding member database: '+ensemblePath)


  def diagnose(self):

    self.logger.info('Populating binning metadata')

    # initialize netcdf4 Dataset's
    ReferenceState = mu.getNCData(self.ReferenceStateFile)
    MPASState = mu.getNCData(self.MPASStateFile)

    # read "metadata" fields that are used for binning and
    #  create reshaped copies to fit 1D and 2D arrays
    pressure = mu.varRead('pressure', ReferenceState)
    nCells = pressure.shape[0]
    nVertLevelsP1 = pressure.shape[1]+1
    nVertLevels = nVertLevelsP1-1

    dbVals = {}
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)] = {}

    # modVarLev (unique for variables with different numbers of levels)
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)][vu.modVarLev] = np.arange(1, nn+1)

    # modVarLat, modVarLon (same for all level-distributions)
    grid = mu.readGrid(gridFile=self.ReferenceStateFile)
    for nn in [1, nVertLevels, nVertLevelsP1]:
      dbVals[(nn,)][vu.modVarLat] = grid['latitude']
      dbVals[(nn,)][vu.modVarLon] = grid['longitude']

    # add pressure as a binning variable
    # TODO: enable binning by altitude
    # TODO: enable binning by pressure
    # TODO: enable pressure on W levels, if desired
    #dbVals[(nVertLevels,)][vu.modVarPrs] = np.array(pressure)
    #dbVals[(1,)][vu.modVarPrs] = np.array(pressure[:,0]).flatten()


    ###############################################
    ## Extract constructor info about the DiagSpace
    ###############################################

    spaceKey = 'mpas'
    DiagSpaceInfo  = conf.DiagSpaceConfig[spaceKey]
    DiagSpaceGrp   = DiagSpaceInfo['DiagSpaceGrp']
    binVarConfigs = DiagSpaceInfo.get('binVarConfigs',{})

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

            config['osName'] = spaceKey
            config['fileFormat'] = 'model'

            binMethods[(binVarKey, binMethodName)] = bu.BinMethod(config)


    ######################################
    ## Collect statistics for all modelVars
    ######################################

    if self.nprocs > 1:
      workers = mp.Pool(processes = self.nprocs)
    else:
      workers = None

    self.logger.info('Calculating diagnostic statistics')
    self.logger.info("with "+str(self.nprocs)+" out of "+str(mp.cpu_count())+" processors")

    subStats = []
    for varName in modelVars: 

      for diagName in mu.variableSpecificDiagnostics(varName):
        # TODO(JJG): extend to ACC diagnostic, e.g., for 500mb geopotential height

        diagFunction = mu.diagnosticFunctions[diagName]
        modelVarName = mu.aggModelVariable(varName)
        diagnostic = np.asarray(diagFunction(modelVarName, ReferenceState, MPASState))

        dShape = diagnostic.shape
        nDims = len(dShape)
        assert (nDims==1 or nDims==2), 'Number of dimensions must be 1 or 2'

        nn = 1
        if nDims==2:
          nn = dShape[1]
        dbValsNN = dbVals[(nn,)]

        minLevel = mu.aggMinLevel(varName)
        maxLevel = mu.aggMaxLevel(varName, nn)

        # TODO: move this binning to a generic binMethod instead of creating
        #       unique level-binned variables for, e.g., qv

        if len(dShape)==2:
          # mask levels < minLevel
          mask = bu.lessBound(dbValsNN[vu.modVarLev], minLevel)
          diagnostic[:,mask] = np.NaN

          # mask levels > maxLevel
          mask = bu.greatBound(dbValsNN[vu.modVarLev], maxLevel)
          diagnostic[:,mask] = np.NaN

        if diagnostic.size-np.isnan(diagnostic).sum() == 0:
          self.logger.warning('All missing values for diagnostic: '+diagName)
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

    ReferenceState.close()
    MPASState.close()

    self.logger.info('Writing statistics file')

    stats = su.BinnedStatisticsFile(statSpace=spaceKey)
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

    #self.logger.info('  starting '+varShort+', '+diagName+', '+binVarKey+', '+binMethodName)

    varShort, varUnits = vu.varAttributes(varName)
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
