#!/usr/bin/env python3

import DiagnoseObsStatisticsArgs

import binning_utils as bu
import predefined_configs as pconf
import config as conf
from copy import deepcopy
import diag_utils as du
from JediDB import JediDB
from JediDBArgs import obsFKey
import logging
import logsetup
import multiprocessing as mp
import numpy as np
import os
import stat_utils as su
import var_utils as vu

_logger = logging.getLogger(__name__)

class DiagnoseObsStatistics:
  '''
  Diagnose observation-space statistics
  Driven by
    - static selections in conf
    - command-line arguments in DiagnoseObsStatisticsArgs
  '''
  def __init__(self):
    self.name = 'DiagnoseObsStatistics'
    self.args = DiagnoseObsStatisticsArgs.args
    self.logger = logging.getLogger(self.name)
    self.nprocs = min(mp.cpu_count(), self.args.nprocs)

    # construct mean DB into 0th member slot
    self.logger.info('mean database: '+self.args.meanPath)
    self.jdbs = {vu.mean: JediDB(self.args.meanPath)}
    self.osKeys = sorted(self.jdbs[vu.mean].Files.keys())

    # construct ens member DBs into subsequent slots (when available)
    for member in list(range(1, self.args.nMembers+1)):
        ensemblePath = str(self.args.ensemblePath).format(member)
        self.logger.info('adding member database: '+ensemblePath)
        self.jdbs[vu.ensSuffix(member)] = JediDB(ensemblePath)

  def diagnose(self, workers = None):
    '''
    conducts diagnoseObsSpace across multiple ObsSpaces in parallel
    '''
    # Loop over all experiment+observation combinations (keys) alphabetically
    for osKey in self.osKeys:
      self.logger.info(osKey)
      if workers is None:
        self.__diagnoseObsSpace(self.jdbs, osKey)
      else:
        self.nprocs = 1
        res = workers.apply_async(self.__diagnoseObsSpace, args=(self.jdbs, osKey))

  def __diagnoseObsSpace(self, jdbs, osKey):
    #  osKey - key of jdbs members to reference
    logger = logging.getLogger(self.name+'.diagnoseObsSpace('+osKey+')')
    nMembers = len(jdbs)-1

    # initialize mean db file handles
    jdbs[vu.mean].initHandles(osKey)

    ###############################################
    ## Extract constructor info about the ObsSpace
    ###############################################

    ObsSpaceName  = jdbs[vu.mean].ObsSpaceName[osKey]
    ObsSpaceInfo  = conf.DiagSpaceConfig[ObsSpaceName]
    ObsSpaceGrp   = ObsSpaceInfo['DiagSpaceGrp']
    binVarConfigs = ObsSpaceInfo.get('binVarConfigs',{})
    selectDiagNames = ObsSpaceInfo.get('diagNames',{})

    # create observed variable list by selecting those variables in the
    # obs feedback files (obsFKey) with the proper suffix
    if self.args.jediAppName == 'variational':
        markerGroup = vu.ombgGroup
    elif self.args.jediAppName == 'hofx':
        markerGroup = vu.hofxGroup
    else:
        logger.error('JEDI Application is not supported:: '+self.args.jediAppName)
    obsVars = jdbs[vu.mean].varList(osKey, obsFKey, markerGroup)

    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    logger.info('Initializing binMethods')

    binMethods = {}

    for binVarKey, binMethodNames in binVarConfigs.items():
        binVarConfig = pconf.binVarConfigs.get(binVarKey,pconf.nullBinVarConfig)
        for binMethodName in binMethodNames:
            config = binVarConfig.get(binMethodName,pconf.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['dsName'] = ObsSpaceName
            config['fileFormat'] = jdbs[vu.mean].fileFormat(osKey, obsFKey)

            binMethods[(binVarKey, binMethodName)] = bu.BinMethod(config)


    ######################################
    ## Construct diagnostic configurations
    ######################################

    logger.info('Initializing diagnosticConfigs')

    diagnosticConfigs = du.diagnosticConfigs(
        selectDiagNames, ObsSpaceName,
        includeEnsembleDiagnostics = (nMembers > 1),
        fileFormat = jdbs[vu.mean].fileFormat(osKey, obsFKey))


    #####################################################
    ## Generate comprehensive dict of required variables
    #####################################################

    logger.info('Initializing dbVars')

    meanDBVars = []
    ensDBVars = []
    dbVars = {vu.mean: [], vu.ensemble: []}
    for varName in obsVars:
        for diagName, diagnosticConfig in diagnosticConfigs.items():
            if 'BinFunction' not in diagnosticConfig: continue

            # variables for diagnostics
            for grpVar in diagnosticConfig['BinFunction'].dbVars(
                varName, diagnosticConfig['outerIter']):
                for memberType in dbVars.keys():
                    if diagnosticConfig[memberType]:
                        dbVars[memberType].append(grpVar)

            # variables for binning
            # TODO: anIter grpVar's are not needed for all applications
            #       can save some reading time+memory by checking all diagnosticConfigs
            #       for required iterations before appending to dbVars[vu.mean] below
            for key, binMethod in binMethods.items():
                for grpVar in binMethod.dbVars(
                    varName, diagnosticConfig['outerIter']):
                    dbVars[vu.mean].append(grpVar)


    #####################################
    ## Read required variables from jdbs
    #####################################

    logger.info('Reading dbVals')

    # read mean database variable values into memory
    dbVals = jdbs[vu.mean].readVars(osKey, dbVars[vu.mean], self.nprocs)

    # destroy mean file handles
    jdbs[vu.mean].destroyHandles(osKey)

    # now for ensemble members
    # TODO: read members in parallel
    for memStr, jdb in jdbs.items():
        if memStr == vu.mean: continue

        # initialize member db file handles
        jdb.initHandles(osKey)

        # read database variable values into memory
        memberDBVals = jdb.readVars(osKey, dbVars[vu.ensemble], self.nprocs)
        for dbVar, vals in memberDBVals.items():
            dbVals[dbVar+memStr] = vals.copy()

        # destroy file handles
        jdb.destroyHandles(osKey)


    ######################################
    ## Collect statistics for all obsVars
    ######################################

    # Initialize a dictionary to contain all statistical info for this osKey

    if self.nprocs > 1:
        workers = mp.Pool(processes = self.nprocs)
    else:
        workers = None

    logger.info('Calculating diagnostic statistics')
    logger.info("with "+str(self.nprocs)+" out of "+str(mp.cpu_count())+" processors")

    subStats = []
    for varName in obsVars:
      # collect stats for all diagnosticConfigs
      for diagName, diagnosticConfig in sorted(diagnosticConfigs.items()):
        if 'BinFunction' not in diagnosticConfig: continue

        if workers is None:
          subStats.append(self._processBinMethods(
            dbVals,
            ObsSpaceGrp,
            varName,
            diagName, diagnosticConfig,
            binMethods,
            logger,
          ))
        else:
          subStats.append(workers.apply_async(self._processBinMethods,
            args=(
              dbVals,
              ObsSpaceGrp,
              varName,
              diagName, diagnosticConfig,
              binMethods,
              logger,
            )
          ))
      #END diagnosticConfigs LOOP
    #END obsVars LOOP

    statsDict = {}
    for attribName in su.fileStatAttributes:
        statsDict[attribName] = []
    for statName in su.allFileStats:
        statsDict[statName] = []

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

    del dbVals, binMethods, diagnosticConfigs

    ## Create a new stats file for osKey
    logger.info('Writing statistics file')
    stats = su.BinnedStatisticsFile(statSpace=osKey)
    stats.write(statsDict)

    logger.info('Finished')


  #@staticmethod
  def _processBinMethods(self,
    dbVals,
    ObsSpaceGrp,
    varName,
    diagName, diagnosticConfig,
    binMethods,
    logger,
  ):

    varShort, varUnits = vu.varAttributes(varName)

    diagFunction = diagnosticConfig['BinFunction']
    outerIter = diagnosticConfig['outerIter']

    diagFunction.evaluate(dbVals, varName, outerIter)
    diagValues = diagFunction.result

    if len(diagValues)-np.isnan(diagValues).sum() == 0:
      logger.warning('All missing values for diagnostic: '+diagName)

    statsDict = {}
    for aa in su.fileStatAttributes:
      statsDict[aa] = []
    for ss in su.allFileStats:
      statsDict[ss] = []

    for (binVarKey, binMethodName), binMethod in binMethods.items():
      if binMethod.excludeDiag(diagName): continue
      if binMethod.excludeVariable(varName): continue

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
        for statName in su.allFileStats:
          statsDict[statName].append(statsVal[statName])

      #END binMethod.values LOOP

      # store metadata common to all bins
      statsDict['DiagSpaceGrp'] += [ObsSpaceGrp]*nBins
      statsDict['varName'] += [varShort]*nBins
      statsDict['varUnits'] += [varUnits]*nBins
      statsDict['diagName'] += [diagName]*nBins
      statsDict['binMethod'] += [binMethodName]*nBins
      statsDict['binVar'] += [binVarShort]*nBins
      statsDict['binUnits'] += [binVarUnits]*nBins

      logger.info('  completed '+varShort+', '+diagName+', '+binVarKey+', '+binMethodName)

    #END binMethods tuple LOOP


    return statsDict


#=========================================================================
# main program
def main():
  _logger.info('Starting '+__name__)

  statistics = DiagnoseObsStatistics()
  statistics.diagnose()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
