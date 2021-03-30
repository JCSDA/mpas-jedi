#!/usr/bin/env python3

import DiagnoseObsStatisticsArgs

import binning_utils as bu
import predefined_configs as pconf
import config as conf
from copy import deepcopy
import diag_utils as du
import fnmatch
import glob
from JediDB import JediDB
from JediDBArgs import obsFKey
import logging
import logsetup
import multiprocessing as mp
from netCDF4 import Dataset
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

    # construct mean DB into 0th member slot
    self.logger.info('mean database: '+self.args.meanPath)
    self.jdbs = {vu.mean: JediDB(self.args.meanPath, 'nc4')}
    self.osKeys = sorted(self.jdbs[vu.mean].Files.keys())

    # construct ens member DBs into subsequent slots (when available)
    for member in list(range(1, self.args.nMembers+1)):
        ensemblePath = str(self.args.ensemblePath).format(member)
        self.logger.info('adding member database: '+ensemblePath)
        self.jdbs[vu.ensSuffix(member)] = JediDB(ensemblePath, 'nc4')

  def diagnose(self, workers = None):
    '''
    conducts diagnoseObsSpace across multiple ObsSpaces in parallel
    '''
    # Loop over all experiment+observation combinations (keys) alphabetically
    for osKey in self.osKeys:
      self.logger.info(osKey)
      if workers is None:
        self.diagnoseObsSpace(self.jdbs, osKey)
      else:
        res = workers.apply_async(self.diagnoseObsSpace, args=(self.jdbs, osKey))

  def diagnoseObsSpace(self, jdbs, osKey):
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
        markerSuffix = vu.depbgGroup
    elif self.args.jediAppName == 'hofx':
        markerSuffix = vu.hofxGroup
    else:
        logger.error('JEDI Application is not supported:: '+self.args.jediAppName)
    obsVars = jdbs[vu.mean].varList(osKey, obsFKey, markerSuffix)

    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    binMethods = {}

    for binVarKey, binMethodKeys in binVarConfigs.items():
        binVarConfig = pconf.binVarConfigs.get(binVarKey,pconf.nullBinVarConfig)
        for binMethodKey in binMethodKeys:
            config = binVarConfig.get(binMethodKey,pconf.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['osName'] = ObsSpaceName
            binMethods[(binVarKey,binMethodKey)] = bu.BinMethod(config)


    ######################################
    ## Construct diagnostic configurations
    ######################################

    diagnosticConfigs = du.diagnosticConfigs(
        selectDiagNames, ObsSpaceName,
        includeEnsembleDiagnostics = (nMembers > 1))


    #####################################################
    ## Generate comprehensive dict of required variables
    #####################################################

    meanDBVars = []
    ensDBVars = []
    dbVars = {vu.mean: [], vu.ensemble: []}
    for varName in obsVars:
        for diagName, diagnosticConfig in diagnosticConfigs.items():
            if 'ObsFunction' not in diagnosticConfig: continue

            # variables for diagnostics
            for varGrp in diagnosticConfig['ObsFunction'].dbVars(
                varName, diagnosticConfig['outerIter']):
                for memberType in dbVars.keys():
                    if diagnosticConfig[memberType]:
                        dbVars[memberType].append(varGrp)

            # variables for binning
            # TODO: anIter varGrp's are not needed for all applications
            #       can save some reading time+memory by checking all diagnosticConfigs
            #       for required iterations before appending to dbVars[vu.mean] below
            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                for varGrp in binMethod.dbVars(
                    varName, diagnosticConfig['outerIter']):
                    dbVars[vu.mean].append(varGrp)


    #####################################
    ## Read required variables from jdbs
    #####################################

    # read mean database variable values into memory
    dbVals = jdbs[vu.mean].readVars(osKey, dbVars[vu.mean])

    # destroy mean file handles
    jdbs[vu.mean].destroyHandles(osKey)

    # now for ensemble members
    for memStr, jdb in jdbs.items():
        if memStr == vu.mean: continue

        # initialize member db file handles
        jdb.initHandles(osKey)

        # read database variable values into memory
        memberDBVals = jdb.readVars(osKey, dbVars[vu.ensemble])
        for dbVar, vals in memberDBVals.items():
            dbVals[dbVar+memStr] = vals.copy()

        # destroy file handles
        jdb.destroyHandles(osKey)


    ######################################
    ## Collect statistics for all obsVars
    ######################################

    # Initialize a dictionary to contain all statistical info for this osKey
    statsDict = {}
    for attribName in su.fileStatAttributes:
        statsDict[attribName] = []
    for statName in su.allFileStats:
        statsDict[statName] = []

    # collect stats for all diagnosticConfigs
    for diagName, diagnosticConfig in sorted(diagnosticConfigs.items()):
        if 'ObsFunction' not in diagnosticConfig: continue

        logger.info('Calculating/writing diagnostic stats for:')
        logger.info('DIAG = '+diagName)
        Diagnostic = diagnosticConfig['ObsFunction']
        outerIter = diagnosticConfig['outerIter']

        for varName in obsVars:
            logger.info('VARIABLE = '+varName)

            varShort, varUnits = vu.varAttributes(varName)

            Diagnostic.evaluate(dbVals, varName, outerIter)
            diagValues = Diagnostic.result

            if len(diagValues)-np.isnan(diagValues).sum() == 0:
                logger.warning('All missing values for diagnostic: '+diagName)

            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                if diagName in binMethod.excludeDiags: continue

                binVarName, binGrpName = vu.splitObsVarGrp(binVarKey)
                binVarShort, binVarUnits = vu.varAttributes(binVarName)

                # initialize binMethod filter function result
                # NOTE: binning is performed using mean values
                #       and not ensemble member values
                binMethod.evaluate(dbVals, varName, outerIter)

                for binVal in binMethod.values:
                    # apply binMethod filters for binVal
                    binnedDiagnostic = binMethod.apply(diagValues,diagName,binVal)

                    # store value and statistics associated with this bin
                    statsDict['binVal'].append(binVal)
                    statsVal = su.calcStats(binnedDiagnostic)
                    for statName in su.allFileStats:
                        statsDict[statName].append(statsVal[statName])

                    # store metadata common to all bins
                    statsDict['DiagSpaceGrp'].append(ObsSpaceGrp)
                    statsDict['varName'].append(varShort)
                    statsDict['varUnits'].append(varUnits)
                    statsDict['diagName'].append(diagName)
                    statsDict['binMethod'].append(binMethodKey)
                    statsDict['binVar'].append(binVarShort)
                    statsDict['binUnits'].append(binVarUnits)

                #END binMethod.values LOOP
            #END binMethods tuple LOOP
        #END obsVars LOOP
    #END diagnosticConfigs LOOP

    ## Create a new stats file for osKey
    su.write_stats_nc(osKey,statsDict)

#=========================================================================
# main program
def main():
  _logger.info('Starting '+__name__)

  statistics = DiagnoseObsStatistics()

  # create pool of workers
  workers = mp.Pool(processes = statistics.args.nprocs)

  # diagnose statistics
  statistics.diagnose(workers)

  # wait for workers to finish
  workers.close()
  workers.join()

  _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
