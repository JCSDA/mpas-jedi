#!/usr/bin/env python3

import argparse
import binning_utils as bu
import binning_configs as bcs
import config as conf
from copy import deepcopy
import diag_utils as du
import fnmatch
import glob
#import math
import logging
import logsetup
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from multiprocessing import Pool
from netCDF4 import Dataset
import numpy as np
import os
import stat_utils as su
import JediDB
import var_utils as vu

_logger = logging.getLogger(__name__)

def writediagstats_obsspace(jdbs, osKey):
    #  jdbs  - list of JediDB objects
    #  osKey - key of jdbs members to reference

    logger = logging.getLogger(__name__+'.'+osKey)

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
    # obs feedback files (JediDB.obsFKey) with the suffix vu.depbgGroup
    obsVars = jdbs[vu.mean].varList(osKey,JediDB.obsFKey,vu.depbgGroup)


    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    binMethods = {}

    for binVarKey, binMethodKeys in binVarConfigs.items():
        binVarConfig = bcs.binVarConfigs.get(binVarKey,bcs.nullBinVarConfig)
        for binMethodKey in binMethodKeys:
            config = binVarConfig.get(binMethodKey,bcs.nullBinMethod).copy()

            if (len(config['values']) < 1 or
                len(config['filters']) < 1): continue

            config['osName'] = ObsSpaceName
            binMethods[(binVarKey,binMethodKey)] = bu.BinMethod(config)


    ######################################
    ## Construct diagnostic configurations
    ######################################

    diagnosticConfigs = du.diagnosticConfigs(
        selectDiagNames, ObsSpaceName, nMembers)


    #####################################################
    ## Generate comprehensive dict of required variables
    #####################################################

    meanDBVars = []
    ensDBVars = []
    dbVars = {vu.mean: [], vu.ensemble: []}
    for varName in obsVars:
        # variables for diagnostics
        for diagName, diagnosticConfig in diagnosticConfigs.items():
            if 'ObsFunction' not in diagnosticConfig: continue
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
                varName, [vu.bgIter, vu.anIter]):
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
    for diagName, diagnosticConfig in diagnosticConfigs.items():
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
#=========================================================================

def main():
    '''
    Wrapper that conducts writediagstats_obsspace across multiple ObsSpaces
    '''
    _logger.info('Starting '+__name__)

    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs", default = 1, type=int,
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-p", "--meanPath", default = JediDB.default_path,
                    help="Path to deterministic or mean state UFO feedback files, default = "
                         +JediDB.default_path)
    ap.add_argument("-o", "--oPrefix", default = JediDB.default_fPrefixes[JediDB.obsFKey],
                    help="prefix for ObsSpace files")
    ap.add_argument("-g", "--gPrefix", default = JediDB.default_fPrefixes[JediDB.geoFKey],
                    help="prefix for GeoVaLs files")
    ap.add_argument("-d", "--dPrefix", default = JediDB.default_fPrefixes[JediDB.diagFKey],
                    help="prefix for ObsDiagnostics files")
    ap.add_argument("-m", "--nMembers", default = 0, type=int,
                    help="number of ensemble members; must be >1 to produce ensemble diagnostic stats")
    ap.add_argument("-e", "--ensemblePath", default = JediDB.default_path+"/mem{:03d}", type=str,
                    help="Path to ensemble member UFO feedback files; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")

    MyArgs = ap.parse_args()

    nprocs = MyArgs.nprocs

    argfPrefixes = {
        JediDB.obsFKey:  MyArgs.oPrefix,
        JediDB.geoFKey:  MyArgs.gPrefix,
        JediDB.diagFKey: MyArgs.dPrefix,
    }

    # construct mean into 0th member slot
    _logger.info('mean database: '+MyArgs.meanPath)
    jdbs = {vu.mean: JediDB.JediDB(MyArgs.meanPath, 'nc4', argfPrefixes)}

    # construct ens members into subsequent slots (when available)
    for member in list(range(1, MyArgs.nMembers+1)):
        ensemblePath = str(MyArgs.ensemblePath).format(member)
        _logger.info('adding member database: '+ensemblePath)
        jdbs[vu.ensSuffix(member)] = JediDB.JediDB(ensemblePath, 'nc4', argfPrefixes)

    # Loop over all experiment+observation combinations (keys) alphabetically
    ospool = Pool(processes=nprocs)

    for osKey in sorted(jdbs[vu.mean].Files):
        res = ospool.apply_async(writediagstats_obsspace, args=(jdbs, osKey))

        # FOR DEBUGGING
#        writediagstats_obsspace(jdbs, osKey)

    ospool.close()
    ospool.join()

    _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
