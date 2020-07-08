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
import ufo_file_utils as ufofu
import var_utils as vu

#analysis outer iteration
anIter=str(du.NOUTER)

_logger = logging.getLogger(__name__)

def writediagstats_obsspace(database, osKey):
    #  database - ufofu.FeedbackFiles object
    #  osKey    - key of database members to reference

    logger = logging.getLogger(__name__+'.'+osKey)

    database.initHandles(osKey)

    ###############################################
    ## Extract constructor info about the ObsSpace
    ###############################################

    ObsSpaceName  = database.ObsSpaceName[osKey]
    ObsSpaceInfo  = conf.DiagSpaceConfig[ObsSpaceName]
    ObsSpaceGrp   = ObsSpaceInfo['DiagSpaceGrp']
    binVarConfigs = ObsSpaceInfo.get('binVarConfigs',{})
    selectDiagNames = ObsSpaceInfo.get('diagNames',{})

    # create observed variable list by selecting those variables in the
    # obs feedback files (ufofu.obsFKey) with the suffix vu.depbgGroup
    obsVars = database.varList(osKey,ufofu.obsFKey,vu.depbgGroup)


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


    ##################################
    ## Construct diagnostic functions
    ##################################

    diagFunctions = {}
    for diagName in selectDiagNames:
        config = deepcopy(du.availableDiagnostics.get(diagName,None))
        if config is None:
            logger.error('diagName is undefined: '+diagName)

        osNames = config.get('osNames',[])
        if (len(osNames) > 0 and
            ObsSpaceName not in osNames): continue

        config['osName'] = ObsSpaceName

        outerIter = '0'
        outerIterStr = config.get('iter',None)
        if outerIterStr is not None:
            if outerIterStr == 'bg':
                outerIter = vu.bgIter
            elif outerIterStr == 'an':
                outerIter = anIter
            elif pu.isint(outerIterStr):
                outerIter = outerIterStr
            else:
                logger.error('outerIter is undefined: '+outerIterStr)

        diagFunctions[diagName] = {
            'ObsFunction': bu.ObsFunctionWrapper(config),
            'outerIter': outerIter,
        }

    ############################################
    ## Generate comprehensive list of required
    ## variables, then read from database
    ############################################

    dbVars = []

    for varName in obsVars:
        # variables for diagnostics
        for diagName, diagFunction in diagFunctions.items():
            for varGrp in diagFunction['ObsFunction'].dbVars(
                varName, diagFunction['outerIter']):
                dbVars.append(varGrp)

        # variables for binning
        for (binVarKey,binMethodKey), binMethod in binMethods.items():
            for varGrp in binMethod.dbVars(
                varName, [vu.bgIter,anIter]):
                dbVars.append(varGrp)

    # read the database values into memory
    dbVals = database.readVars(osKey,dbVars)


    ######################################
    ## Collect statistics for all obsVars
    ######################################

    # Initialize a dictionary to contain all statistical info for this osKey
    statsDict = {}
    for attribName in su.fileStatAttributes:
        statsDict[attribName] = []
    for statName in su.allFileStats:
        statsDict[statName] = []

    logger.info('Calculating/writing diagnostic stats for:')
    for varName in obsVars:
        logger.info('VARIABLE = '+varName)

        varShort, varUnits = vu.varAttributes(varName)

        # collect stats for all diagFunctions
        for diagName, diagFunction in diagFunctions.items():
            logger.info('DIAG = '+diagName)

            func = diagFunction['ObsFunction']
            outerIter = diagFunction['outerIter']
            func.evaluate(dbVals,varName,outerIter)
            Diagnostic = func.result

            for (binVarKey,binMethodKey), binMethod in binMethods.items():
                if diagName in binMethod.excludeDiags: continue

                binVarName, binGrpName = vu.splitObsVarGrp(binVarKey)
                binVarShort, binVarUnits = vu.varAttributes(binVarName)

                # initialize binMethod filter function result
                binMethod.evaluate(dbVals,varName,outerIter)

                for binVal in binMethod.values:
                    # store metadata common to all bins
                    statsDict['DiagSpaceGrp'].append(ObsSpaceGrp)
                    statsDict['varName'].append(varShort)
                    statsDict['varUnits'].append(varUnits)
                    statsDict['diagName'].append(diagName)
                    statsDict['binMethod'].append(binMethodKey)
                    statsDict['binVar'].append(binVarShort)
                    statsDict['binUnits'].append(binVarUnits)

                    # apply binMethod filters for binVal
                    binnedDiagnostic = binMethod.apply(Diagnostic,diagName,binVal)

                    # store value and statistics associated with this bin
                    statsDict['binVal'].append(binVal)
                    statsVal = su.calcStats(binnedDiagnostic)
                    for statName in su.allFileStats:
                        statsDict[statName].append(statsVal[statName])

                #END binMethod.values LOOP
            #END binMethods tuple LOOP
        #END selectDiagNames LOOP
    #END obsVars LOOP

    ## Create a new stats file for osKey
    su.write_stats_nc(osKey,statsDict)

    ## Destroy UFO file handles
    database.destroyHandles(osKey)


#=========================================================================
#=========================================================================

def main():
    '''
    Wrapper that conducts writediagstats_obsspace across multiple ObsSpaces
    '''
    _logger.info('Starting '+__name__)

    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-p", "--data_path",
                    help="Path to UFO feedback files, default = "
                         +ufofu.default_path)
    ap.add_argument("-o", "--oPrefix",
                    help="prefix for ObsSpace files")
    ap.add_argument("-g", "--gPrefix",
                    help="prefix for GeoVaLs files")
    ap.add_argument("-d", "--dPrefix",
                    help="prefix for ObsDiagnostics files")

    MyArgs = ap.parse_args()

    if MyArgs.nprocs:
        nprocs = int(MyArgs.nprocs)
    else:
        nprocs = 1

    data_path = ufofu.default_path
    if MyArgs.data_path: data_path = MyArgs.data_path

    argfPrefixes = {
        ufofu.obsFKey:  MyArgs.oPrefix,
        ufofu.geoFKey:  MyArgs.gPrefix,
        ufofu.diagFKey: MyArgs.dPrefix,
    }
    for key, prefix in ufofu.default_fPrefixes.items():
        if not argfPrefixes.get(key,False):
            argfPrefixes[key] = ufofu.default_fPrefixes[key]

    database = ufofu.FeedbackFiles(data_path,'nc4',argfPrefixes)

    # Loop over all experiment+observation combinations (keys) alphabetically
    ospool = Pool(processes=nprocs)

    for osKey in sorted(database.Files):
        res = ospool.apply_async(writediagstats_obsspace, args=(database,osKey))

        # FOR DEBUGGING
#        writediagstats_obsspace(database,osKey)

    ospool.close()
    ospool.join()

    _logger.info('Finished '+__name__+' successfully')

if __name__ == '__main__': main()
