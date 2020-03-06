import argparse
import binning_utils as bu
import config as conf
from copy import deepcopy
import fnmatch
import glob
#import math
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from multiprocessing import Pool
from netCDF4 import Dataset
import numpy as np
import os
import par_utils as paru
import stat_utils as su
import sys
import ufo_file_utils as ufofu
import var_utils as vu

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See par_utils.par_args for more information.

#Select the diagnostics for which statistics are calculated
# options: 'omb','oma','obs','bak','ana'
selectDiagNames = ['omb','oma']

this_program = os.path.basename(sys.argv[0])

#analysis outer iteration
NOUTER=os.getenv('NOUTER',1) #set equal to number of outer iterations
anIter=str(NOUTER)

def writediagstats_obsspace(database, osKey):
    #  database - ufofu.FeedbackFiles object
    #  osKey    - key of database members to reference

    LOGPREFIX = this_program+" : "+osKey+" : "

    database.initHandles(osKey)


    ##############################################
    ## Extract constructor info about the ObsSpace
    ##############################################

    ObsSpaceName  = database.ObsSpaceName[osKey]
    ObsSpaceInfo  = conf.DiagSpaceConfig[ObsSpaceName]
    ObsSpaceGrp   = ObsSpaceInfo['DiagSpaceGrp']
    binVarConfigs = ObsSpaceInfo.get('binVarConfigs',{})

    # create observed variable list by selecting those variables in the
    # obs feedback files (ufofu.obsFKey) with the suffix vu.depbgGroup
    obsVars = database.varList(osKey,ufofu.obsFKey,vu.depbgGroup)


    ########################################################
    ## Construct dictionary of binMethods for this ObsSpace
    ########################################################

    binMethods = {}

    for binVarKey, binMethodKeys in binVarConfigs.items():
        binVarConfig = bu.binVarConfigs.get(binVarKey,bu.nullBinVarConfig)
        for binMethodKey in binMethodKeys:
            binMethodConfig = binVarConfig.get(binMethodKey,bu.nullBinMethod)

            if (len(binMethodConfig['values']) < 1 or
                len(binMethodConfig['filters']) < 1): continue

            binMethods[(binVarKey,binMethodKey)] = bu.BinMethod(binMethodConfig,ObsSpaceName)


    ############################################
    ## Generate comprehensive list of variables
    ##  and then read from database
    ############################################

    dbVars = []

    for varName in obsVars:
        # add obs and departure variables
        for grpName in [vu.obsGroup, vu.depbgGroup, vu.depanGroup]:
            varGrp   = varName+'@'+grpName
            dbVars.append(varGrp)

        # add all variables needed for binning
        for (binVarKey,binMethodKey), binMethod in binMethods.items():
            for varGrp in binMethod.dbVars(
                varName,[vu.bgIter,anIter]):
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

    print(LOGPREFIX+"Calculating/writing diagnostic stats for:")
    for varName in obsVars:
        print(LOGPREFIX+"VARIABLE = "+varName)
        obsName   = varName+'@'+vu.obsGroup
        depbgName = varName+'@'+vu.depbgGroup
        depanName = varName+'@'+vu.depanGroup

        varShort, varUnits = vu.varAttributes(varName)

        # collect stats for all selectDiagNames
        for diagName in selectDiagNames:
            print(LOGPREFIX+"DIAG = "+diagName)
            outerIter = vu.bgIter
            # TODO: calculate diagnostics using *FilterFunc classes
            #  --> if more complex diagnostics are needed
            if diagName == 'omb':
                Diagnostic = np.negative(dbVals[depbgName])
            elif diagName == 'oma':
                Diagnostic = np.negative(dbVals[depanName])
                outerIter = anIter
            elif diagName == 'obs':
                Diagnostic = deepcopy(dbVals[obsName])
            elif diagName == 'bak':
                Diagnostic = np.add(dbVals[depbgName], dbVals[obsName])
            elif diagName == 'ana':
                Diagnostic = np.add(dbVals[depanName], dbVals[obsName])
                outerIter = anIter
            else:
                print('\n\nERROR: diagName is undefined: '+diagName)
                os._exit(1)

            if ((diagName == 'omb' or diagName == 'oma')
               and 'gnssro' in ObsSpaceName):
                Diagnostic = np.multiply(np.divide(Diagnostic,dbVals[obsName]),100.)

            for (binVarKey,binMethodKey), binMethod in binMethods.items():
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

#        # FOR DEBUGGING
#        writediagstats_obsspace(database,osKey)

    ospool.close()
    ospool.join()

if __name__ == '__main__': main()
