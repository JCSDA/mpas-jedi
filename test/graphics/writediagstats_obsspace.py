import argparse
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
import re
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

#quality control
qcIterBG='0'
NOUTER=os.getenv('NOUTER',1) #set equal to number of outer iterations
qcIterAN=str(NOUTER)

def writediagstats_obsspace(database, osKey):
    #  database - ufofu.FeedbackFiles object
    #  osKey    - key of database members to reference

    LOGPREFIX = this_program+" : "+osKey+" : "

    database.initHandles(osKey)

    ##############################################
    ## Extract constructor info about the ObsSpace
    ##############################################

    ObsSpaceName = database.ObsSpaceName[osKey]
    ObsSpaceInfo = vu.DiagSpaceDict[ObsSpaceName]
    ObsSpaceGrp = ObsSpaceInfo['DiagSpaceGrp']
    binGrps = ObsSpaceInfo.get('binGrps',{})


    ###########################################
    ## Generate comprehensive list of variables
    ##  to read from database
    ###########################################

    # compile list of unique variables needed for binning
    uniqBinVarGrps = []
    # TODO: initialize uniqBinVarGrps with selected binMethods, binArgs,
    #       and binFunc and not 'variables' dictionary member
    for binGrpKey, binMethods in binGrps.items():
        grpDict = vu.binGrpDict.get(binGrpKey,{'variables': []})
        binVarGrps = grpDict['variables']
        for binVarGrp in binVarGrps:
            if binVarGrp not in uniqBinVarGrps:
                uniqBinVarGrps.append(binVarGrp)

    # add iteration numbers to QC varGrp
    for ivar, varGrp in enumerate(uniqBinVarGrps):
        if varGrp==vu.selfQCValue:
            uniqBinVarGrps[ivar] = varGrp+qcIterBG
            uniqBinVarGrps.insert(ivar+1,varGrp+qcIterAN)

    # create observed variable list by selecting those variables in the
    # obs feedback files (ufofu.obsFKey) with the suffix vu.depbgGroup
    varlist, chanlist = database.varList(osKey,ufofu.obsFKey,vu.depbgGroup)

    # replace binning variables that
    # (1) contain vu.vNameStr
    # (2) contain vu.vChanStr
    for ivar, varGrp in enumerate(uniqBinVarGrps):
        binVarName, binGrpName = vu.splitObsVarGrp(varGrp)
        if (binVarName == vu.vNameStr
            or vu.vChanStr in binVarName.split("_")):
            for jvar, varName in enumerate(varlist):
                varGrpSub = re.sub(vu.vNameStr,varName,varGrp)
                varGrpSub = re.sub(vu.vChanStr,chanlist[jvar],varGrpSub)
                uniqBinVarGrps.insert(ivar+1,varGrpSub)
            uniqBinVarGrps.remove(varGrp)

    # add unique obs and departure variables to combined list
    uniqFileVarGrps = deepcopy(uniqBinVarGrps)
    for varName in varlist:
        for grpName in [vu.obsGroup, vu.depbgGroup, vu.depanGroup]:
            varGrp   = varName+'@'+grpName
            if varGrp not in uniqFileVarGrps:
                uniqFileVarGrps.append(varGrp)
    #print(uniqFileVarGrps)

    # read the variables into memory
    fileVarsVals = database.readVars(osKey,uniqFileVarGrps)


    #################################################
    ## Collect statistics for each varName in varlist
    #################################################

    # Initialize a dictionary to contain all statistical info for this osKey
    statsDict = {}
    for attribName in su.fileStatAttributes:
        statsDict[attribName] = []
    for statName in su.allFileStats:
        statsDict[statName] = []

    print(LOGPREFIX+"Calculating/writing diagnostic stats for:")
    for ivar, varName in enumerate(varlist):
        print(LOGPREFIX+"VARIABLE = "+varName)
        obsName   = varName+'@'+vu.obsGroup
        depbgName = varName+'@'+vu.depbgGroup
        depanName = varName+'@'+vu.depanGroup

        chanStr = chanlist[ivar]
        if ObsSpaceGrp == vu.radiance_s:
            dictName = '_'.join(varName.split("_")[:-1])
            varVal = vu.varDictObs.get(dictName,['',dictName])
            varShort = varVal[1]+'_'+chanStr
        else:
            varVal = vu.varDictObs.get(varName,['',varName])
            varShort = varVal[1]
        varUnits = varVal[0]


        # add binning arguments that are functions of
        # variables already contained within fileVarsVals:
        # TODO: move this to bin variables constructor stage
        #       using a new binFunction class and binFunction.argVars member
        # (1) varName@vu.bakGroup
        varGrp = varName+'@'+vu.bakGroup
        if varGrp in uniqBinVarGrps:
            fileVarsVals[varGrp] = np.add(
                fileVarsVals[depbgName],
                fileVarsVals[obsName])


        # collect stats for all selectDiagNames
        for diagName in selectDiagNames:
            print(LOGPREFIX+"DIAG = "+diagName)
            if diagName == 'omb':
                Diagnostic = np.negative(fileVarsVals[depbgName])
            elif diagName == 'oma':
                Diagnostic = np.negative(fileVarsVals[depanName])
            elif diagName == 'obs':
                Diagnostic = deepcopy(fileVarsVals[obsName])
            elif diagName == 'bak':
                Diagnostic = np.add(fileVarsVals[depbgName], fileVarsVals[obsName])
            elif diagName == 'ana':
                Diagnostic = np.add(fileVarsVals[depanName], fileVarsVals[obsName])
            else:
                print('\n\nERROR: diagName is undefined: '+diagName)
                os._exit(1)

            if ((diagName == 'omb' or diagName == 'oma')
               and 'gnssro' in ObsSpaceName):
                Diagnostic = np.multiply(np.divide(Diagnostic,fileVarsVals[obsName]),100.)


            # filter and store Diagnostic statistics for each binGrp
            for binGrpKey, binMethods in binGrps.items():
                grpDict = vu.binGrpDict.get(binGrpKey,{'variables': [vu.miss_s]})
                if grpDict['variables'][0] == vu.miss_s: continue

                binVarName = grpDict.get('binVarName',vu.miss_s)
                binVar, binGrp = vu.splitObsVarGrp(binVarName)
                binVarDict  = vu.varDictObs.get(binVar,[vu.miss_s,binVar])
                binVarUnits = binVarDict[0]
                binVarShort = binVarDict[1]

                for ikey, binMethod in enumerate(binMethods):
                    keyDesc = grpDict.get(binMethod,vu.nullBinMethod)
                    binVals = keyDesc['labels']
                    if binVals[0] == vu.miss_s: continue

                    binFilters = keyDesc['filters']
                    for ibin, binVal in enumerate(binVals):
                        # filter data that meets bin criteria
                        binnedDiagnostic = deepcopy(Diagnostic)
                        for binFilter in binFilters:
                            applyToDiags = binFilter.get('apply_to',vu.allDiags)
                            if not (diagName in applyToDiags): continue

                            binWhere   = binFilter['where']
                            binFunc    = binFilter.get('binFunc',vu.firstKey)
                            binArgs    = binFilter['binArgs']
                            binArgsMod = deepcopy(binArgs)
                            # substitute relevant varName and channel as needed
                            for jvar, varGrp in enumerate(binArgsMod):
                                varGrpSub = re.sub(vu.vNameStr,varName,varGrp)
                                varGrpSub = re.sub(vu.vChanStr,chanStr,varGrpSub)
                                if (diagName == 'oma' or diagName == 'ana'):
                                    varGrpSub = re.sub(vu.qcGroup,vu.qcGroup+qcIterAN,varGrpSub)
                                else:
                                    varGrpSub = re.sub(vu.qcGroup,vu.qcGroup+qcIterBG,varGrpSub)
                                binArgsMod[jvar] = varGrpSub

                            binArgsVals = {}
                            for iarg, arg in enumerate(binArgs):
                                binArgsVals[arg] = fileVarsVals[binArgsMod[iarg]]

                            binBound = (binFilter['bounds'])[ibin]

                            # blacklist locations where the mask is True
                            mask = binWhere(binFunc(binArgsVals),binBound)
                            maskValue = binFilter.get('mask_value',np.NaN)
                            if len(mask) == len(binnedDiagnostic):
                                binnedDiagnostic[mask] = maskValue
                            else:
                                print('\n\nERROR: mask is incorrectly defined for '
                                              +diagName+" "+binGrpKey+" "+binMethod+" "+binVal)
                                os._exit(1)
                        #END binFilters LOOP

                        statsDict['DiagSpaceGrp'].append(ObsSpaceGrp)
                        statsDict['varName'].append(varShort)
                        statsDict['varUnits'].append(varUnits)
                        statsDict['diagName'].append(diagName)
                        statsDict['binMethod'].append(binMethod)
                        statsDict['binVar'].append(binVarShort)
                        statsDict['binUnits'].append(binVarUnits)
                        statsDict['binVal'].append(binVal)

                        statsVal = su.calcStats(binnedDiagnostic)
                        for statName in su.allFileStats:
                            statsDict[statName].append(statsVal[statName])

                    #END binVals LOOP
                #END binMethods LOOP
            #END binGrps LOOP
        #END selectDiagNames LOOP
    #END varlist LOOP


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

    ospool = Pool(processes=nprocs)

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
    for osKey in sorted(database.Files):
        res = ospool.apply_async(writediagstats_obsspace, args=(database,osKey))

#        # FOR DEBUGGING
#        writediagstats_obsspace(database,osKey)

    ospool.close()
    ospool.join()

if __name__ == '__main__': main()
