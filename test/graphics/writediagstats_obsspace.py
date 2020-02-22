from copy import deepcopy
import fnmatch
import glob
#import math
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from netCDF4 import Dataset
import numpy as np
import os
import par_utils as paru
import re
import stat_utils as su
import sys
import var_utils as vu

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See par_utils.par_args for more information.

#Select the diagnostics for which statistics are calculated
# options: 'omb','oma','obs','bak','ana'
selectDiagNames = ['omb','oma']

# Diagnostic oma/omb files located in diagDir
diagDir           = '../Data/'
obsOutFilePrefix  = 'obsout_'
obsDiagFilePrefix = 'ydiags_'
fileExt           = 'nc4'
this_program      = 'writediagstats_obsspace'

#quality control
qcIterBG='0'
NOUTER=os.getenv('NOUTER',1) #set equal to number of outer iterations
qcIterAN=str(NOUTER)

ObsSpaceFileKey = 'ObsSpaceFiles'
ObsDiagFileKey  = 'ObsDiagFiles'
geovalFileKey   = 'geovalFiles'

filePrefixes = {}
filePrefixes[ObsSpaceFileKey] = 'obsout_'
filePrefixes[ObsDiagFileKey]  = 'ydiags_'
filePrefixes[geovalFileKey]   = 'geoval_'

def writediagstats_obsspace(nproc=1, myproc=0):
#  nproc - total number of processors
#  myproc - processor rank, starting at 0
#  Note: these arguments are used when the calling program has multiple processors available.

    # catalog all relevant files
    allFiles = {}
    for key, prefix in filePrefixes.items():
        allFiles[key] = []
        for File in glob.glob(diagDir+prefix+'*.'+fileExt):
            allFiles[key].append(File)

    # Group files by experiment-ObsSpace combination
    #  Must fit the format filePrefixes[ObsSpaceFileKey]+"_"+expObsSpaceCombo+"_"+PE+"."+fileExt,
    #  where each group contains files from all PE's,
    #  which is the 4-digit processor rank [0000 to XXXX]
    #  examples:
    #   +      obsout_3dvar_sondes_0000.nc4, where
    #   expObsSpaceCombo = "3dvar_sondes"
    #
    #   +      obsout_3dvar_aircraft_0000.nc4, where
    #   expObsSpaceCombo = "3dvar_aircraft"
    #
    #   +      obsout_3dvar_bumcov_amsua_n19_0000.nc4, where
    #   expObsSpaceCombo = "3dvar_bumcov_amsua_n19"
    #
    #   +      obsout_3dvar_bumcov_bumploc_amsua_n19_0000.nc4, where
    #   expObsSpaceCombo = "3dvar_bumcov_bumploc_amsua_n19"
    #
    #   +      obsout_3denvar_abi_g16_0000.nc4, where
    #   expObsSpaceCombo = "3denvar_abi_g16"
    allExpObsFiles = {}
    for fileType, files in allFiles.items():
        for fileName in files:
            # expObsSpaceCombo excludes everything outside the first/final '_'
            expObsSpaceCombo =  '_'.join(fileName.split("_")[1:][:-1])
            if expObsSpaceCombo not in allExpObsFiles:
                allExpObsFiles[expObsSpaceCombo] = {}
            if fileType not in allExpObsFiles[expObsSpaceCombo]:
                allExpObsFiles[expObsSpaceCombo][fileType] = []
            allExpObsFiles[expObsSpaceCombo][fileType].append(fileName)

    # TODO: force crash when there are no ObsSpaceFiles available for an expObsSpaceCombo
    # TODO: force crash when # of ObsDiag or geoval files is > 1 but ~= # of ObsSpaceFiles

    # sort expObsSpaceCombos alphabetically to ensure
    # all processors exhibit identical behavior
    expObsSpaceCombos = list(allExpObsFiles.keys())
    expObsSpaceCombos.sort()

    # sort files by PE
    for expObsSpaceCombo in expObsSpaceCombos:
        expObsFiles = allExpObsFiles[expObsSpaceCombo]
        for fileType, files in expObsFiles.items():
            if len(files) > 0:
                PEs = []
                for fileName in files:
                    PE = ''.join(fileName.split("_")[-1].split(".")[0])
                    PEs.append(int(PE))
                indices = list(range(len(files)))
                indices.sort(key=PEs.__getitem__)
                allExpObsFiles[expObsSpaceCombo][fileType] = \
                    list(map(files.__getitem__, indices))

    # Loop over all experiment+observation combinations
    jj=-1
    for expObsSpaceCombo in expObsSpaceCombos:
        # determine ObsSpace from expObsSpaceCombo string
        # and only process valid ObsSpaceNames
        expt_parts = expObsSpaceCombo.split("_")
        nstr = len(expt_parts)
        ObsSpaceInfo = vu.nullDiagSpaceInfo
        for i in range(0,nstr):
            ObsSpaceName_ = '_'.join(expt_parts[i:nstr])
            ObsSpaceInfo_ = vu.DiagSpaceDict.get( ObsSpaceName_,vu.nullDiagSpaceInfo)
            if ObsSpaceInfo_['process']:
                ObsSpaceName = deepcopy(ObsSpaceName_)
                ObsSpaceInfo = deepcopy(ObsSpaceInfo_)
        if not ObsSpaceInfo['process']: continue
        if ObsSpaceInfo.get('DiagSpaceGrp',vu.model_s) == vu.model_s: continue

        # process valid ObsSpaceNames round-robin on available processors
        jj = jj + 1
        if jj%nproc != myproc: continue

        paru.proc_print(nproc,myproc,this_program+": EXPERIMENT_OBSSPACE =  "+expObsSpaceCombo)

        ## Extract constructor info about this ObsSpace
        ObsSpaceGrp = ObsSpaceInfo['DiagSpaceGrp']
        binGrps = ObsSpaceInfo.get('binGrps',{})

        # compile list of unique variables needed for binning
        uniqBinVarGrps = []
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

        # extract relevant file lists
        expObsFiles = allExpObsFiles[expObsSpaceCombo]
        firstObsSpaceFile = expObsFiles[ObsSpaceFileKey][0]

        # select departure variables with the suffix vu.depbgGroup
        #  from firstObsSpaceFile
        nc = Dataset(firstObsSpaceFile, 'r')
        varlist = nc.variables.keys()
        varlist = [''.join(var.split("@")[:-1]) for var in varlist
                    if (var[-len(vu.depbgGroup):] == vu.depbgGroup)]

        # sort varlist alphabetically
        indices = list(range(len(varlist)))
        if ObsSpaceGrp == vu.radiance_s:
            # Sort by channel number (int) for radiances
            chlist = [int(''.join(varName.split("_")[-1:])) for varName in varlist]
            indices.sort(key=chlist.__getitem__)
        else:
            indices.sort(key=varlist.__getitem__)
        varlist = list(map(varlist.__getitem__, indices))

        # compile list of unique obs and departure variables and
        #  binning variables that are unique to varName
        uniqObsVarGrps = []
        for varName in varlist:
            # Collect obs and departure variables
            uniqObsVarGrps.append(varName+'@'+vu.obsGroup)
            uniqObsVarGrps.append(varName+'@'+vu.depbgGroup)
            uniqObsVarGrps.append(varName+'@'+vu.depanGroup)

            # Collect binning variables that
            # (1) contain vu.vNameStr
            # (2) contain vu.vChanStr
            if ObsSpaceGrp == vu.radiance_s:
                chanStr = ''.join(varName.split("_")[-1:])
            else:
                chanStr = ''

            for ivar, varGrp in enumerate(uniqBinVarGrps):
                binVarName = ''.join(varGrp.split("@")[:-1])
                if (binVarName == vu.vNameStr
                    or vu.vChanStr in binVarName.split("_")):
                    varGrpSub = re.sub(vu.vNameStr,varName,varGrp)
                    varGrpSub = re.sub(vu.vChanStr,chanStr,varGrpSub)
                    uniqBinVarGrps.insert(ivar+1,varGrpSub)

        # combine bin and obs variable lists
        uniqFileVarGrps = deepcopy(uniqBinVarGrps)
        for varGrp in uniqObsVarGrps:
           if varGrp not in uniqFileVarGrps:
               uniqFileVarGrps.append(varGrp)

        paru.proc_print(nproc,myproc,this_program+": Reading all variables from file(s)...")
        #print(uniqFileVarGrps)
        fileVarsVals = readVarsNC(uniqFileVarGrps,expObsFiles)

        # convert pressure from Pa to hPa if needed (kludge)
        for varGrp, vals in fileVarsVals.items():
            if (vu.obsVarPrs in varGrp and np.max(vals) > 10000.0):
                fileVarsVals[varGrp] = np.divide(vals,100.0)

        # Initialize a dictionary to contain all statistical info for this expObsSpaceCombo
        statsDict = {}
        for attribName in su.fileStatAttributes:
            statsDict[attribName] = []
        for statName in su.allFileStats:
            statsDict[statName] = []

        # Collect statistics for each varName in varlist
        paru.proc_print(nproc,myproc,this_program+": Calculating/writing diagnostic stats for:")
        for varName in varlist:
            paru.proc_print(nproc,myproc,this_program+": VARIABLE = "+varName)
            obsName   = varName+'@'+vu.obsGroup
            depbgName = varName+'@'+vu.depbgGroup
            depanName = varName+'@'+vu.depanGroup

            if ObsSpaceGrp == vu.radiance_s:
                chanStr = ''.join(varName.split("_")[-1:])
                dictName = '_'.join(varName.split("_")[:-1])
                varVal = vu.varDictObs.get(dictName,['',dictName])
                varShort = varVal[1]+'_'+chanStr
            else:
                varVal = vu.varDictObs.get(varName,['',varName])
                varShort = varVal[1]
                chanStr = ''
            varUnits = varVal[0]


            # add binning variables that are functions of
            # variables already contained within fileVarsVals:
            # (1) varName@vu.bakGroup
            varGrp = varName+'@'+vu.bakGroup
            if varGrp in uniqBinVarGrps:
                fileVarsVals[varGrp] = np.add(
                    fileVarsVals[depbgName],
                    fileVarsVals[obsName])

            # Collect and write stats for all selectDiagNames
            for diagName in selectDiagNames:
                paru.proc_print(nproc,myproc,this_program+": DIAG = "+diagName)
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
                    paru.proc_print(nproc,myproc,'\n\nERROR: diagName is undefined: '+diagName)
                    os._exit(1)

                if ((diagName == 'omb' or diagName == 'oma')
                   and 'gnssro' in ObsSpaceName):
                    Diagnostic = np.multiply(np.divide(Diagnostic,fileVarsVals[obsName]),100.)

                # Filter and write Diagnostic statistics for each binGrp
                for binGrpKey, binMethods in binGrps.items():
                    grpDict = vu.binGrpDict.get(binGrpKey,{'variables': [vu.miss_s]})
                    if grpDict['variables'][0] == vu.miss_s: continue

                    binVarName = grpDict.get('binVarName',vu.miss_s)
                    if "@" in binVarName:
                        binVar = ''.join(binVarName.split("@")[:-1])
                    else:
                        binVar = binVarName
                    binVarDict = vu.varDictObs.get(binVar,['null',binVar])
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
                                for ivar, varGrp in enumerate(binArgsMod):
                                    varGrpSub = re.sub(vu.vNameStr,varName,varGrp)
                                    varGrpSub = re.sub(vu.vChanStr,chanStr,varGrpSub)
                                    if (diagName == 'oma' or diagName == 'ana'):
                                        varGrpSub = re.sub(vu.qcGroup,vu.qcGroup+qcIterAN,varGrpSub)
                                    else:
                                        varGrpSub = re.sub(vu.qcGroup,vu.qcGroup+qcIterBG,varGrpSub)
                                    binArgsMod[ivar] = varGrpSub

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
                                    paru.proc_print(nproc,myproc,'\n\nERROR: mask is incorrectly defined for '
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

        #Create a new stats file for expObsSpaceCombo
        su.write_stats_nc(expObsSpaceCombo,statsDict)

MAXINT32  = np.int32(1e9)
MAXFLOAT  = np.float32(1.e12)
MAXDOUBLE = np.float64(1.e12)

def readVarsNC(varNames,ncFileDict):
# varNames (string)       - list of variables desired from the obsNCFiles
# ncFileDict (dictionary) - contains a list of NC files for each *FileKey at the top of this script

    obsNCFiles  = ncFileDict.get(ObsSpaceFileKey,[])
    diagNCFiles = ncFileDict.get(ObsDiagFileKey,[])
    geoNCFiles  = ncFileDict.get(geovalFileKey,[])

    obsNCs  = []
    diagNCs = []
    geoNCs  = []

    # Open all NC files as Dataset's
    diagExists = True
    geoExists = True
    for iFile, obsFile in enumerate(obsNCFiles):
        obsNCs.append(Dataset(obsFile, 'r'))
        obsNCs[iFile].set_auto_mask(False)

        if len(diagNCFiles) == len(obsNCFiles) and diagExists:
            diagFile = diagNCFiles[iFile]
            if os.path.exists(diagFile):
                diagNCs.append(Dataset(diagFile, 'r'))
                diagNCs[iFile].set_auto_mask(False)
            else:
                diagExists = False
        else:
            diagExists = False

        if len(geoNCFiles) == len(obsNCFiles) and geoExists:
            geoFile = geoNCFiles[iFile]
            if os.path.exists(geoFile):
                geoNCs.append(Dataset(geoFile, 'r'))
                geoNCs[iFile].set_auto_mask(False)
            else:
                geoExists = False
        else:
            geoExists = False

    # Construct output dictionary
    varsVals = {}
    iNC = 0
    obsNC = obsNCs[iNC]
    for varGrp in varNames:
        varName = ''.join(varGrp.split("@")[:-1])
        grpName = ''.join(varGrp.split("@")[-1])
        if varGrp in obsNC.variables:
            varsVals[varGrp] = np.asarray([])
        elif vu.diagGroup in grpName and diagExists:
            diagNC = diagNCs[iNC]
            if varName in diagNC.variables:
                varsVals[varGrp] = np.asarray([])
        elif vu.geovalGroup in grpName and geoExists:
            geoNC = geoNCs[iNC]
            if varName in geoNC.variables:
                varsVals[varGrp] = np.asarray([])

    # Read variables across all Dataset's
    diagLev = 0
    geoLev = 0

    for iNC, obsNC in enumerate(obsNCs):
        for varGrp in varNames:
            varName = ''.join(varGrp.split("@")[:-1])
            grpName = ''.join(varGrp.split("@")[-1])
            varFound = True
            if varGrp in obsNC.variables:
                varsVals[varGrp] = np.append( varsVals[varGrp], obsNC.variables[varGrp] )
                dtype = obsNC.variables[varGrp].datatype
            elif vu.diagGroup in grpName and diagExists:
                diagNC = diagNCs[iNC]
                if varName in diagNC.variables:
                    varsVals[varGrp] = np.append( varsVals[varGrp],
                            np.transpose(np.asarray(diagNC.variables[varName][:]))[diagLev][:] )
                    dtype = diagNC.variables[varName].datatype
            elif vu.geovalGroup in grpName and geoExists:
                geoNC = geoNCs[iNC]
                if varName in geoNC.variables:
                    varsVals[varGrp] = np.append( varsVals[varGrp],
                            np.transpose(np.asarray(geoNC.variables[varName][:]))[geoLev][:] )
                    dtype = diagNC.variables[varName].datatype
            else:
                varFound = False

            if varFound:
                if 'int32' in dtype.name:
                    missing = np.greater(np.abs(varsVals[varGrp]), MAXINT32)
                elif 'float32' in dtype.name:
                    missing = np.greater(np.abs(varsVals[varGrp]), MAXFLOAT)
                elif 'float64' in dtype.name:
                    missing = np.greater(np.abs(varsVals[varGrp]), MAXDOUBLE)
                else:
                    missing = np.full_like(varsVals[varGrp],False,dtype=bool)

                varsVals[varGrp][missing] = np.NaN

    return varsVals

#================================
#================================

def main():
    nproc, myproc = paru.par_args(sys.argv[:])
    writediagstats_obsspace(nproc, myproc)

if __name__ == '__main__': main()
