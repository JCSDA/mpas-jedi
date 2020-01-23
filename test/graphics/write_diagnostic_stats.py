from copy import deepcopy
import fnmatch
#import math
import matplotlib.pyplot as plt
import matplotlib.axes as maxes
from netCDF4 import Dataset
import numpy as np
import os
import plot_utils as pu
import re
import sys

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See plot_utils.par_args for more information.

diagNames = pu.allDiags

# Diagnostic oma/omb files located in diagdir
diagdir    = '../Data/'
diagprefix = 'obsout_'
diagsuffix = '_*.nc4'
# * in diagsuffix is 4-digit processor rank [0000 to XXXX]
#npedigits = 4

# OOPS suffixes for required nc variables
# observations
obs_var = 'ObsValue'

# departures
depbg_var = 'depbg'
depan_var = 'depan'

#quality control
qcIterBG='0'
NOUTER=os.getenv('NOUTER',1) #set equal to number of outer iterations
qcIterAN=str(NOUTER)

qcbg_var  = 'EffectiveQC0' #EffectiveQCi, where i is the iteration for depbg_var
qcan_var  = 'EffectiveQC2' #EffectiveQCi, where i is the iteration for depan_var

def write_diag_stats(nproc, myproc):
#  nproc - total number of processors
#  myproc - processor rank, starting at 0
#  Note: these arguments are used when an external program has multiple processors
#        available.  The default values are 1 and 0, respectively.

    # Assign processors round-robin to each ObsSpace name
    ObsSpaceDict = {}
    for ii, (key,baseval) in enumerate(pu.DiagSpaceDict.items()):
        if ii%nproc != myproc or baseval['DiagSpaceGrp'] == pu.model_s: continue
        ObsSpaceDict[key] = deepcopy(baseval)

    obsoutfiles = []
    for files in os.listdir(diagdir):
        #print(files)
        if fnmatch.fnmatch(files, diagprefix+'*'+diagsuffix):   # 1tile
            obsoutfiles.append(diagdir+files)
    #print(obsoutfiles)

    # Group files by experiment-ObsSpace combination
    #  (e.g., 3dvar_aircraft), where each group 
    #  contains files from all PE's
    exob_groups = [[]]
    for j, fileName in enumerate(obsoutfiles):
        # exob_group_name excludes everything outside the first/final '_'
        exob_group_name =  '_'.join(fileName.split("_")[1:][:-1])
        for i, exob_group in enumerate(exob_groups):
            if exob_group_name in exob_group:
                # If exob_group with exob_group_name exists, add new file to it
                update_group = exob_group
                update_group.append(fileName)
                exob_groups[i][:] = update_group
                break
            elif i == len(exob_groups)-1:
                # If group with exob_group_name does not exist, add one
                new_group = [exob_group_name]
                new_group.append(fileName)
                exob_groups.append(new_group)
                break
    exob_groups = exob_groups[1:][:]

    # Loop over unique experiment-ObsSpace groups
    for exob_group in exob_groups:
        expt_obs = exob_group[0]
        ObsSpaceFiles = exob_group[1:]
        # Determine ObsSpace from expt_ObsSpace string
        expt_parts = expt_obs.split("_")
        nstr = len(expt_parts)
        ObsSpaceInfo = pu.nullDiagSpaceInfo
        for i in range(0,nstr):
            ObsSpaceName_ = '_'.join(expt_parts[i:nstr])
            ObsSpaceInfo_ = ObsSpaceDict.get( ObsSpaceName_,pu.nullDiagSpaceInfo)
            if ObsSpaceInfo_['process']:
                ObsSpaceName = ObsSpaceName_
                ObsSpaceInfo = ObsSpaceInfo_

        if not ObsSpaceInfo['process']: continue

        pu.proc_print(nproc,myproc,"write_diag_stats: EXPERIMENT_OBSSPACE =  "+expt_obs)

        ObsSpaceGrp = ObsSpaceInfo['DiagSpaceGrp']
        binGrps = ObsSpaceInfo['binGrps']

        # compile list of unique variables needed for binning
        uniqBinVarGrps = []
        for binGrp in binGrps:
            binGrpKey = binGrp[0]
            if binGrpKey == pu.miss_s: continue
            grpDict = pu.binGrpDict.get(binGrpKey,{'variables': [pu.miss_s]})
            binVarGrps = grpDict['variables']
            for binVarGrp in binVarGrps:
                if binVarGrp == pu.miss_s: continue
                uniqBinVarGrps.append(binVarGrp)
        uniqBinVarGrps = pu.uniqueMembers(uniqBinVarGrps)

        #Add appropriate QC interation numbers
        for ivar, varGrp in enumerate(uniqBinVarGrps):
            if varGrp==pu.varQC:
                uniqBinVarGrps[ivar] = varGrp+qcIterBG
                uniqBinVarGrps.append(varGrp+qcIterAN)

        uniqBinVars = []
        uniqBinGrps = []
        for varGrp in uniqBinVarGrps:
            uniqBinVars.append(''.join(varGrp.split("@")[:-1]))
            uniqBinGrps.append(''.join(varGrp.split("@")[-1]))

        #Extract binning variables (e.g., vertical coordinate, latitude) 
        fileBinVarsVals = collectVarNC(uniqBinVarGrps,ObsSpaceFiles)
        #Convert air_pressure from Pa to hPa if needed
        for ivar, varName in enumerate(uniqBinVars):
            if varName == 'air_pressure' and np.max(fileBinVarsVals[ivar]) > 10000.0:
                fileBinVarsVals[ivar] = np.divide(fileBinVarsVals[ivar],100.0)

        # Select variables with the suffix depbg_var from the first file
        nc = Dataset(ObsSpaceFiles[0], 'r')
        varlist = nc.variables.keys()
        bglist = [var for var in varlist if (var[-len(depbg_var):] == depbg_var)]

        # Sort bglist alphabetically
        indices = list(range(len(bglist)))
        if ObsSpaceGrp == pu.radiance_s:
            # Sort by channel number (int) for radiances
            chlist = []
            for bgname in bglist:
                var = ''.join(bgname.split("@")[:-1])
                chlist.append(int(''.join(var.split("_")[-1:])))
            indices.sort(key=chlist.__getitem__)
        else:
            indices.sort(key=bglist.__getitem__)
        bglist = list(map(bglist.__getitem__, indices))

        #Create a new stats file for exob_group
        newFile = True

        statsDict = {}
        for attribName in pu.fileStatAttributes:
            statsDict[attribName] = []
        for statName in pu.allFileStats:
            statsDict[statName] = []

        # Loop over variables with omb suffix
        for ivar, depbgName in enumerate(bglist):
            varName = ''.join(depbgName.split("@")[:-1])
            pu.proc_print(nproc,myproc,"write_diag_stats: VARIABLE = "+varName)
            obsName=''.join(depbgName.split("@")[:-1])+'@'+obs_var
            depanName=''.join(depbgName.split("@")[:-1])+'@'+depan_var
            #print(nproc,myproc,"obs="+obsName+"depbgName="+depbgName+"depan="+depanName)

            tmp = collectVarNC([obsName, depbgName, depanName],ObsSpaceFiles)
            obsVals = tmp[0]
            ombVals = np.negative(tmp[1]) # omb = (-) depbg
            omaVals = np.negative(tmp[2]) # oma = (-) depan

            #gather non-QC varNam@* binning variables for this varName
            for ivar, varGrp in enumerate(uniqBinVarGrps):
                if uniqBinVars[ivar]==pu.vNameStr:
                    varGrpSub = re.sub(pu.vNameStr,varName,varGrp)
                    tmp = collectVarNC([varGrpSub],ObsSpaceFiles)
                    fileBinVarsVals[ivar] = tmp[0]

            if ObsSpaceGrp == pu.radiance_s:
                ch = ''.join(varName.split("_")[-1:])
                dictName = '_'.join(varName.split("_")[:-1])
                varVal = pu.varDict.get(dictName,['',dictName])
                varShortName = varVal[1]+'_'+ch
            else:
                varVal = pu.varDict.get(varName,['',varName])
                varShortName = varVal[1]
            varUnits = varVal[0]

            # Collect and write stats for all diagNames
            for diagName in diagNames:
                pu.proc_print(nproc,myproc,"write_diag_stats: DIAG = "+diagName)
                if diagName == 'omb':
                    Diagnostic = deepcopy(ombVals)
                elif diagName == 'oma':
                    Diagnostic = deepcopy(omaVals)
                elif diagName == 'obs':
                    Diagnostic = deepcopy(obsVals)
                elif diagName == 'bak':
                    Diagnostic = np.subtract(obsVals, ombVals)
                elif diagName == 'ana':
                    Diagnostic = np.subtract(obsVals, omaVals)
                else:
                    pu.proc_print(nproc,myproc,'\n\nERROR: diagName is undefined: '+diagName)
                    os._exit(1)

                if (diagName == 'omb' or diagName == 'oma') \
                   and  ObsSpaceName[:6] == 'gnssro':
                    Diagnostic = (Diagnostic / obsVals) * 100.

                # Filter and write Diagnostic statistics for each binGrp
                for binGrp in binGrps:
                    binGrpKey = binGrp[0]
                    if binGrpKey == pu.miss_s: continue
                    grpDict = pu.binGrpDict.get(binGrpKey,{'variables': [pu.miss_s]})
                    binVars = grpDict['variables']
                    if binVars[0] == pu.miss_s: continue
                    nonQCBinVars = [var for var in binVars if var!=pu.varQC]
                    binKeys = binGrp[1]

                    # collect only the variables associated with this binGrp and diagName
                    iterBinVars = deepcopy(binVars)
                    if (diagName == 'oma' or diagName == 'ana'):
                        #use AN iteration for varNam@EffectiveQC
                        for ivar, varGrp in enumerate(binVars):
                            if varGrp==pu.varQC:
                                iterBinVars[ivar] = varGrp+qcIterAN
                    else:
                        #use BG iteration for varNam@EffectiveQC
                        for ivar, varGrp in enumerate(binVars):
                            if varGrp==pu.varQC:
                                iterBinVars[ivar] = varGrp+qcIterBG
                    binVarsVals = [fileBinVarsVals[uniqBinVarGrps.index(binVar)] for binVar in iterBinVars]

                    nlocs = len(binVarsVals[0])
                    for ikey, binKey in enumerate(binKeys):
                        keyDesc = grpDict.get(binKey,pu.nullBinDesc)
                        binVals = keyDesc['labels']
                        if binVals[0] == pu.miss_s: continue

                        binFilters = keyDesc['filters']
                        if len(nonQCBinVars) == 1 and pu.isfloat(binVals[0]):
                            binVar = ''.join(nonQCBinVars[0].split("@")[:-1])
                            binVarDict = pu.varDict.get(binVar,['null',binVar])
                            binVar = binVarDict[1]
                            binUnits = binVarDict[0]
                        elif len(binVars) == 1 and binVars[0]==pu.varQC:
                            binVar = 'QCFlag'
                            binUnits = pu.miss_s
                        else:
                            binVar = pu.miss_s
                            binUnits = pu.miss_s

                        for ibin, binVal in enumerate(binVals):
                            # filter data that meets bin criteria
                            binnedDiagnostic = deepcopy(Diagnostic)
                            for binFilter in binFilters:
                                applyToDiags = binFilter.get('apply_to',pu.allDiags)
                                if not (diagName in applyToDiags): continue
                                binFunc  = binFilter['where']
                                binArgs  = binFilter['args']
                                binBound = (binFilter['bounds'])[ibin]
                                binArgsVals = [binVarsVals[binVars.index(arg)] for arg in binArgs]
                                mask = binFunc(binArgsVals,binBound)
                                maskValue = binFilter.get('mask_value',np.NaN)
                                if len(mask) == len(binnedDiagnostic):
                                    binnedDiagnostic[mask] = maskValue
                                else:
                                    pu.proc_print(nproc,myproc,'\n\nERROR: mask is incorrectly defined for '\
                                                  +diagName+" "+binGrpKey+" "+binKey+" "+binVal)
                                    os._exit(1)
                            #END binFilters LOOP

                            statsDict['DiagSpaceGrp'].append(ObsSpaceGrp)
                            statsDict['varName'].append(varShortName)
                            statsDict['varUnits'].append(varUnits)
                            statsDict['diagName'].append(diagName)
                            statsDict['binVar'].append(binVar)
                            statsDict['binVal'].append(binVal)
                            statsDict['binUnits'].append(binUnits)

                            statsVal = pu.calcStats(binnedDiagnostic)
                            for statName in pu.allFileStats:
                                statsDict[statName].append(statsVal[statName])

                        #END binVals LOOP
                    #END binKeys LOOP
                #END binGrps LOOP
            #END diagNames LOOP
        #END bglist LOOP

        pu.write_stats_nc(expt_obs,statsDict)
        

def collectVarNC(varNames,ncFiles):
# varNames (string) - list of variables desired from the ncFiles
# ncFiles (string)  - list of netcdf files to be read (e.g., from multiple processor/experiment outputs)

    varsVals = [np.asarray([])]*len(varNames)

    # Build up array in loop over ncFiles
    for fileName in ncFiles:
        nc = Dataset(fileName, 'r')
        nc.set_auto_mask(False)
        for ivar, varName in enumerate(varNames):
            if varName in nc.variables:
                varsVals[ivar] = np.append( varsVals[ivar], nc.variables[varName] )
    return varsVals

#================================
#================================

def main():
    nproc, myproc = pu.par_args(sys.argv[:])
    write_diag_stats(nproc, myproc)

if __name__ == '__main__': main()
