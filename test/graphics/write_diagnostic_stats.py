from netCDF4 import Dataset
import os, sys
#import getopt
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch
import math
import plot_utils as pu

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See plot_utils.par_args for more information.

diagNames = ['omb','oma','obs','bak','ana']

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
qcbg_var  = 'EffectiveQC0' #EffectiveQCi, where i is the iteration for depbg_var
qcan_var  = 'EffectiveQC2' #EffectiveQCi, where i is the iteration for depan_var

def write_diag_stats(nproc, myproc):
#  nproc - total number of processors
#  myproc - processor rank, starting at 0
#  Note: these arguments are used when an external program has multiple processors
#        available.  The default values are 1 and 0, respectively.

    # Assign processors round-robin to each ObsSpace name
    ObsSpaceDict = {}
    for ii, (key,baseval) in enumerate(pu.ObsSpaceDict_base.items()):
        val = deepcopy(baseval)
        if ii%nproc != myproc: continue
        ObsSpaceDict[key] = val

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
    for j, file_name in enumerate(obsoutfiles):
        # exob_group_name excludes everything outside the first/final '_'
        exob_group_name =  '_'.join(file_name.split("_")[1:][:-1])
        for i, exob_group in enumerate(exob_groups):
            if exob_group_name in exob_group:
                # If exob_group with exob_group_name exists, add new file to it
                update_group = exob_group
                update_group.append(file_name)
                exob_groups[i][:] = update_group
                break
            elif i == len(exob_groups)-1:
                # If group with exob_group_name does not exist, add one
                new_group = [exob_group_name]
                new_group.append(file_name)
                exob_groups.append(new_group)
                break
    exob_groups = exob_groups[1:][:]

    # Loop over unique experiment-ObsSpace groups
    for exob_group in exob_groups:
        expt_obs = exob_group[0]
        # Determine ObsSpace from expt_ObsSpace string
        expt_parts = expt_obs.split("_")
        nstr = len(expt_parts)
        ObsSpaceName = pu.miss_s
        for i in range(0,nstr):
            ObsSpaceName_ = '_'.join(expt_parts[i:nstr])
            ObsSpaceInfo_ = ObsSpaceDict.get( ObsSpaceName_,[pu.miss_s, 0, pu.nullBinKeys])
            if ObsSpaceInfo_[0] != pu.miss_s and ObsSpaceInfo_[1] > 0:
                ObsSpaceName = ObsSpaceName_
                ObsSpaceInfo = ObsSpaceInfo_

        if ObsSpaceName == pu.miss_s:
           continue

        pu.proc_print(nproc,myproc,"write_diag_stats: EXPERIMENT_OBSSPACE =  "+expt_obs)

        ObsSpaceGrp = ObsSpaceInfo[0]
        binGrps = ObsSpaceInfo[2]

        # compile list of unique variables needed for binning
        uniqBinVars = []
        for binGrp in binGrps:
            binGrpKey = binGrp[0]
            if binGrpKey == pu.miss_s: continue
            grpDict = pu.binGrpDict.get(binGrpKey,{'variables': [pu.miss_s]})
            binVars = grpDict['variables']
            if binVars[0] == pu.miss_s: continue
            for binVar in binVars:
                uniqBinVars.append(binVar)
        uniqBinVars = pu.uniqueMembers(uniqBinVars)

        #Extract binning variables (e.g., vertical coordinate, latitude)
        # in loop over exob_group, 
        # excluding category in exob_group[0]
        fileBinVarsVals = [np.asarray([])]*len(uniqBinVars)
        for file_name in exob_group[1:]:
            nc = Dataset(file_name, 'r')
            nc.set_auto_mask(False)

            for ivar, binVarGrp in enumerate(uniqBinVars):
                if binVarGrp != pu.miss_s and binVarGrp != '':
                    binVar = ''.join(binVarGrp.split("@")[:-1])

                    tmp = nc.variables[binVarGrp]
                    #Convert air_pressure from Pa to hPa if needed
                    if binVar == 'air_pressure' and np.max(tmp) > 10000.0:
                        tmp = np.divide(tmp,100.0)
                    fileBinVarsVals[ivar] = np.append( fileBinVarsVals[ivar], tmp )


        # Select variables with the suffix depbg_var (required for OMB)
        nc = Dataset(exob_group[1], 'r')
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

        # Loop over variables with omb suffix
        for ivar, depbg in enumerate(bglist):
            varName = ''.join(depbg.split("@")[:-1])
            pu.proc_print(nproc,myproc,"write_diag_stats: VARIABLE = "+varName)
            obs=''.join(depbg.split("@")[:-1])+'@'+obs_var
            depan=''.join(depbg.split("@")[:-1])+'@'+depan_var
            qcb = ''.join(depbg.split("@")[:-1])+'@'+qcbg_var
            qca = ''.join(depbg.split("@")[:-1])+'@'+qcan_var
            #print(nproc,myproc,"obs="+obs+"depbg="+depbg+"depan="+depan)

            obsnc = np.asarray([])
            ombnc = np.asarray([])
            omanc = np.asarray([])
            qcbnc = np.asarray([])
            qcanc = np.asarray([])

            # Build up arrays in loop over exob_group, 
            # excluding category in exob_group[0]
            for file_name in exob_group[1:]:
                nc = Dataset(file_name, 'r')
                nc.set_auto_mask(False)
                #file_rank = file_name[-(4+npedigits):-4]

                obsnc = np.append( obsnc, nc.variables[obs] )
                ombnc = np.append( ombnc, np.negative( nc.variables[depbg] ) ) # omb = (-) depbg
                qcbnc = np.append( qcbnc, nc.variables[qcb]  )
                if depan in nc.variables:
                    omanc = np.append( omanc, np.negative( nc.variables[depan] ) ) # oma = (-) depan
                if qca in nc.variables:
                    qcanc  = np.append( qcanc,  nc.variables[qca]  )

            #@EffectiveQC, 1: missing; 0: good; 10: rejected by first-guess check
            #keep data @EffectiveQC=0
            obsnc[np.less(obsnc,-1.0e+15)] = np.NaN
            obsnc[qcbnc != 0]              = np.NaN

            ombnc[np.less(ombnc,-1.0e+15)] = np.NaN
            ombnc[np.isnan(obsnc)]         = np.NaN
            ombnc[qcbnc != 0]              = np.NaN

            omanc[np.less(omanc,-1.0e+15)] = np.NaN
            omanc[np.isnan(obsnc)]         = np.NaN
            omanc[qcanc != 0]              = np.NaN

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
                if diagName == 'omb':
                    Diagnostic = deepcopy(ombnc)
                elif diagName == 'oma':
                    Diagnostic = deepcopy(omanc)
                elif diagName == 'obs':
                    Diagnostic = deepcopy(obsnc)
                elif diagName == 'bak':
                    Diagnostic = np.subtract(obsnc, ombnc)
                elif diagName == 'ana':
                    Diagnostic = np.subtract(obsnc, omanc)
                else:
                    pu.proc_print(nproc,myproc,'\n\nERROR: diagName is undefined: '+diagName)
                    os._exit(1)

                if (diagName == 'omb' or diagName == 'oma') \
                   and  ObsSpaceName[:6] == 'gnssro':
                    Diagnostic = (Diagnostic / obsnc) * 100.

                # Write un-binned Diagnostic
                newFile = write_stats(expt_obs, ObsSpaceGrp, varShortName, varUnits, diagName, \
                                      pu.miss_s, pu.allBins, pu.miss_s, Diagnostic, newFile)

                # Write binned Diagnostic
                for binGrp in binGrps:
                    binGrpKey = binGrp[0]
                    if binGrpKey == pu.miss_s: continue
                    grpDict = pu.binGrpDict.get(binGrpKey,{'variables': [pu.miss_s]})
                    binVars = grpDict['variables']
                    if binVars[0] == pu.miss_s: continue

                    selectKeys = binGrp[1]
                    binVarsVals = [fileBinVarsVals[uniqBinVars.index(binVar)] for binVar in binVars]

                    nlocs = len(binVarsVals[0])
                    for ikey, selectKey in enumerate(selectKeys):
                        keyDesc = grpDict.get(selectKey,pu.nullBinDesc)
                        binVals = keyDesc['values']
                        if binVals[0] == pu.miss_s: continue

                        binFilters = keyDesc['filters']
                        for ibin, binVal in enumerate(binVals):

                            if len(binVars) == 1 and pu.isfloat(binVal):
                                binVar = ''.join(binVars[0].split("@")[:-1])
                                binVarDict = pu.varDict.get(binVar,['null',binVar])
                                binVar = binVarDict[1]
                                binUnits = binVarDict[0]
                            else:
                                binVar = pu.miss_s
                                binUnits = pu.miss_s

                            # filter/remove data that is outside this bin
                            blackList = np.full((nlocs),False)
                            for binFilter in binFilters:
                                binFunc = binFilter['where']
                                binArgs = binFilter['args']
                                binBound = (binFilter['bounds'])[ibin]
                                binArgsVals = [binVarsVals[i] for i in binArgs]
                                blackList = np.logical_or(blackList, binFunc(binArgsVals,binBound))

                            binnedDiagnostic = deepcopy(Diagnostic)
                            binnedDiagnostic[blackList] = np.NaN
                            
                            newFile = write_stats(expt_obs, ObsSpaceGrp, varShortName, varUnits, diagName, \
                                                  binVar, binVal, binUnits, binnedDiagnostic, newFile)


def write_stats(ObsSpace, ObsSpaceGrp, varName, varUnits, diagName, \
                binVar, binVal, binUnits, array_f, newFile):
# ObsSpace (string)      - ObsSpace name (e.g., sonde)
# ObsSpaceGrp (string)   - group of the ObsSpace (e.g., profile, radar, radiance)
# varName (string)       - observed variable associated with the diagnostic
# varUnits (string)      - units of the variable
# diagName (string)      - name of the diagnostic (e.g., 'omb')
# binVar (string)        - binning variable (e.g., altitude)
# binVal (string)        - name of this bin
# binUnits (string)      - units of binning variable
# array_f (float(:))     - 1-d array of float values for which statistics should be calculated
# newFile (logical)      - whether a new file should be created

    statsFile = 'stats_'+ObsSpace+'.txt'
    if newFile and os.path.exists(statsFile):
        os.remove(statsFile)

    fp = open(statsFile, 'a')

    #Only include non-NaN values in statistics
    STATS = {}
    STATS['Count']  = len(array_f)-np.isnan(array_f).sum()
    if STATS['Count'] > 0:
        STATS['Mean'] = np.nanmean(array_f)
        STATS['RMS']  = np.sqrt(np.nanmean(array_f**2))
        STATS['STD']  = np.nanstd(array_f)
        STATS['Min']  = np.nanmin(array_f)
        STATS['Max']  = np.nanmax(array_f)
    else:
        STATS['Mean'] = np.NaN
        STATS['RMS']  = np.NaN
        STATS['STD']  = np.NaN
        STATS['Min']  = np.NaN
        STATS['Max']  = np.NaN

    #TODO: use binary/netcdf or other non-ASCII format
    line = ObsSpaceGrp          \
           +pu.csvSEP+varName   \
           +pu.csvSEP+varUnits  \
           +pu.csvSEP+diagName  \
           +pu.csvSEP+binVar    \
           +pu.csvSEP+binVal    \
           +pu.csvSEP+binUnits

    for statName in pu.allFileStats:
        if statName == "Count":
            line = line+pu.csvSEP+"{:d}".format(STATS[statName])
        else:
            line = line+pu.csvSEP+"{:.8e}".format(STATS[statName])

    fp.write( line+"\n")

    fp.close()

    if os.path.exists(statsFile):
        return False
    else:
        return True

#================================
#================================

def main():
    nproc, myproc = pu.par_args(sys.argv[:])
    write_diag_stats(nproc, myproc)

if __name__ == '__main__': main()
