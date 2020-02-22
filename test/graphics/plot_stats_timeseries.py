import basic_plot_functions as bpf
from copy import deepcopy
import datetime as dt
import glob
import numpy as np
import pandas as pd
#from pandas.plotting import register_matplotlib_converters
#register_matplotlib_converters()
import par_utils as paru
import plot_utils as pu
import re
import os
import stat_utils as su
import sys
import var_utils as vu

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See par_utils.par_args for more information.

#Select the type of plots
# options: 'SNAPSHOTCY', 'AGGREGATEFC','SNAPSHOTCY2D', 'AGGREGATEFC2D', 'SNAPSHOTFC'
# SNAPSHOTCY - creates a timeseries figure between firstCycleDTime and lastCycleDTime
#               for each forecast length between fcTDeltaFirst and fcTDeltaLast
#             - x-axis: cycle initial time
#             -    line: per experiment
#             - subplot: per diag space variable
#             -    file: combination of FC lead time, statistic, and bin (if applicable)
# AGGREGATEFC - creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing
#               aggregated statistics for the period between firstCycleDTime and lastCycleDTime
#             - x-axis: forecast duration
#             -    line: per experiment
#             - subplot: per diag space variable
#             -    file: combination of statistic and bin
# SNAPSHOTCY2D/AGGREGATEFC2D creates contour maps similar to above with vertical info on y-axis
#             - only applicable to binned diagnostics (e.g., vertical dimension, latitude)
#             - subplot: column by experiment, row by diag space variable
#             SNAPSHOTCY2D
#             -    file: combination of FC lead time and statistic
#             AGGREGATEFC2D
#             -    file: per statistic
# AGGREGATEFCPROFILE similar to AGGREGATEFC2D, except
#             - each vertical column of data is plotted as a profile on
#               a separate set of axes
#             - therefore this is a valid plot even for a single forecast length (omb)
#             -    line: per experiment
#             - subplot: column by lead time, row by diag space variable
#             -    file: per statistic
#             - MAX_FC_SUBFIGS determines number of FC lead times to include
MAX_FC_SUBFIGS = 3
# SNAPSHOTFC  similar to SNAPSHOTCY, except
#             -    line: per FC lead time
#             - subplot: per diag space variable
#             -    file: combination of experiment, statistic, and bin
#             - MAX_FC_LINES determines number of FC lead times to include
MAX_FC_LINES = 5
# SNAPSHOTCY_LAT1 similar to SNAPSHOTCY, except
#             -    line: named latitude bins
#             - subplot: column by experiment, row by diag space variable
#             -    file: combination of statistic and forecast length
#             - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_LAT2 similar to SNAPSHOTCY, except
#             -    line: per experiment
#             - subplot: column by LAT bin, row by diag space variable
#             -    file: combination of statistic and forecast length
#             - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_QC similar to SNAPSHOTCY, except
#             - only plots 'Count' statistic
#             -    line: named QC bins
#             - subplot: column by experiment, row by diag space variable
#             -    file: forecast length
# AGGREGATEFC_DIFFCI/AGGREGATEFCPROFILE_DIFFCI similar to their per experiment counterparts
#             - shows difference between experiment(s) and control
#             - control is selected using cntrlExpInd
cntrlExpInd = 0
#             - statistics are selected with bootStrapStats
bootStrapStats = []
for x in su.sampleableAggStats:
    if x != 'Count': bootStrapStats.append(x)

#TODO: AGGREGATEFC_LAT1/2 named latitude figures w/ FC lead time x-axis

# NOTE: all non-PROFILE *FC* figures require non-zero forecast length
# NOTE: all *CY* figures require > 1 analysis cycle


#Select the stats for plotting
# options: see su.allFileStats
statNames = ['Count','Mean','RMS','STD']


#Select the variable for 2D figures
# options: 'P','alt','lat'
# NOTE: radiances can only be binned by 'lat' so far
binVars2D = {
    'P':       {'profileFunc': bpf.plotProfile}
  , 'alt':     {'profileFunc': bpf.plotProfile}
  , 'lat':     {'profileFunc': bpf.plotProfile}
  , 'zen':     {'profileFunc': bpf.plotSeries}
  , 'cldfrac': {'profileFunc': bpf.plotSeries}
  , 'Ca':      {'profileFunc': bpf.plotSeries}
}

user = os.getenv('USER','jban')

singleFCLen = "singleFCLen"
multiFCLen = "multiFCLen"

plotGroup = singleFCLen
#plotGroup = multiFCLen

figureFileType = 'png'
#figureFileType = 'pdf'

# multiFCLen ASCII statistics file example for cycling run (on cheyenne):
#statFile = '/glade/scratch/user/pandac/DA_3dvar/2018041500/0/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.nc'
#            |                         |        |          | |                      |             |         |
#                       ^                   ^        ^      ^           ^                  ^           ^
#                  expDirectory          expName  cyDTime fcTDelta statsFilePrefix     DAMethod   DiagSpaceName

# singleFCLen ASCII statistics file example for cycling run (on cheyenne):
#statFile = '/glade/scratch/user/pandac/DA_3dvar/2018041500/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.nc'
#            |                         |        |          |                      |             |         |
#                       ^                   ^        ^                ^                  ^           ^
#                  expDirectory          expName  cyDTime      statsFilePrefix     DAMethod   DiagSpaceName


plotTypes = []

if plotGroup == singleFCLen:
    #Select the diagnostics for plotting
    # options: 'omb','oma','obs','bak','ana'
    diagNames = ['omb','oma']
    #TODO (maybe): multiple diagNames on same subplot
    #diagGroups_to_plot = [['omb'],['omb','oma'],['obs','bak','ana']]

    user = 'jban'

    expLongNames = [
                 'conv60it-2node/DADIAG',
                 'amua60it-2node/DADIAG',
                 ]

    expNames = [
                 'conv',
                 'conv+amsua',
                 ]

    DAMethods = [
                 '3denvar_bumploc',
                 '3denvar_bumploc',
                 ]

    #First and Last CYCLE dates
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,5,14,0,0,0)
    cyTimeInc  = dt.timedelta(hours=6)

    #First and Last FORECAST durations
    fcTDeltaFirst = dt.timedelta(days=0)
    fcTDeltaLast  = dt.timedelta(days=0)
    fcTimeInc  = dt.timedelta(hours=24)

    ## plotTypes for considering single forecast length
    ## ------------------------------------------------------------------------
    plotTypes.append('SNAPSHOTCY')
    plotTypes.append('SNAPSHOTCY_LAT2')
    plotTypes.append('SNAPSHOTCY_QC')
    # plotTypes.append('SNAPSHOTCY2D')


if plotGroup == multiFCLen:
    #Select the diagnostics for plotting
    # options: 'omb','obs','bak'
    diagNames = ['omb']

    user = 'jban'

    expLongNames = [
                 'conv60it-2node/OMF',
                 'amua60it-2node/OMF',
                ]

    expNames = [
                 'conv',
                 'conv+amsua',
                ]

    DAMethods = [
                 '3dvar',
                 '3dvar',
                ]

    #First and Last CYCLE dates
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,5,4,0,0,0)
    cyTimeInc  = dt.timedelta(hours=24)

    #First and Last FORECAST durations
    fcTDeltaFirst = dt.timedelta(days=0)
    fcTDeltaLast  = dt.timedelta(days=10)
    fcTimeInc  = dt.timedelta(hours=24)
    #TODO: define FC directory names in terms of seconds or d_hh-mm-ss
    #      in order to allow for increments < 1 day --> modify workflow scripts

    # plotTypes for considering multiple forecast lengths
    # --------------------------------------------------------------------------
    plotTypes.append('AGGREGATEFC')
    plotTypes.append('SNAPSHOTFC')
    if len(expNames) > 1:
        plotTypes.append('AGGREGATEFC_DIFFCI')
    # plotTypes.append('AGGREGATEFC2D')



nExp  = len(expNames)

#AGGREGATEFCPROFILE* figures work for any forecast duration
plotTypes.append('AGGREGATEFCPROFILE')
if nExp > 1:
    plotTypes.append('AGGREGATEFCPROFILE_DIFFCI')

    cntrlName = expNames[min([cntrlExpInd,nExp-1])]
    noncntrlExpNames = [x for x in expNames if x != cntrlName]
    print('\nControl Experiment: '+cntrlName)
    print('\nNon-control Experiment(s): ')
    print(noncntrlExpNames)

expDirectory = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac/')

statsFilePrefix = 'diagnostic_stats/'+su.statsFilePrefix

#plot settings
interiorLabels = True

sciticks = []
for statName in statNames:
    if statName == 'Count': sciticks.append(True)
    else: sciticks.append(False)

#
# primary script to be called by main()
#

def plot_stats_timeseries(nproc=1, myproc=0):
#  nproc - total number of processors
#  myproc - processor rank, starting at 0
#  Note: these arguments are used when the calling program has multiple processors available.

    # Assign processors round-robin to each DiagSpace name
    DiagSpaceDict = {}
    jj=-1
    for ii, (key,baseval) in enumerate(vu.DiagSpaceDict.items()):
        if not baseval['process']: continue
        jj = jj + 1
        if jj%nproc != myproc: continue
        DiagSpaceDict[key] = deepcopy(baseval)

    if fcTDeltaFirst == fcTDeltaLast and firstCycleDTime == lastCycleDTime:
        print("\n\n===============================================================================")
        print("===============================================================================")
        print("\nWARNING: Only PROFILE plotTypes can be generated without a time difference")
        print("         (either forecast or multiple cycles)\n")
        print("===============================================================================")
        print("===============================================================================\n\n")

    fcTDeltas_dir = []
    fcTDeltas = []
    dumTimeDelta = fcTDeltaFirst
    while dumTimeDelta <= fcTDeltaLast:
        #TODO: define FC directory names with in terms of seconds or d_hh-mm-ss
        #fc_date_str = str(dumTimeDelta.total_seconds())
        #fc_date_str = dumTimeDelta.__str__() #[d days, h:mm:ss] --> need to parse
        fc_date_str = str(dumTimeDelta.days)
        fcTDeltas_dir.append(fc_date_str)
        fcTDeltas.append(dumTimeDelta)
        dumTimeDelta = dumTimeDelta + fcTimeInc
    nFC = len(fcTDeltas)

    cyDTimes_dir = []
    cyDTimes = []
    dumDateTime = firstCycleDTime
    while dumDateTime <= lastCycleDTime:
        cy_date_str = "{:04d}".format(dumDateTime.year)  \
                    + "{:02d}".format(dumDateTime.month) \
                    + "{:02d}".format(dumDateTime.day)   \
                    + "{:02d}".format(dumDateTime.hour)
        cyDTimes_dir.append(cy_date_str)
        cyDTimes.append(dumDateTime)
        dumDateTime = dumDateTime + cyTimeInc
    nCY = len(cyDTimes)

    # Retrieve list of DiagSpaceNames from all experiments
    expsDiagSpaceNames = []
    for expName in expLongNames:
        dateDir = cyDTimes_dir[0]
        if plotGroup == multiFCLen:
            dateDir = dateDir+'/'+fcTDeltas_dir[0]

        FILEPREFIX0 = expDirectory + expName +'/'+dateDir+'/' \
                      +statsFilePrefix+DAMethods[0]+"_"

        DiagSpaceNames = []
        for File in glob.glob(FILEPREFIX0+'*.nc'):
           DiagSpaceName = re.sub(".nc","",re.sub(FILEPREFIX0,"",File))
           DiagSpaceInfo_ = DiagSpaceDict.get(DiagSpaceName,vu.nullDiagSpaceInfo)
           if DiagSpaceInfo_['process']:
               DiagSpaceNames.append(DiagSpaceName)
        expsDiagSpaceNames.append(DiagSpaceNames)

    # Remove DiagSpaceNames that are not common to all experiments
    DiagSpaceNames = deepcopy(expsDiagSpaceNames[0])
    if len(expsDiagSpaceNames) > 1:
        for expDiagSpaceNames in expsDiagSpaceNames[1:]:
            for DiagSpaceName in expDiagSpaceNames:
                if not (DiagSpaceName in expDiagSpaceNames):
                    DiagSpaceNames.remove(DiagSpaceName)
    DiagSpaceNames.sort()

    print("")
    print("Processing data for these DiagSpaceNames:")
    print("---------------------------------------")
    print(DiagSpaceNames)

    # Read stats and make figures for all DiagSpaceNames
    for DiagSpaceName in DiagSpaceNames:
        print("\n==============================")
        print("DiagSpaceName = "+DiagSpaceName)
        print("==============================")

        print("\nReading NC intermediate files into common Pandas DataFrame...")

        dsDict = {}
        dsDict['expName'] = np.asarray([])
        dsDict['fcTDelta'] = np.asarray([])
        dsDict['cyDTime'] = np.asarray([])
        for attribName in su.fileStatAttributes:
            dsDict[attribName] = np.asarray([])
        for statName in su.allFileStats:
            dsDict[statName] = np.asarray([])

        for iexp, expName in enumerate(expNames):
            expPrefix = expDirectory + expLongNames[iexp] +'/'
            ncStatsFile = statsFilePrefix+DAMethods[iexp]+'_'+DiagSpaceName+'.nc'
            for ifc, fcTDelta in enumerate(fcTDeltas):
                fcstr=fcTDelta.__str__()
                for icy, cyDTime in enumerate(cyDTimes):
                    #Read all stats/attributes from NC file for DiagSpaceName, ExpName, fcTDelta, cyDTime
                    dateDir = cyDTimes_dir[icy]
                    if plotGroup == multiFCLen:
                        dateDir = dateDir+'/'+fcTDeltas_dir[ifc]
                    cyStatsFile = expPrefix+dateDir+'/'+ncStatsFile

                    statsDict = su.read_stats_nc(cyStatsFile)
                    nrows = len(statsDict[su.fileStatAttributes[0]])
                    dsDict['expName'] = \
                        np.append(dsDict['expName'], [expName] * nrows)
                    dsDict['fcTDelta'] = \
                        np.append(dsDict['fcTDelta'], [fcTDelta] * nrows)
                    dsDict['cyDTime'] = \
                        np.append(dsDict['cyDTime'], [cyDTime] * nrows)

                    for attribName in su.fileStatAttributes:
                        dsDict[attribName] = \
                            np.append(dsDict[attribName],statsDict[attribName])
                    for statName in su.allFileStats:
                        dsDict[statName] = \
                            np.append(dsDict[statName],statsDict[statName])

        #Convert dsDict to DataFrame
        dsDF = pd.DataFrame.from_dict(dsDict)

        indexNames = ['expName','fcTDelta','cyDTime','DiagSpaceGrp',
                      'varName','diagName','binVar','binVal','binMethod']

        dsDF.set_index(indexNames,inplace=True)
        dsDF.sort_index(inplace=True)

        ##  diagspace
        DiagSpaceGrp = dsDF.index.levels[indexNames.index('DiagSpaceGrp')]

        ##  variables
        # get varNames and sort alphabetically
        varNames = dsDF.index.levels[indexNames.index('varName')]
        nVars = len(varNames)
        indices = list(range(nVars))
        if DiagSpaceGrp[0] == vu.radiance_s:
            # sort by channel number (int) for radiances
            chlist = []
            for varName in varNames:
                chlist.append(int(''.join(varName.split("_")[-1:])))
            indices.sort(key=chlist.__getitem__)
        else:
            indices.sort(key=varNames.__getitem__)
        varNames = list(map(varNames.__getitem__, indices))

        ##  bin method (not used yet)
        # binMethods = dsDF.index.levels[indexNames.index('binMethod')].tolist()

        ##  bin variable (not used yet)
        # binVars = dsDF.index.levels[indexNames.index('binVar')].tolist()

        ##  bins (numerical values get sorted independently for each binVar)
        allBinVals = dsDF.index.levels[indexNames.index('binVal')].tolist()

        # extract units for all varNames from varUnits DF column
        varsLoc = (expNames[0],fcTDeltas[0],cyDTimes[0],DiagSpaceGrp,
                   varNames,diagNames[0],vu.obsVarQC,[vu.goodFlagName],slice(None))
        varUnitss = dsDF.loc[varsLoc,'varUnits'].tolist()

        # convert allBinVals to numeric type that can be used as axes values
        binNumVals = []
        binNumVals2DasStr = []
        for binVal in allBinVals:
            if pu.isfloat(binVal):
                binNumVals.append(float(binVal))
                binNumVals2DasStr.append(binVal)
            elif pu.isint(binVal):
                binNumVals.append(int(binVal))
                binNumVals2DasStr.append(binVal)
            else:
                binNumVals.append(np.NaN)

        #==================================
        # Create figures for all diagNames
        #==================================

        OSFigPath = "./"+DiagSpaceName+"_figs"
        if not os.path.exists(OSFigPath):
            os.mkdir(OSFigPath)

        for diagName in diagNames:
            print("\nCreate figures for diagName = "+diagName)
            print("---------------------------------------")

            fcDiagName = diagName
            if diagName == 'omb':
                fcDiagName = 'omf'
            elif diagName == 'bak':
                fcDiagName = 'fc'

            data_labels = []
            fcdata_labels = []
            for statName in statNames:
                if statName == 'Count':
                    data_labels.append(statName)
                    fcdata_labels.append(statName)
                else:
                    data_labels.append(statName+"("+diagName+")")
                    fcdata_labels.append(statName+"("+fcDiagName+")")


            signdef = []
            for statName in statNames:
                #This only applies to the unbounded quantities (omb, oma, ana/bak for velocity)
                if (statName == 'Mean' and
                    (diagName == 'omb' or diagName == 'oma')):
                    signdef.append(False)
                else:
                    signdef.append(True)

            # reduce the index space by two dimensions
            #            expName     fcTDelta    cyDTime                    varName              binVar      binVal      binMethod
            diagLoc = (slice(None),slice(None),slice(None),DiagSpaceGrp[0],slice(None),diagName,slice(None),slice(None),slice(None))
            diagDF = dsDF.xs(diagLoc)


            #=============
            # 1-D figures
            #=============
            if ("SNAPSHOTCY" in plotTypes or
                "AGGREGATEFC" in plotTypes or
                "AGGREGATEFC_DIFFCI" in plotTypes or
                "SNAPSHOTFC" in plotTypes):

                ## make one figure for all good/unbinned data
                selectBinVar    = vu.obsVarQC
                selectBinVals1D = [vu.goodFlagName]
                selectBinMethod = vu.defaultBinMethod

                ## make one figure for each non-float and non-int bin
                #selectBinVals1D = []
                #for binVal in allBinVals:
                #    if not pu.isfloat(binVal) and not pu.isint(binVal):
                #        selectBinVals1D.append(binVal)

                ## make one figure for each bin
                #selectBinVals1D = []
                #for binVal in allBinVals:
                #    if not (binVal in [vu.goodFlagName]):
                #        selectBinVals1D.append(binVal)

                ## make one figure for unbinned data and one for each bin
                #selectBinVals1D = allBinVals

                # top-level 1-D plot settings
                nsubplots = nVars
                nx = np.int(np.ceil(np.sqrt(nsubplots)))
                ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

                binLoc = (slice(None),slice(None),slice(None),slice(None),
                          selectBinVar,selectBinVals1D[0],selectBinMethod)
                binUnits = pu.uniqueMembers(diagDF.loc[binLoc,'binUnits'].tolist())[0]

            if ("SNAPSHOTCY" in plotTypes and
                len(selectBinVals1D) > 0 and
                nCY > 1):
                print("\nGenerating SNAPSHOTCY figures...")
                plotTypePath = OSFigPath+'/SNAPSHOTCY'
                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = allBinVals.index(binVal)
                        titlebin = " @ "+selectBinVar+"="+binVal
                        filebin = "_"+binVal
                        if binUnits != vu.miss_s:
                            titlebin = titlebin+" "+binUnits
                            filebin = filebin+binUnits
                    elif binVal in [vu.goodFlagName]:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #file loop 2
                    for ifc, fcTDelta in enumerate(fcTDeltas):
                        #file loop 3
                        for istat, statName in enumerate(statNames):
                            # establish a new figure
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,
                                          selectBinVar,selectBinVals1D,selectBinMethod)
                                dmin = diagDF.loc[varLoc,statName].dropna().min()
                                dmax = diagDF.loc[varLoc,statName].dropna().max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,
                                               selectBinVar,binVal,selectBinMethod)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel(
                                               ['expName','fcTDelta','varName','binVar','binVal','binMethod'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != vu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotTimeSeries(
                                    fig,
                                    cyDTimes, linesVals, expNames,
                                    title, data_labels[istat],
                                    sciticks[istat], signdef[istat],
                                    ny, nx, nsubplots, ivar,
                                    dmin = dmin, dmax = dmax,
                                    interiorLabels = interiorLabels)

                            # save each figure
                            if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                            filename = plotTypePath+'/TSeries_%sday_%s_%s_%s'%(
                                       fcTDeltas_dir[ifc],DiagSpaceName,
                                       diagName,statName)+filebin

                            pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                # end binName loop

            # end SNAPSHOTCY


            if ("AGGREGATEFC" in plotTypes and
                len(selectBinVals1D) > 0 and
                nFC > 1):
                print("\nGenerating AGGREGATEFC figures...")
                plotTypePath = OSFigPath+'/AGGREGATEFC'
                subplot_size = 1.9
                aspect = 0.6

                # apply aggregation over cyDTime via including all other free indices in groupby
                aggDF = diagDF.groupby(
                        ['expName','fcTDelta','varName','binVar','binVal','binMethod']).apply(su.aggStatsDF)

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = allBinVals.index(binVal)
                        titlebin = " @ "+selectBinVar+"="+binVal
                        filebin = "_"+binVal
                        if binUnits != vu.miss_s:
                            titlebin = titlebin+" "+binUnits
                            filebin = filebin+binUnits
                    elif binVal in [vu.goodFlagName]:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #file loop 2
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        #subplot loop
                        for ivar, varName in enumerate(varNames):
                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,
                                      selectBinVar,selectBinVals1D,selectBinMethod)
                            dmin = aggDF.loc[varLoc,statName].dropna().min()
                            dmax = aggDF.loc[varLoc,statName].dropna().max()

                            #collect aggregated statNames, varying across fcTDelta
                            linesVals = []
                            for expName in expNames:
                                #                                   fcTDelta
                                linesVals.append(aggDF.loc[(expName,slice(None),varName,
                                                            selectBinVar,binVal,selectBinMethod),
                                                            statName].to_numpy())

                            # define subplot title
                            title = varName
                            if varUnitss[ivar] != vu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"
                            title = title+titlebin

                            # perform subplot agnostic plotting (all expNames)
                            bpf.plotTimeSeries(
                                fig,
                                fcTDeltas, linesVals, expNames,
                                title, fcdata_labels[istat],
                                sciticks[istat], signdef[istat],
                                ny, nx, nsubplots, ivar,
                                dmin = dmin, dmax = dmax,
                                interiorLabels = interiorLabels)

                        # save each figure
                        if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                        filename = plotTypePath+'/TSeries_%s-%sday_%s_%s_%s'%(
                                   fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                   fcDiagName,statName)+filebin

                        pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                    # end statName loop

                # end binName loop

            # end AGGREGATEFC


            if ("AGGREGATEFC_DIFFCI" in plotTypes and
                len(selectBinVals1D) > 0 and
                nFC > 1 and nExp > 1):
                print("\nGenerating AGGREGATEFC_DIFFCI figures...")
                plotTypePath = OSFigPath+'/AGGREGATEFC_DIFFCI'
                subplot_size = 1.9
                aspect = 0.6

                # nx = 1
                # ny = nVars
                # nsubplots = nx * ny

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = allBinVals.index(binVal)
                        titlebin = " @ "+selectBinVar+"="+binVal
                        filebin = "_"+binVal
                        if binUnits != vu.miss_s:
                            titlebin = titlebin+" "+binUnits
                            filebin = filebin+binUnits
                    elif binVal in [vu.goodFlagName]:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #figure loop 2
                    for statName in bootStrapStats:

                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            # define subplot title
                            title = varName
                            if varUnitss[ivar] != vu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"
                            title = title+titlebin

                            linesVals = {}
                            for trait in su.ciTraits: linesVals[trait] = []

                            for expName in noncntrlExpNames:

                                lineVals = {}
                                for trait in su.ciTraits: lineVals[trait] = []

                                for fcTDelta in fcTDeltas:
                                    #                              cyDTime
                                    cntrlLoc = (cntrlName,fcTDelta,slice(None),varName,
                                                selectBinVar,binVal,selectBinMethod)

                                    #                          cyDTime
                                    expLoc = (expName,fcTDelta,slice(None),varName,
                                              selectBinVar,binVal,selectBinMethod)

                                    ciVals = su.bootStrapClusterFunc(
                                                 X = diagDF.loc[expLoc,:],
                                                 Y = diagDF.loc[cntrlLoc,:],
                                                 n_samples = 10000,
                                                 statNames = [statName])

                                    for trait in su.ciTraits:
                                        lineVals[trait].append(ciVals[statName][trait][0])

                                for trait in su.ciTraits:
                                    linesVals[trait].append(
                                        lineVals[trait] )

                            # use specific y-axis limits for each varName
                            dmin = np.nanmin(linesVals[su.cimin])
                            dmax = np.nanmax(linesVals[su.cimax])

                            # perform subplot agnostic plotting (all expNames)
                            bpf.plotTimeSeries(
                                fig,
                                fcTDeltas, linesVals[su.cimean],
                                noncntrlExpNames,
                                title,
                                statName+"("+fcDiagName+"): [EXP - CTRL]",
                                False, False,
                                ny, nx, nsubplots, ivar,
                                linesValsMinCI = linesVals[su.cimin],
                                linesValsMaxCI = linesVals[su.cimax],
                                dmin = dmin, dmax = dmax,
                                lineAttribOffset = 1,
                                interiorLabels = interiorLabels)

                        # end varName loop

                        # save each figure
                        if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                        filename = plotTypePath+'/TSeries_%s-%sday_%s_%s_%s'%(
                                   fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                   fcDiagName,statName)+filebin

                        pu.finalize_fig(fig, filename, figureFileType, interiorLabels)
                    # end statName loop

                # end binName loop

            # end AGGREGATEFC

            if ("SNAPSHOTFC" in plotTypes and
                len(selectBinVals1D) > 0 and
                nFC > 1 and nCY > 1):
                print("\nGenerating SNAPSHOTFC figures...")
                plotTypePath = OSFigPath+'/SNAPSHOTFC'
                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = allBinVals.index(binVal)
                        titlebin = " @ "+selectBinVar+"="+binVal
                        filebin = "_"+binVal
                        if binUnits != vu.miss_s:
                            titlebin = titlebin+" "+binUnits
                            filebin = filebin+binUnits
                    elif binVal in [vu.goodFlagName]:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #file loop 2
                    for iexp, expName in enumerate(expNames):
                        #file loop 3
                        for istat, statName in enumerate(statNames):
                            # establish a new figure
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,
                                          selectBinVar,selectBinVals1D,selectBinMethod)
                                dmin = diagDF.loc[varLoc,statName].dropna().min()
                                dmax = diagDF.loc[varLoc,statName].dropna().max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                xsVals = []
                                linesVals = []
                                fcTDeltas_labels = []
                                for ifc, fcTDelta in enumerate(fcTDeltas):
                                    xVals = []
                                    for cyDTime in cyDTimes:
                                        xVals.append(cyDTime+fcTDelta)
                                    xsVals.append(xVals)

                                    #Setting to avoid over-crowding
                                    if ifc > (MAX_FC_LINES-1): continue

                                    fcTDelta_sec = pu.TDeltas2Seconds([fcTDelta])
                                    fcTDeltas_labels.append(pu.timeTicks(fcTDelta_sec[0],0))
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,
                                               selectBinVar,binVal,selectBinMethod)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel(
                                               ['expName','fcTDelta','varName','binVar','binVal','binMethod'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != vu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotTimeSeries(
                                    fig,
                                    xsVals, linesVals, fcTDeltas_labels,
                                    title, data_labels[istat],
                                    sciticks[istat], signdef[istat],
                                    ny, nx, nsubplots, ivar,
                                    dmin = dmin, dmax = dmax,
                                    interiorLabels = interiorLabels)


                            expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                            if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                            filename = plotTypePath+'/TSeries_%s_%s_%s_%s'%(
                                       expName_file,DiagSpaceName,
                                       fcDiagName,statName)+filebin

                            pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                # end binName loop

            # end SNAPSHOTFC

            selectBinVar = 'lat'
            selectBinVals1D = vu.namedLatBands['labels']
            selectBinMethod = 'NAMED'
            if ("SNAPSHOTCY_LAT1" in plotTypes and
                len(selectBinVals1D) > 0 and
                nCY > 1 and
                DiagSpaceGrp[0] != vu.radiance_s):

                print("\nGenerating SNAPSHOTCY_LAT1 figures...")
                plotTypePath = OSFigPath+'/SNAPSHOTCY_LAT1'

                subplot_size = 1.9
                aspect = 0.75

                nsubplots = nExp * nVars
                nx = nExp
                ny = nVars

                #file loop 1
                for ifc, fcTDelta in enumerate(fcTDeltas):
                    #file loop 2
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        iplot = 0

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            #subplot loop 2
                            for expName in expNames:

                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,
                                          selectBinVar,selectBinVals1D,selectBinMethod)
                                dmin = diagDF.loc[varLoc,statName].dropna().min()
                                dmax = diagDF.loc[varLoc,statName].dropna().max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for binVal in selectBinVals1D:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,
                                               selectBinVar,binVal,selectBinMethod)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel(
                                               ['expName','fcTDelta','varName','binVar','binVal','binMethod'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = expName+"\n"+varName
                                if varUnitss[ivar] != vu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotTimeSeries(
                                    fig,
                                    cyDTimes, linesVals, selectBinVals1D,
                                    title, data_labels[istat],
                                    sciticks[istat], signdef[istat],
                                    ny, nx, nsubplots, iplot,
                                    dmin = dmin, dmax = dmax,
                                    interiorLabels = interiorLabels)

                                iplot = iplot + 1
                        if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                        filename = plotTypePath+'/TSeries_%sday_%s_%s_%s'%(
                                   fcTDeltas_dir[ifc],DiagSpaceName,
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                # end binName loop

            # end SNAPSHOTCY_LAT1

            selectBinVar = 'lat'
            selectBinVals1D = vu.namedLatBands['labels']
            selectBinMethod = 'NAMED'
            if ("SNAPSHOTCY_LAT2" in plotTypes and
                len(selectBinVals1D) > 0 and
                nCY > 1 and
                DiagSpaceGrp[0] != vu.radiance_s):
                print("\nGenerating SNAPSHOTCY_LAT2 figures...")
                plotTypePath = OSFigPath+'/SNAPSHOTCY_LAT2'
                nsubplots = len(selectBinVals1D)*nVars
                nx = len(selectBinVals1D)
                ny = nVars

                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for ifc, fcTDelta in enumerate(fcTDeltas):
                    #file loop 2
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        iplot = 0
                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):


                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),slice(None),varName,
                                      selectBinVar,selectBinVals1D,selectBinMethod)
                            dmin = diagDF.loc[varLoc,statName].dropna().min()
                            dmax = diagDF.loc[varLoc,statName].dropna().max()

                            #subplot loop 2
                            for binVal in selectBinVals1D:

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,
                                               selectBinVar,binVal,selectBinMethod)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel(
                                               ['expName','fcTDelta','varName','binVar','binVal','binMethod'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = binVal+"\n"+varName
                                if varUnitss[ivar] != vu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotTimeSeries(
                                    fig,
                                    cyDTimes, linesVals, expNames,
                                    title, data_labels[istat],
                                    sciticks[istat], signdef[istat],
                                    ny, nx, nsubplots, iplot,
                                    dmin = dmin, dmax = dmax,
                                    interiorLabels = interiorLabels)

                                iplot = iplot + 1

                        #expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                        if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                        filename = plotTypePath+'/TSeries_%sday_%s_%s_%s'%(
                                   fcTDeltas_dir[ifc],DiagSpaceName,
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                # end binName loop

            # end SNAPSHOTCY_LAT2

            selectBinVar = vu.obsVarQC
            selectBinVals1D = [vu.goodFlagName] + vu.badFlagNames
            selectBinMethods = [vu.defaultBinMethod,'bad']
            if ("SNAPSHOTCY_QC" in plotTypes and
                len(selectBinVals1D) > 0 and
                nCY > 1):
                print("\nGenerating SNAPSHOTCY_QC figures...")
                plotTypePath = OSFigPath+'/SNAPSHOTCY_QC'
                subplot_size = 1.9
                aspect = 0.75

                nsubplots = nExp * nVars
                nx = nExp
                ny = nVars

                #file loop 1
                for ifc, fcTDelta in enumerate(fcTDeltas):
                    statName = 'Count'
                    istat = statNames.index(statName)

                    # establish a new figure
                    fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                    iplot = 0

                    #subplot loop 1
                    for ivar, varName in enumerate(varNames):
                        #subplot loop 2
                        for expName in expNames:

                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),slice(None),varName,
                                      selectBinVar,selectBinVals1D,selectBinMethods)
                            dmin = diagDF.loc[varLoc,statName].dropna().min()
                            dmax = diagDF.loc[varLoc,statName].dropna().max()

                            # collect statName for all lines on this subplot, letting cyDTime vary
                            linesVals = []
                            for binVal in selectBinVals1D:
                                #                           cyDTime
                                dataLoc = (expName,fcTDelta,slice(None),varName,
                                           selectBinVar,binVal,selectBinMethods)
                                dataDF = diagDF.loc[dataLoc,statName].droplevel(
                                           ['expName','fcTDelta','varName','binVar','binVal'])
                                dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                lineVals = np.empty(nCY)
                                lineVals[:] = np.NaN
                                for cyDTime in dataCYDTimes:
                                    icy = cyDTimes.index(cyDTime)
                                    lineVals[icy] = dataDF.loc[cyDTime]
                                linesVals.append(lineVals)

                            # define subplot title
                            title = expName+"\n"+varName
                            if varUnitss[ivar] != vu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"

                            # perform subplot agnostic plotting (all expNames)
                            bpf.plotTimeSeries(
                                fig,
                                cyDTimes, linesVals, selectBinVals1D,
                                title, data_labels[istat],
                                sciticks[istat], signdef[istat],
                                ny, nx, nsubplots, iplot,
                                dmin = dmin, dmax = dmax,
                                legend_inside = False,
                                interiorLabels = interiorLabels)

                            iplot = iplot + 1

                    if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                    filename = plotTypePath+'/TSeries_%sday_%s_%s_%s'%(
                               fcTDeltas_dir[ifc],DiagSpaceName,
                               diagName,statName)

                    pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                # end binName loop

            # end SNAPSHOTCY_QC


            #=============
            # 2-D figures
            #=============
            diagMI = diagDF.index.names
            binVars = pu.uniqueMembers(
                          diagDF.index.get_level_values(
                              diagMI.index('binVar') ).tolist() )

            PLOT_ACROSS_BINS = (
                "SNAPSHOTCY2D" in plotTypes or
                "AGGREGATEFC2D" in plotTypes or
                "AGGREGATEFCPROFILE" in plotTypes or
                "AGGREGATEFCPROFILE_DIFFCI" in plotTypes)

            if PLOT_ACROSS_BINS:
                for binVar2D, binVarDict in binVars2D.items():
                    if (binVar2D not in binVars): continue

                    #Get all float/int binVals associated with binVar2D
                    binVarLoc = (slice(None),slice(None),slice(None),slice(None),
                                 binVar2D,binNumVals2DasStr,slice(None))
                    binVarDF = diagDF.loc[binVarLoc,:].droplevel(['binVar'])

                    # apply aggregation over cyDTime via including all other free indices in groupby
                    aggVarDF = binVarDF.groupby(
                            ['expName','fcTDelta','varName','binVal','binMethod']).apply(su.aggStatsDF)

                    # determine binMethods for this slice
                    binVarMI = binVarDF.index.names
                    selectBinMethods = pu.uniqueMembers(
                                          binVarDF.index.get_level_values(
                                              binVarMI.index('binMethod') ).tolist() )

                    for binMethod in selectBinMethods:
                        binMethodLoc = (slice(None),slice(None),slice(None),slice(None),
                                        slice(None),binMethod)
                        binMethodDF = binVarDF.loc[binMethodLoc,:].droplevel(['binMethod'])
                        aggMethodLoc = (slice(None),slice(None),slice(None),
                                        slice(None),binMethod)
                        aggMethodDF = aggVarDF.loc[aggMethodLoc,:].droplevel(['binMethod'])

                        # determine binVals for this slice
                        binMethodMI = binMethodDF.index.names
                        selectBinVals2D = pu.uniqueMembers(
                                              binMethodDF.index.get_level_values(
                                                  binMethodMI.index('binVal') ).tolist() )

                        # determine binUnits for this slice
                        binUnits = pu.uniqueMembers(binMethodDF.loc[:,'binUnits'].tolist())[0]

                        binMethodFile = ''
                        if binMethod != vu.defaultBinMethod: binMethodFile = '_'+binMethod

                        # bin info
                        selectBinNumVals2D = []
                        for binVal in selectBinVals2D:
                            ibin = allBinVals.index(binVal)
                            selectBinNumVals2D.append(binNumVals[ibin])

                            # assume all bins represent same variable/units
                            indepLabel = binVar2D
                            if binUnits != vu.miss_s:
                                indepLabel = indepLabel+" ("+binUnits+")"

                            # invert independent variable axis for pressure bins
                            pressure_dict = vu.varDictObs.get(vu.obsVarPrs,['',''])
                            invert_ind_axis = (pressure_dict[1] == binVar2D)

                        # sort bins by numeric value
                        indices = list(range(len(selectBinNumVals2D)))
                        indices.sort(key=selectBinNumVals2D.__getitem__)
                        selectBinNumVals2D = list(map(selectBinNumVals2D.__getitem__, indices))
                        selectBinVals2D = list(map(selectBinVals2D.__getitem__, indices))


                        if ("SNAPSHOTCY2D" in plotTypes and
                            len(selectBinVals2D) > 0 and
                            nCY > 1):

                            print("\nGenerating SNAPSHOTCY2D figures across "+binVar2D+" and binMethod=>"+binMethod)
                            plotTypePath = OSFigPath+'/SNAPSHOTCY2D'
                            # top-level 2-D plot settings
                            subplot_size = 2.4
                            aspect = 0.65
                            nx = nExp
                            ny = nVars
                            nsubplots = nx * ny

                            #file loop 1
                            for ifc, fcTDelta in enumerate(fcTDeltas):
                                #print("fcTDelta = "+fcTDeltas_dir[ifc])

                                #file loop 2
                                for istat, statName in enumerate(statNames):
                                    # establish a new figure
                                    fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                    iplot = 0

                                    #subplot loop 1
                                    for ivar, varName in enumerate(varNames):
                                        # use specific c-axis limits for each varName
                                        varLoc = (slice(None),slice(None),slice(None),varName,
                                                  selectBinVals2D,slice(None))
                                        dmin = binVarDF.loc[varLoc,statName].dropna().min()
                                        dmax = binVarDF.loc[varLoc,statName].dropna().max()

                                        #subplot loop 2
                                        # letting cyDTime and binVal vary
                                        for expName in expNames:
                                            contourVals = np.empty((len(selectBinVals2D),nCY))
                                            contourVals[:,:] = np.NaN
                                            #                           cyDTime
                                            dataLoc = (expName,fcTDelta,slice(None),varName,selectBinVals2D)
                                            dataDF = binMethodDF.loc[dataLoc,statName].droplevel(
                                                       ['expName','fcTDelta','varName'])
                                            dataCYDTimes = dataDF.loc[(slice(None),selectBinVals2D[0])
                                                                     ].index.get_level_values('cyDTime')
                                            for ibin, binVal in enumerate(selectBinVals2D):
                                                tmp = dataDF.loc[(slice(None),binVal)].to_numpy()
                                                for jcy, cyDTime in enumerate(dataCYDTimes):
                                                    icy = cyDTimes.index(cyDTime)
                                                    contourVals[ibin,icy] = tmp[jcy]

                                            # define subplot title
                                            title = expName+"\n"+varName
                                            if varUnitss[ivar] != vu.miss_s:
                                                title = title+" ("+varUnitss[ivar]+")"

                                            # perform subplot agnostic plotting (all expNames)
                                            bpf.plotTimeSeries2D(
                                                fig,
                                                cyDTimes, selectBinNumVals2D, contourVals,
                                                title, data_labels[istat],
                                                sciticks[istat], signdef[istat],
                                                indepLabel, invert_ind_axis,
                                                ny, nx, nsubplots, iplot,
                                                dmin = dmin, dmax = dmax,
                                                interiorLabels = interiorLabels)

                                            iplot = iplot + 1

                                    if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                                    filename = plotTypePath+'/%s%s_TSeries_%sday_%s_%s_%s'%(
                                               binVar2D,binMethodFile,fcTDeltas_dir[ifc],DiagSpaceName,
                                               diagName,statName)

                                    pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                                # end statName loop

                            # end fcTDelta loop

                        # end SNAPSHOTCY2D


                        if ("AGGREGATEFC2D" in plotTypes and
                            len(selectBinVals2D) > 0 and
                            nFC > 1):

                            print("\nGenerating AGGREGATEFC2D figures across "+binVar2D+" and binMethod=>"+binMethod)
                            plotTypePath = OSFigPath+'/AGGREGATEFC2D'
                            # top-level 2-D plot settings
                            subplot_size = 2.4
                            aspect = 0.55
                            nx = nExp
                            ny = nVars
                            nsubplots = nx * ny

                            #file loop 1
                            for istat, statName in enumerate(statNames):
                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                iplot = 0
                                #subplot loop
                                for ivar, varName in enumerate(varNames):
                                    # use specific c-axis limits for each varName
                                    varLoc = (slice(None),slice(None),varName,selectBinVals2D,slice(None))
                                    dmin = aggVarDF.loc[varLoc,statName].dropna().min()
                                    dmax = aggVarDF.loc[varLoc,statName].dropna().max()

                                    #collect aggregated statName, varying across fcTDelta+binVal
                                    for expName in expNames:
                                        contourVals = np.empty((len(selectBinVals2D),nFC))
                                        contourVals[:,:] = np.NaN
                                        #                  fcTDelta
                                        dataLoc = (expName,slice(None),varName,selectBinVals2D)
                                        dataDF = aggMethodDF.loc[dataLoc,statName].droplevel(
                                                   ['expName','varName'])
                                        dataFCTDeltas = dataDF.loc[(slice(None),selectBinVals2D[0])
                                                                  ].index.get_level_values('fcTDelta')
                                        for ibin, binVal in enumerate(selectBinVals2D):
                                            tmp = dataDF.loc[(slice(None),binVal)].to_numpy()
                                            for jfc, fcTDelta in enumerate(dataFCTDeltas):
                                                ifc = fcTDeltas.index(fcTDelta)
                                                contourVals[ibin,ifc] = tmp[jfc]

                                        # define subplot title
                                        title = expName+"\n"+varName
                                        if varUnitss[ivar] != vu.miss_s:
                                            title = title+" ("+varUnitss[ivar]+")"

                                        # perform subplot agnostic plotting (all expNames)
                                        bpf.plotTimeSeries2D(
                                            fig,
                                            fcTDeltas, selectBinNumVals2D, contourVals,
                                            title, fcdata_labels[istat],
                                            sciticks[istat], signdef[istat],
                                            indepLabel, invert_ind_axis,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)

                                        iplot = iplot + 1

                                # save each figure
                                if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                                filename = plotTypePath+'/%s%s_TSeries_%s-%sday_%s_%s_%s'%(
                                           binVar2D,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName)

                                pu.finalize_fig(fig, filename, figureFileType, interiorLabels)

                            # end statName loop

                        # end AGGREGATEFC2D


                        if ("AGGREGATEFCPROFILE" in plotTypes and
                            len(selectBinVals2D) > 0):
                            print("\nGenerating AGGREGATEFCPROFILE figures across "+binVar2D+" and binMethod=>"+binMethod)
                            plotTypePath = OSFigPath+'/AGGREGATEFCPROFILE'
                            # top-level 2-D plot settings
                            subplot_size = 1.2
                            aspect = 1.3

                            if nFC > 1:
                                nx = min([nFC,MAX_FC_SUBFIGS])
                                ny = nVars
                                nsubplots = nx * ny
                            else:
                                nsubplots = nVars
                                nx = np.int(np.ceil(np.sqrt(nsubplots)))
                                ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

                            #file loop 1
                            for istat, statName in enumerate(statNames):
                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                iplot = 0

                                #subplot loop 1
                                for ivar, varName in enumerate(varNames):
                                    # use specific x-axis limits for each varName
                                    varLoc = (slice(None),slice(None),varName,selectBinVals2D,slice(None))
                                    dmin = aggVarDF.loc[varLoc,statName].dropna().min()
                                    dmax = aggVarDF.loc[varLoc,statName].dropna().max()

                                    #subplot loop 2
                                    for ifc, fcTDelta in enumerate(fcTDeltas):

                                        #Setting to avoid over-crowding
                                        if ifc > (MAX_FC_SUBFIGS-1): continue

                                        #collect aggregated statNames, varying across fcTDelta
                                        linesVals = []
                                        for expName in expNames:
                                            expVals = []
                                            for binVal in selectBinVals2D:
                                                binLoc = (expName,fcTDelta,varName,binVal)
                                                expVals.append(aggMethodDF.loc[binLoc,statName])
                                            linesVals.append(expVals)

                                        # define subplot title
                                        title = varName
                                        if varUnitss[ivar] != vu.miss_s:
                                            title = title+" ("+varUnitss[ivar]+")"
                                        title = title+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                        # perform subplot agnostic plotting (all expNames)
                                        binVarDict['profileFunc'](
                                            fig,
                                            linesVals, selectBinNumVals2D,
                                            expNames,
                                            title, fcdata_labels[istat],
                                            sciticks[istat], signdef[istat],
                                            indepLabel, invert_ind_axis,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)


                                        iplot = iplot + 1

                                # save each figure
                                if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                                filename = plotTypePath+'/%s%s_TSeries_%s-%sday_%s_%s_%s'%(
                                           binVar2D,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName)

                                pu.finalize_fig(fig, filename, figureFileType, interiorLabels, True)

                            # end statName loop

                        # end AGGREGATEFCPROFILE


                        if ("AGGREGATEFCPROFILE_DIFFCI" in plotTypes and
                            len(selectBinVals2D) > 0 and
                            nExp > 1):
                            print("\nGenerating AGGREGATEFCPROFILE_DIFFCI figures across "+binVar2D+" and binMethod=>"+binMethod)
                            plotTypePath = OSFigPath+'/AGGREGATEFCPROFILE_DIFFCI'
                            # top-level 2-D plot settings
                            subplot_size = 1.2
                            aspect = 1.3

                            if nFC > 1:
                                nx = min([nFC,MAX_FC_SUBFIGS])
                                ny = nVars
                                nsubplots = nx * ny
                            else:
                                nsubplots = nVars
                                nx = np.int(np.ceil(np.sqrt(nsubplots)))
                                ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

                            #file loop 1
                            for statName in bootStrapStats:
                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                iplot = 0

                                #subplot loop 1
                                for ivar, varName in enumerate(varNames):
                                    #subplot loop 2
                                    for ifc, fcTDelta in enumerate(fcTDeltas):

                                        #Setting to avoid over-crowding
                                        if ifc > (MAX_FC_SUBFIGS-1): continue

                                        linesVals = {}
                                        for trait in su.ciTraits: linesVals[trait] = []
                                        for expName in noncntrlExpNames:
                                            lineVals = {}
                                            for trait in su.ciTraits: lineVals[trait] = []

                                            for binVal in selectBinVals2D:
                                                #                              cyDTime
                                                cntrlLoc = (cntrlName,fcTDelta,slice(None),varName,binVal)

                                                #                          cyDTime
                                                expLoc = (expName,fcTDelta,slice(None),varName,binVal)

                                                ciVals = su.bootStrapClusterFunc(
                                                             X = binMethodDF.loc[expLoc,:],
                                                             Y = binMethodDF.loc[cntrlLoc,:],
                                                             n_samples = 10000,
                                                             statNames = [statName])

                                                for trait in su.ciTraits:
                                                    lineVals[trait].append(ciVals[statName][trait][0])

                                            for trait in su.ciTraits:
                                                linesVals[trait].append(lineVals[trait])

                                        # define subplot title
                                        title = varName
                                        if varUnitss[ivar] != vu.miss_s:
                                            title = title+" ("+varUnitss[ivar]+")"
                                        title = title+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                        # use specific y-axis limits for each varName
                                        dmin = np.nanmin(linesVals[su.cimin])
                                        dmax = np.nanmax(linesVals[su.cimax])

                                        # perform subplot agnostic plotting (all expNames)
                                        binVarDict['profileFunc'](
                                            fig,
                                            linesVals[su.cimean], selectBinNumVals2D,
                                            noncntrlExpNames,
                                            title,
                                            statName+"("+fcDiagName+"): [EXP - CTRL]",
                                            False, False,
                                            indepLabel, invert_ind_axis,
                                            ny, nx, nsubplots, iplot,
                                            linesValsMinCI = linesVals[su.cimin],
                                            linesValsMaxCI = linesVals[su.cimax],
                                            dmin = dmin, dmax = dmax,
                                            lineAttribOffset = 1,
                                            interiorLabels = interiorLabels)

                                        iplot = iplot + 1

                                # save each figure
                                if not os.path.exists(plotTypePath): os.mkdir(plotTypePath)
                                filename = plotTypePath+'/%s%s_TSeries_%s-%sday_%s_%s_%s'%(
                                           binVar2D,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName)

                                pu.finalize_fig(fig, filename, figureFileType, interiorLabels, True)

                            # end statName loop

                        # end AGGREGATEFCPROFILE_DIFFCI

                    # end binMethods loop

                # end binVars2D loop

            # end PLOT_ACROSS_BINS

        # end diagName loop

    # end DiagSpaceName loop

def main():
    nproc, myproc = paru.par_args(sys.argv[:])
    plot_stats_timeseries(nproc, myproc)

if __name__ == '__main__': main()
