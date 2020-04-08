import basic_plot_functions as bpf
import binning_utils as bu
import binning_configs as bcs
from collections.abc import Iterable
import config as conf
from copy import deepcopy
import collections
import datetime as dt
import glob
import numpy as np
import par_utils as paru
from pathlib import Path
import plot_utils as pu
import re
import os
import stat_utils as su
import StatsDB as sdb
from StatsDB import DFWrapper, StatsDataBase
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
# SNAPSHOTCY2D/AGGREGATEFC2D creates raster maps similar to above with metadata variable on y-axis
#             - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
#             - subplot: column by experiment, row by diag space variable
#             SNAPSHOTCY2D
#             -    file: combination of FC lead time and statistic
#             AGGREGATEFC2D
#             -    file: per statistic
# AGGREGATEFC_Profile similar to AGGREGATEFC2D, except
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
# SNAPSHOTCY_LatLines similar to SNAPSHOTCY, except
#             -    line: named latitude bins
#             - subplot: column by experiment, row by diag space variable
#             -    file: combination of statistic and forecast length
#             - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_QCLines similar to SNAPSHOTCY, except
#             - only plots 'Count' statistic
#             -    line: named QC bins
#             - subplot: column by experiment, row by diag space variable
#             -    file: forecast length
# AGGREGATEFC_DiffCI/AGGREGATEFC_Profile_DiffCI similar to their per experiment counterparts
#             - shows difference between experiment(s) and control
#             - control is selected using cntrlExpInd
cntrlExpInd = 0
#             - statistics are selected with bootStrapStats
bootStrapStats = []
for x in su.sampleableAggStats:
    if x != 'Count': bootStrapStats.append(x)

#TODO: AGGREGATEFC_LatLines figures w/ FC lead time x-axis (low priority)

# NOTE: all non-Profile *FC* figures require non-zero forecast length
# NOTE: all *CY* figures require > 1 analysis cycle

#Select the stats for plotting
# options: see su.allFileStats
statNames = ['Count','Mean','RMS','STD']

plotGroup = sdb.singleFCLen
#plotGroup = sdb.multiFCLen

#TODO: move the following configuration step to a wrapper function
#      possibly take arguments from the command line
#TODO: use multiprocessing module instead of requiring GNU parallel,
#      similar to writediagstats_obsspace.py
#TODO: move independent plotLoops (1-D, 2-D, QCLines, LatLines, PDF, statsComposite) to
#      individual classes that take a StatsDataBase object as an argument
#TODO: multiple diagNames on same subplot
#      e.g., diagNameGroups = [['omb','oma'],['obs','bak','ana']]

#Note: refer to the StatsDataBase class to configure the
#      directory structure below

plotTypes = []

if plotGroup == sdb.singleFCLen:
    user = 'guerrett'

    expLongNames = [
#                 'cycling_30km_omb_conv_clear_abi/DAdiag',
                 'cycling_30km_omb_conv_cloud_abi/DAdiag',
                 ]

    expNames = [
#                 'clrsky-CRTM',
                 'allsky-CRTM',
                 ]

    DAMethods = [
#                 'omb',
                 'omb',
                 ]

    #First and Last CYCLE dates
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,4,22,0,0,0)
    cyTimeInc  = dt.timedelta(hours=6)

    #First and Last FORECAST durations
    fcTDeltaFirst = dt.timedelta(days=0)
    fcTDeltaLast  = dt.timedelta(days=0)
    fcTimeInc  = dt.timedelta(hours=24)

    ## plotTypes for considering single forecast length
    ## ------------------------------------------------------------------------
    plotTypes.append('SNAPSHOTCY')
    plotTypes.append('SNAPSHOTCY_QCLines')
    # plotTypes.append('SNAPSHOTCY2D')


if plotGroup == sdb.multiFCLen:
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

    # plotTypes for considering multiple forecast lengths
    # --------------------------------------------------------------------------
    plotTypes.append('AGGREGATEFC')
    plotTypes.append('SNAPSHOTFC')
    # plotTypes.append('AGGREGATEFC2D')

nExp  = len(expNames)

#AGGREGATEFC_Profile* figures work for any forecast duration
plotTypes.append('AGGREGATEFC_Profile')
plotTypes.append('AGGREGATEFC_PDF')
#plotTypes.append('AGGREGATEFC_StatsComposite')
plotTypes.append('CalcGrossStats')
if nExp > 1:
    plotTypes.append('AGGREGATEFC_Profile_DiffCI')

    cntrlName = expNames[min([cntrlExpInd,nExp-1])]
    noncntrlExpNames = [x for x in expNames if x != cntrlName]
    print('\nControl Experiment: '+cntrlName)
    print('\nNon-control Experiment(s): ')
    print(noncntrlExpNames)

expDirectory = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac/')

statsFilePrefix = 'diagnostic_stats/'+su.statsFilePrefix


## plot settings
figureFileType = 'pdf' #['pdf','png']

interiorLabels = True

sciTickss = []
for statName in statNames:
    if statName == 'Count': sciTickss.append(True)
    else: sciTickss.append(False)

#
# primary script to be called by main()
#

def plot_stats_timeseries(nproc=1, myproc=0):
#  nproc - total number of processors
#  myproc - processor rank, starting at 0
#  Note: these arguments are used when the calling program has multiple processors available.

    # Assign processors round-robin to each DiagSpace name
    DiagSpaceConfig = {}
    jj=-1
    for key, baseval in conf.DiagSpaceConfig.items():
        if not baseval['process']: continue
        jj = jj + 1
        if jj%nproc != myproc: continue
        DiagSpaceConfig[key] = deepcopy(baseval)

    if fcTDeltaFirst == fcTDeltaLast and firstCycleDTime == lastCycleDTime:
        print("\n\n===============================================================================")
        print("===============================================================================")
        print("\nWARNING: Only Profile plotTypes can be generated without a time difference")
        print("         (either forecast or multiple cycles)\n")
        print("===============================================================================")
        print("===============================================================================\n\n")

    ## Configure the statistics database
    dbConf = {}

    # cycle DateTimes
    dbConf['firstCycleDTime'] = firstCycleDTime
    dbConf['lastCycleDTime'] = lastCycleDTime
    dbConf['cyTimeInc'] = cyTimeInc

    # forecast TimeDeltas
    dbConf['fcTDeltaFirst'] = fcTDeltaFirst
    dbConf['fcTDeltaLast'] = fcTDeltaLast
    dbConf['fcTimeInc'] = fcTimeInc

    # experiment info
    dbConf['expDirectory'] = expDirectory
    dbConf['expLongNames'] = expLongNames
    dbConf['expNames'] = expNames
    dbConf['DAMethods'] = DAMethods

    # group of selected plots (single/multi)
    dbConf['plotGroup'] = plotGroup

    # Read stats and make figures for all DiagSpaceNames
    for DiagSpaceName in sorted(DiagSpaceConfig):
        dbConf['DiagSpaceName'] = DiagSpaceName
        diagNames = DiagSpaceConfig[DiagSpaceName]['diagNames']

        ## construct and initialize the statistical database
        statsDB = StatsDataBase(dbConf)

        statsDB.initialize()

        if not statsDB.available: continue

        ## Extract useful variables from the database
        fcTDeltas = statsDB.fcTDeltas
        fcTDeltas_dir = statsDB.fcTDeltas_dir
        fcMap = list(zip(fcTDeltas,fcTDeltas_dir))
        nFC = len(fcTDeltas)

        cyDTimes = statsDB.cyDTimes
        nCY = len(cyDTimes)

        varNames = statsDB.varNames
        nVars = len(varNames)

        varLabels = []
        for (varName, varUnits) in zip(varNames,statsDB.varUnitss):
            label = varName
            if varUnits != vu.miss_s:
                label = label+" ("+varUnits+")"
            varLabels.append(label)
        varMap = list(zip(varNames,varLabels))

        allBinVals = statsDB.allBinVals
        allDiagNames = statsDB.allDiagNames
        binNumVals = statsDB.binNumVals
        binNumVals2DasStr = statsDB.binNumVals2DasStr

        #=======================================
        # Create figures for all diagNames
        #=======================================

        OSFigPath = Path("./"+DiagSpaceName+"_figs")
        OSFigPath.mkdir(parents=True, exist_ok=True)

        for diagName in diagNames:
            if diagName not in allDiagNames: continue

            print("\n\n---------------------------------------")
            print("Create figures for diagName = "+diagName)
            print("---------------------------------------")

            fcDiagName = diagName
            if nFC > 1 or fcTDeltaFirst > dt.timedelta(0):
                if diagName == 'omb':
                    fcDiagName = 'omf'
                elif diagName == 'bak':
                    fcDiagName = 'fc'

            statDiagLabels = []
            fcstatDiagLabels = []
            for statName in statNames:
                if statName == 'Count':
                    statDiagLabels.append(statName)
                    fcstatDiagLabels.append(statName)
                else:
                    statDiagLabels.append(statName+"("+diagName+")")
                    fcstatDiagLabels.append(statName+"("+fcDiagName+")")


            signDefinites = []
            for statName in statNames:
                #This only applies to the unbounded quantities (omb, oma, ana/bak for velocity)
                if (statName == 'Mean' and
                    (diagName == 'omb' or diagName == 'oma')):
                    signDefinites.append(False)
                else:
                    signDefinites.append(True)

            StatsMap = list(zip(statNames, statDiagLabels, sciTickss, signDefinites))
            FCStatsMap = list(zip(statNames, fcstatDiagLabels, sciTickss, signDefinites))

            diagLoc = {'diagName':diagName}
            diagDFW = statsDB.loc(diagLoc)
            diagBinVars = diagDFW.levels('binVar')
            diagBinMethods = diagDFW.levels('binMethod')

            # apply aggregation over cyDTime
            aggCYDFW = DFWrapper.fromAggStats(diagDFW,['cyDTime'])

            #=============
            # 1-D figures
            #=============
            PLOT_EXP_LINES = (
                "SNAPSHOTCY" in plotTypes or
                "AGGREGATEFC" in plotTypes or
                "AGGREGATEFC_DiffCI" in plotTypes or
                "SNAPSHOTFC" in plotTypes or
                "CalcGrossStats" in plotTypes )

            if PLOT_EXP_LINES:
                binVarsExpLines = {
                    (vu.obsVarQC, bu.goodQCMethod):{
                        'CalcGrossStats': True,
                        'OnlyPlotStats': statNames,
                    },
                    (vu.obsVarLat, bu.latbandsMethod):{
                        'CalcGrossStats': False,
                        'OnlyPlotStats': statNames,
                    },
                    (vu.obsVarPrs, bu.PjetMethod):{
                        'CalcGrossStats': False,
                        'OnlyPlotStats': statNames,
                    },
                    (vu.obsVarAlt, bu.altjetMethod):{
                        'CalcGrossStats': False,
                        'OnlyPlotStats': statNames,
                    },
                    (vu.obsVarCldFrac, bu.cloudbandsMethod):{
                        'CalcGrossStats': True,
                        'OnlyPlotStats': statNames,
                    },
                    (vu.obsVarLandFrac, bu.surfbandsMethod):{
                        'CalcGrossStats': False,
                        'OnlyPlotStats': statNames,
                    },
                }

                for (selectBinVar,binMethod), options in binVarsExpLines.items():
                    binVar = vu.varDictObs[selectBinVar][1]
                    if (binVar not in diagBinVars or
                        binMethod not in diagBinMethods): continue

                    selectStatNames = options.get('OnlyPlotStats',statNames)

                    binLoc = {}
                    binLoc['diagName'] = diagName
                    binLoc['binVar'] = binVar
                    binVarLoc = deepcopy(binLoc)
                    binLoc['binMethod'] = binMethod

                    dbSelect1DBinVals = diagDFW.levels('binVal',binLoc)
                    binUnitss = diagDFW.uniquevals('binUnits',binLoc)
                    if (len(binUnitss) == 0 or
                        len(dbSelect1DBinVals) < 1): continue

                    binUnits = binUnitss[0]

                    binMethodFile = ''
                    if binMethod != bu.goodQCMethod: binMethodFile = binMethod+'_'

                    # reorder select1DBinVals to match binMethod definition
                    tmp = deepcopy(bcs.binVarConfigs.get(
                        selectBinVar,{}).get(
                        binMethod,{}).get(
                        'values',dbSelect1DBinVals))
                    select1DBinVals = []
                    if (not isinstance(tmp, Iterable) or
                        isinstance(tmp,str)):
                        select1DBinVals += [tmp]
                    else:
                        select1DBinVals += tmp
                    for Bin in dbSelect1DBinVals:
                        if Bin not in select1DBinVals:
                            select1DBinVals.append(Bin)
                    for Bin in list(select1DBinVals):
                        if Bin not in dbSelect1DBinVals:
                            select1DBinVals.remove(Bin)

                    binTitles = []
                    for binVal in select1DBinVals:
                        if pu.isfloat(binVal) or pu.isint(binVal):
                            t = " @ "+binVar+"="+binVal
                            if binUnits != vu.miss_s:
                                t = t+" "+binUnits
                        elif binVal in [bcs.goodFlagName]:
                            t = ""
                        else:
                            t = " @ "+binVal
                        binTitles.append(t)

                    binMap = list(zip(select1DBinVals,binTitles))

                    # subplot configuration
                    if len(select1DBinVals) > 1:
                        nsubplots = len(select1DBinVals)*nVars
                        nx = len(select1DBinVals)
                        ny = nVars
                    else:
                        nsubplots = nVars
                        nx = np.int(np.ceil(np.sqrt(nsubplots)))
                        ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

                    binMessage = ' across '+binVar+' and binMethod=>'+binMethod

                    if ("SNAPSHOTCY" in plotTypes and
                        nCY > 1):
                        print("\nGenerating SNAPSHOTCY figures"+binMessage)
                        ptypePath = OSFigPath/'SNAPSHOTCY'/diagName
                        subplot_size = 1.9
                        aspect = 0.75

                        lineLoc = deepcopy(binLoc)
                        axLimLoc = deepcopy(binVarLoc)
                        # axLimLoc['binVal'] = select1DBinVals

                        #file loop 1
                        for (fcTDelta, fcTDelta_dir) in fcMap:
                            lineLoc['fcTDelta'] = fcTDelta

                            #file loop 2
                            for (statName, statDiagLabel, sciTicks, signDefinite) in StatsMap:
                                if statName not in selectStatNames: continue

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                                iplot = 0

                                #subplot loop 1
                                for (varName, varLabel) in varMap:
                                    lineLoc['varName'] = varName

                                    # use specific y-axis limits for each varName
                                    axLimLoc['varName'] = varName
                                    if statName == 'Count':
                                        dmin = 0.
                                    else:
                                        dmin = diagDFW.min(axLimLoc,statName)
                                    dmax = diagDFW.max(axLimLoc,statName)

                                    #subplot loop 2
                                    for binVal, binTitle in binMap:
                                        lineLoc['binVal'] = binVal

                                        # collect statName for all lines on this subplot
                                        linesVals = []
                                        for expName in expNames:
                                            lineLoc['expName'] = expName
                                            lineCYDTimes = diagDFW.levels('cyDTime',lineLoc)

                                            lineVals = np.full(nCY,np.NaN)
                                            cyLoc = deepcopy(lineLoc)
                                            for cyDTime in lineCYDTimes:
                                                icy = cyDTimes.index(cyDTime)
                                                cyLoc['cyDTime'] = cyDTime
                                                lineVals[icy] = diagDFW.loc(cyLoc,statName)
                                            linesVals.append(lineVals)

                                        # define subplot title
                                        title = varLabel+binTitle

                                        # perform subplot agnostic plotting (all expNames)
                                        bpf.plotTimeSeries(
                                            fig,
                                            cyDTimes, linesVals, expNames,
                                            title, statDiagLabel,
                                            sciTicks, signDefinite,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)

                                        iplot = iplot + 1

                                    # end statMap loop

                                # end varMap loop

                                # save each figure
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%sTSeries_%sday_%s_%s_%s'%(
                                           binMethodFile,fcTDelta_dir,DiagSpaceName,
                                           diagName,statName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                            # end statName loop

                        # end fcMap loop

                    # end SNAPSHOTCY


                    if ("AGGREGATEFC" in plotTypes and
                        nFC > 1):
                        print("\nGenerating AGGREGATEFC figures"+binMessage)
                        ptypePath = OSFigPath/'AGGREGATEFC'/fcDiagName
                        subplot_size = 1.9
                        aspect = 0.6

                        lineLoc = deepcopy(binLoc)
                        axLimLoc = deepcopy(binVarLoc)
                        # axLimLoc['binVal'] = select1DBinVals

                        #file loop 1
                        for (statName, statDiagLabel, sciTicks, signDefinite) in FCStatsMap:
                            if statName not in selectStatNames: continue

                            # establish a new figure
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                            iplot = 0

                            #subplot loop 1
                            for (varName, varLabel) in varMap:
                                lineLoc['varName'] = varName

                                # use specific y-axis limits for each varName
                                axLimLoc['varName'] = varName
                                if statName == 'Count':
                                    dmin = 0.
                                else:
                                    dmin = aggCYDFW.min(axLimLoc,statName)
                                dmax = aggCYDFW.max(axLimLoc,statName)

                                #subplot loop 2
                                for binVal, binTitle in binMap:
                                    lineLoc['binVal'] = binVal

                                    #collect aggregated statNames, varying across fcTDelta
                                    linesVals = []
                                    for expName in expNames:
                                        lineLoc['expName'] = expName
                                        lineFCTDeltas = aggCYDFW.levels('fcTDelta',lineLoc)

                                        lineVals = np.full(nFC,np.NaN)
                                        fcLoc = deepcopy(lineLoc)
                                        for fcTDelta in lineFCTDeltas:
                                            ifc = fcTDeltas.index(fcTDelta)
                                            fcLoc['fcTDelta'] = fcTDelta
                                            lineVals[ifc] = aggCYDFW.loc(fcLoc,statName)
                                        linesVals.append(lineVals)

                                    # define subplot title
                                    title = varLabel+binTitle

                                    # perform subplot agnostic plotting (all expNames)
                                    bpf.plotTimeSeries(
                                        fig,
                                        fcTDeltas, linesVals, expNames,
                                        title, statDiagLabel,
                                        sciTicks, signDefinite,
                                        ny, nx, nsubplots, iplot,
                                        dmin = dmin, dmax = dmax,
                                        interiorLabels = interiorLabels)

                                    iplot = iplot + 1

                                # end statMap loop

                            # end varMap loop

                            # save each figure
                            ptypePath.mkdir(parents=True, exist_ok=True)
                            filename = ptypePath/('%sTSeries_%s-%sday_%s_%s_%s'%(
                                       binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                       fcDiagName,statName))

                            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                        # end statName loop

                    # end AGGREGATEFC


                    if ("AGGREGATEFC_DiffCI" in plotTypes and
                        nFC > 1 and nExp > 1):
                        print("\nGenerating AGGREGATEFC_DiffCI figures"+binMessage)
                        ptypePath = OSFigPath/'AGGREGATEFC_DiffCI'/fcDiagName
                        subplot_size = 1.9
                        aspect = 0.6

                        # nx = 1
                        # ny = nVars
                        # nsubplots = nx * ny

                        commonLoc = deepcopy(binLoc)

                        #figure loop 2
                        for statName in bootStrapStats:
                            if statName not in selectStatNames: continue

                            # establish a new figure
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                            iplot = 0

                            #subplot loop 1
                            for (varName, varLabel) in varMap:
                                commonLoc['varName'] = varName

                                #subplot loop 2
                                for binVal, binTitle in binMap:
                                    commonLoc['binVal'] = binVal

                                    # define subplot title
                                    title = varLabel+binTitle

                                    linesVals = {}
                                    for trait in su.ciTraits: linesVals[trait] = []

                                    cntrlLoc = deepcopy(commonLoc)
                                    cntrlLoc['expName'] = cntrlName
                                    expLoc = deepcopy(commonLoc)
                                    for expName in noncntrlExpNames:
                                        expLoc['expName'] = expName

                                        lineVals = {}
                                        for trait in su.ciTraits: lineVals[trait] = []

                                        for fcTDelta in fcTDeltas:
                                            cntrlLoc['fcTDelta'] = fcTDelta
                                            expLoc['fcTDelta'] = fcTDelta

                                            ciVals = su.bootStrapClusterFunc(
                                                         X = diagDFW.loc(expLoc),
                                                         Y = diagDFW.loc(cntrlLoc),
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
                                        ny, nx, nsubplots, iplot,
                                        linesValsMinCI = linesVals[su.cimin],
                                        linesValsMaxCI = linesVals[su.cimax],
                                        dmin = dmin, dmax = dmax,
                                        lineAttribOffset = 1,
                                        interiorLabels = interiorLabels)
                                    iplot = iplot + 1

                                # end binMap loop

                            # end varName loop

                            # save each figure
                            ptypePath.mkdir(parents=True, exist_ok=True)
                            filename = ptypePath/('%sTSeries_%s-%sday_%s_%s_%s'%(
                                       binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                       fcDiagName,statName))

                            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)
                        # end statName loop

                    # end AGGREGATEFC_DiffCI

                    if ("SNAPSHOTFC" in plotTypes and
                        nFC > 1 and nCY > 1):
                        print("\nGenerating SNAPSHOTFC figures"+binMessage)
                        ptypePath = OSFigPath/'SNAPSHOTFC'/fcDiagName
                        subplot_size = 1.9
                        aspect = 0.75

                        lineLoc = deepcopy(binLoc)
                        axLimLoc = deepcopy(binVarLoc)
                        # axLimLoc['binVal'] = select1DBinVals

                        #file loop 1
                        for expName in expNames:
                            lineLoc['expName'] = expName

                            #file loop 2
                            for (statName, statDiagLabel, sciTicks, signDefinite) in StatsMap:
                                if statName not in selectStatNames: continue

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                                iplot = 0

                                #subplot loop 1
                                for (varName, varLabel) in varMap:
                                    lineLoc['varName'] = varName

                                    # use specific y-axis limits for each varName
                                    axLimLoc['varName'] = varName
                                    if statName == 'Count':
                                        dmin = 0.
                                    else:
                                        dmin = diagDFW.min(axLimLoc,statName)
                                    dmax = diagDFW.max(axLimLoc,statName)

                                    #subplot loop 2
                                    for binVal, binTitle in binMap:
                                        lineLoc['binVal'] = binVal

                                        # collect statName for all lines on this subplot, letting cyDTime vary
                                        xsVals = []
                                        linesVals = []
                                        fcTDeltas_labels = []
                                        for fcTDelta in fcTDeltas:
                                            lineLoc['fcTDelta'] = fcTDelta

                                            xVals = []
                                            for cyDTime in cyDTimes:
                                                xVals.append(cyDTime+fcTDelta)
                                            xsVals.append(xVals)

                                            #Setting to avoid over-crowding
                                            if fcTDeltas.index(fcTDelta) > (MAX_FC_LINES-1): continue

                                            fcTDelta_sec = pu.TDeltas2Seconds([fcTDelta])
                                            fcTDeltas_labels.append(pu.timeTicks(fcTDelta_sec[0],0))

                                            lineCYDTimes = diagDFW.levels('cyDTime',lineLoc)

                                            lineVals = np.full(nCY,np.NaN)
                                            cyLoc = deepcopy(lineLoc)
                                            for cyDTime in lineCYDTimes:
                                                icy = cyDTimes.index(cyDTime)
                                                cyLoc['cyDTime'] = cyDTime
                                                lineVals[icy] = diagDFW.loc(cyLoc,statName)

                                            linesVals.append(lineVals)

                                        # define subplot title
                                        title = varLabel+binTitle

                                        # perform subplot agnostic plotting (all expNames)
                                        bpf.plotTimeSeries(
                                            fig,
                                            xsVals, linesVals, fcTDeltas_labels,
                                            title, statDiagLabel,
                                            sciTicks, signDefinite,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)
                                        iplot = iplot + 1

                                    # end binMap loop

                                # end varMap loop

                                expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%sTSeries_%s_%s_%s_%s'%(
                                           binMethodFile,expName_file,DiagSpaceName,
                                           fcDiagName,statName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                            # end statMap loop

                        # end expName loop

                    # end SNAPSHOTFC

                    CalcGrossStats = options.get('CalcGrossStats',False)
                    if ('CalcGrossStats' in plotTypes and CalcGrossStats):
                        print("\nCalculating gross statistics"+binMessage)
                        print(" at FC length ",fcTDeltas[0])
                        # Calculate gross statistics for this binVal
                        statsLoc = deepcopy(binLoc)
                        statsLoc['fcTDelta'] = fcTDeltas[0]
                        for binVal in select1DBinVals:
                            statsLoc['binVal'] = binVal
                            GrossStats = {}
                            for varName in varNames:
                                statsLoc['varName'] = varName
                                for expName in expNames:
                                    statsLoc['expName'] = expName
                                    statsDFW = DFWrapper.fromLoc(aggCYDFW,statsLoc)
                                    for statName in statNames:
                                        GrossStats[(statName,expName,varName)] = statsDFW.var(statName).to_numpy()
                            for expName in expNames:
                                print("Gross statistics for")
                                print("experiment=>"+expName)
                                if len(select1DBinVals) > 1:
                                    print("binVal=>"+binVal)
                                print(" variables: ",varNames)

                                for statName in statNames:
                                    print(statName)
                                    tmp = np.asarray([])
                                    for varName in varNames:
                                        tmp = np.append(tmp,GrossStats[(statName,expName,varName)])
                                    print(tmp)
                    # end CalcGrossStats

                # end binVar+binMethod loop

            # end PLOT_EXP_LINES


            #TODO: LatLines and QCLines both have binVal lines.  Combine these into one loop.
            #      differentiate from above ExpLines plots that have expName lines.
            #      QCLines only plots statName=='Count', so must add selectStatNames
            #      dictionary option similar to 2D figures or plot all statNames...
            selectBinVar = vu.varDictObs[vu.obsVarLat][1]
            select1DBinVals = bcs.namedLatBands['values']
            selectBinMethod = bu.latbandsMethod
            if ("SNAPSHOTCY_LatLines" in plotTypes and
                selectBinVar in diagBinVars and
                len(select1DBinVals) > 0 and
                nCY > 1 and
                statsDB.DiagSpaceGrp[0] != conf.radiance_s):

                print("\nGenerating SNAPSHOTCY_LatLines figures...")
                ptypePath = OSFigPath/'SNAPSHOTCY_LatLines'/diagName

                subplot_size = 1.9
                aspect = 0.75

                nsubplots = nExp * nVars
                nx = nExp
                ny = nVars

                axLimLoc = {}
                axLimLoc['diagName'] = diagName
                axLimLoc['binVar'] = selectBinVar
                axLimLoc['binMethod'] = selectBinMethod
                axLimLoc['binVal'] = select1DBinVals

                lineLoc = deepcopy(axLimLoc)
                #file loop 1
                for (fcTDelta, fcTDelta_dir) in fcMap:
                    lineLoc['fcTDelta'] = fcTDelta

                    #file loop 2
                    for (statName, statDiagLabel, sciTicks, signDefinite) in StatsMap:

                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        iplot = 0

                        #subplot loop 1
                        for (varName, varLabel) in varMap:
                            lineLoc['varName'] = varName

                            # use specific y-axis limits for each varName
                            axLimLoc['varName'] = varName
                            if statName == 'Count':
                                dmin = 0.
                            else:
                                dmin = diagDFW.min(axLimLoc,statName)
                            dmax = diagDFW.max(axLimLoc,statName)

                            #subplot loop 2
                            for expName in expNames:
                                lineLoc['expName'] = expName

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for binVal in select1DBinVals:
                                    lineLoc['binVal'] = binVal
                                    lineCYDTimes = diagDFW.levels('cyDTime',lineLoc)

                                    lineVals = np.full(nCY,np.NaN)
                                    cyLoc = deepcopy(lineLoc)
                                    for cyDTime in lineCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        cyLoc['cyDTime'] = cyDTime
                                        lineVals[icy] = diagDFW.loc(cyLoc,statName)

                                    linesVals.append(lineVals)

                                # define subplot title
                                title = expName+"\n"+varLabel

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotTimeSeries(
                                    fig,
                                    cyDTimes, linesVals, select1DBinVals,
                                    title, statDiagLabel,
                                    sciTicks, signDefinite,
                                    ny, nx, nsubplots, iplot,
                                    dmin = dmin, dmax = dmax,
                                    interiorLabels = interiorLabels)

                                iplot = iplot + 1
                        ptypePath.mkdir(parents=True, exist_ok=True)
                        filename = ptypePath/('TSeries_%sday_%s_%s_%s'%(
                                   fcTDelta_dir,DiagSpaceName,
                                   diagName,statName))

                        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                # end binVal loop

            # end SNAPSHOTCY_LatLines

            selectBinVar = vu.varDictObs[vu.obsVarQC][1]
            select1DBinVals = [bcs.goodFlagName] + bcs.badFlagNames
            selectBinMethods = [bu.goodQCMethod,bu.badQCMethod]
            if ("SNAPSHOTCY_QCLines" in plotTypes and
                selectBinVar in diagBinVars and
                len(select1DBinVals) > 0 and
                nCY > 1):
                print("\nGenerating SNAPSHOTCY_QCLines figures...")
                ptypePath = OSFigPath/'SNAPSHOTCY_QCLines'/diagName
                subplot_size = 1.9
                aspect = 0.75

                nsubplots = nExp * nVars
                nx = nExp
                ny = nVars

                statName = 'Count'
                istat = statNames.index(statName)
                (statName, statDiagLabel, sciTicks, signDefinite) = StatsMap[istat]

                axLimLoc = {}
                axLimLoc['diagName'] = diagName
                axLimLoc['binVar'] = selectBinVar
                # axLimLoc['binMethod'] = selectBinMethods
                # axLimLoc['binVal'] = select1DBinVals

                lineLoc = deepcopy(axLimLoc)
                lineLoc['binMethod'] = selectBinMethods
                lineLoc['binVal'] = select1DBinVals

                #file loop 1
                for (fcTDelta, fcTDelta_dir) in fcMap:
                    lineLoc['fcTDelta'] = fcTDelta

                    # establish a new figure
                    fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                    iplot = 0

                    #subplot loop 1
                    for (varName, varLabel) in varMap:
                        lineLoc['varName'] = varName

                        # use specific y-axis limits for each varName
                        axLimLoc['varName'] = varName
                        if statName == 'Count':
                            dmin = 0.
                        else:
                            dmin = diagDFW.min(axLimLoc,statName)
                        dmax = diagDFW.max(axLimLoc,statName)

                        #subplot loop 2
                        for expName in expNames:
                            lineLoc['expName'] = expName

                            # collect statName for all lines on this subplot, letting cyDTime vary
                            linesVals = []
                            for binVal in select1DBinVals:
                                lineLoc['binVal'] = binVal
                                lineCYDTimes = diagDFW.levels('cyDTime',lineLoc)

                                lineVals = np.full(nCY,np.NaN)
                                cyLoc = deepcopy(lineLoc)
                                for cyDTime in lineCYDTimes:
                                    icy = cyDTimes.index(cyDTime)
                                    cyLoc['cyDTime'] = cyDTime
                                    lineVals[icy] = diagDFW.loc(cyLoc,statName)
                                linesVals.append(lineVals)

                            # define subplot title
                            title = expName+"\n"+varLabel

                            # perform subplot agnostic plotting (all expNames)
                            bpf.plotTimeSeries(
                                fig,
                                cyDTimes, linesVals, select1DBinVals,
                                title, statDiagLabel,
                                sciTicks, signDefinite,
                                ny, nx, nsubplots, iplot,
                                dmin = dmin, dmax = dmax,
                                legend_inside = False,
                                interiorLabels = interiorLabels)

                            iplot = iplot + 1

                    ptypePath.mkdir(parents=True, exist_ok=True)
                    filename = ptypePath/('TSeries_%sday_%s_%s_%s'%(
                               fcTDelta_dir,DiagSpaceName,
                               diagName,statName))

                    pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                # end binVal loop

            # end SNAPSHOTCY_QCLines


            #=============
            # 2-D figures
            #=============
            PLOT_ACROSS_BINS = (
                "SNAPSHOTCY2D" in plotTypes or
                "AGGREGATEFC2D" in plotTypes or
                "AGGREGATEFC_Profile" in plotTypes or
                "AGGREGATEFC_Profile_DiffCI" in plotTypes)

            if PLOT_ACROSS_BINS:
                binVars2D = {
                    vu.obsVarAlt: {'profilePlotFunc': bpf.plotProfile},
                    vu.obsVarACI: {'profilePlotFunc': bpf.plotSeries},
                    vu.obsVarCldFrac: {'profilePlotFunc': bpf.plotSeries},
                    vu.obsVarGlint: {'profilePlotFunc': bpf.plotSeries},
                    vu.obsVarLandFrac: {'profilePlotFunc': bpf.plotSeries},
                    vu.obsVarLat: {'profilePlotFunc': bpf.plotProfile},
                    vu.obsVarLT: {'profilePlotFunc': bpf.plotSeries},
                    vu.obsVarPrs: {'profilePlotFunc': bpf.plotProfile},
                    vu.obsVarSenZen: {'profilePlotFunc': bpf.plotSeries},
                }

                for selectBinVar, options in binVars2D.items():
                    binVar = vu.varDictObs[selectBinVar][1]
                    if (binVar not in diagBinVars): continue


                    #Get all float/int binVals associated with binVar
                    binVarLoc = {}
                    binVarLoc['diagName'] = diagName
                    binVarLoc['binVar'] = binVar
                    binVarLoc['binVal'] = binNumVals2DasStr

                    selectBinMethods = diagDFW.levels('binMethod',binVarLoc)

                    selectStatNames = options.get('stats',statNames)

                    varMethodLoc = deepcopy(binVarLoc)
                    axLimLoc = {}
                    axLimLoc['diagName'] = diagName
                    axLimLoc['binVal'] = binNumVals2DasStr

                    ## all binMethod's with same binVar have same axis limits
                    #  comment out to share limits between binVars
                    axLimLoc['binVar'] = binVar

                    for binMethod in selectBinMethods:
                        ## all binVar's with same binMethod have same axis limits
                        #  comment out to share limits between binMethods
                        axLimLoc['binMethod'] = binMethod

                        varMethodLoc['binMethod'] = binMethod
                        select2DBinVals = diagDFW.levels('binVal',varMethodLoc)

                        binUnits = diagDFW.uniquevals('binUnits',varMethodLoc)[0]

                        # assume all bins represent same variable/units
                        indepLabel = binVar
                        if binUnits != vu.miss_s:
                            indepLabel = indepLabel+" ("+binUnits+")"

                        binMethodFile = ''
                        if binMethod != bu.identityBinMethod: binMethodFile = '_'+binMethod

                        # bin info
                        select2DBinNumVals = []
                        for binVal in select2DBinVals:
                            ibin = allBinVals.index(binVal)
                            select2DBinNumVals.append(binNumVals[ibin])

                            # invert independent variable axis for pressure bins
                            pressure_dict = vu.varDictObs.get(vu.obsVarPrs,['',''])
                            invert_ind_axis = (pressure_dict[1] == binVar)

                        # sort bins by numeric value
                        indices = list(range(len(select2DBinNumVals)))
                        indices.sort(key=select2DBinNumVals.__getitem__)
                        select2DBinNumVals = list(map(select2DBinNumVals.__getitem__, indices))
                        select2DBinVals = list(map(select2DBinVals.__getitem__, indices))

                        if len(select2DBinVals) < 2: continue

                        binMessage = ' across '+binVar+' and binMethod=>'+binMethod

                        if ("SNAPSHOTCY2D" in plotTypes and
                            nCY > 1):

                            print("\nGenerating SNAPSHOTCY2D figures"+binMessage)
                            ptypePath = OSFigPath/'SNAPSHOTCY2D'/diagName
                            # top-level 2-D plot settings
                            subplot_size = 2.4
                            aspect = 0.65
                            nx = nExp
                            ny = nVars
                            nsubplots = nx * ny

                            planeLoc = deepcopy(varMethodLoc)

                            #file loop 1
                            for (fcTDelta, fcTDelta_dir) in fcMap:
                                planeLoc['fcTDelta'] = fcTDelta

                                #file loop 2
                                for (statName, statDiagLabel, sciTicks, signDefinite) in StatsMap:

                                    if statName not in selectStatNames: continue

                                    # establish a new figure
                                    fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                    iplot = 0

                                    #subplot loop 1
                                    for (varName, varLabel) in varMap:
                                        planeLoc['varName'] = varName

                                        # use specific c-axis limits for each varName
                                        axLimLoc['varName'] = varName
                                        if statName == 'Count':
                                            dmin = 0.
                                        else:
                                            dmin = diagDFW.min(axLimLoc,statName)
                                        dmax = diagDFW.max(axLimLoc,statName)

                                        #subplot loop 2
                                        # letting cyDTime and binVal vary
                                        for expName in expNames:
                                            planeLoc['expName'] = expName
                                            planeCYDTimes = diagDFW.levels('cyDTime',planeLoc)

                                            planeVals = np.full((len(select2DBinVals),nCY),np.NaN)
                                            binLoc = deepcopy(planeLoc)
                                            for ibin, binVal in enumerate(select2DBinVals):
                                                binLoc['binVal'] = binVal
                                                tmp = diagDFW.loc(binLoc,statName).to_numpy()
                                                for jcy, cyDTime in enumerate(planeCYDTimes):
                                                    icy = cyDTimes.index(cyDTime)
                                                    planeVals[ibin,icy] = tmp[jcy]

                                            # define subplot title
                                            title = expName+"\n"+varLabel

                                            # perform subplot agnostic plotting (all expNames)
                                            bpf.plotTimeSeries2D(
                                                fig,
                                                cyDTimes, select2DBinNumVals, planeVals,
                                                title, statDiagLabel,
                                                sciTicks, signDefinite,
                                                indepLabel, invert_ind_axis,
                                                ny, nx, nsubplots, iplot,
                                                dmin = dmin, dmax = dmax,
                                                interiorLabels = interiorLabels)

                                            iplot = iplot + 1

                                    ptypePath.mkdir(parents=True, exist_ok=True)
                                    filename = ptypePath/('%s%s_TSeries_%sday_%s_%s_%s'%(
                                               binVar,binMethodFile,fcTDelta_dir,DiagSpaceName,
                                               diagName,statName))

                                    pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                                # end statName loop

                            # end fcTDelta loop

                        # end SNAPSHOTCY2D

                        if ("AGGREGATEFC2D" in plotTypes and
                            nFC > 1):

                            print("\nGenerating AGGREGATEFC2D figures"+binMessage)
                            ptypePath = OSFigPath/'AGGREGATEFC2D'/fcDiagName
                            # top-level 2-D plot settings
                            subplot_size = 2.4
                            aspect = 0.55
                            nx = nExp
                            ny = nVars
                            nsubplots = nx * ny

                            planeLoc = deepcopy(varMethodLoc)

                            #file loop 1
                            for (statName, statDiagLabel, sciTicks, signDefinite) in FCStatsMap:
                                if statName not in selectStatNames: continue

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                                iplot = 0
                                #subplot loop
                                for (varName, varLabel) in varMap:
                                    planeLoc['varName'] = varName

                                    # use specific c-axis limits for each varName
                                    axLimLoc['varName'] = varName
                                    if statName == 'Count':
                                        dmin = 0.
                                    else:
                                        dmin = aggCYDFW.min(axLimLoc,statName)
                                    dmax = aggCYDFW.max(axLimLoc,statName)

                                    #collect aggregated statName, varying across fcTDelta+binVal
                                    for expName in expNames:
                                        planeLoc['expName'] = expName
                                        planeFCTDeltas = aggCYDFW.levels('fcTDelta',planeLoc)

                                        planeVals = np.full((len(select2DBinVals),nFC),np.NaN)
                                        binLoc = deepcopy(planeLoc)
                                        for ibin, binVal in enumerate(select2DBinVals):
                                            binLoc['binVal'] = binVal
                                            tmp = aggCYDFW.loc(binLoc,statName).to_numpy()
                                            for jfc, fcTDelta in enumerate(planeFCTDeltas):
                                                ifc = fcTDeltas.index(fcTDelta)
                                                planeVals[ibin,ifc] = tmp[jfc]

                                        # define subplot title
                                        title = expName+"\n"+varLabel

                                        # perform subplot agnostic plotting (all expNames)
                                        bpf.plotTimeSeries2D(
                                            fig,
                                            fcTDeltas, select2DBinNumVals, planeVals,
                                            title, statDiagLabel,
                                            sciTicks, signDefinite,
                                            indepLabel, invert_ind_axis,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)

                                        iplot = iplot + 1

                                # save each figure
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%s%s_TSeries_%s-%sday_%s_%s_%s'%(
                                           binVar,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

                            # end statName loop

                        # end AGGREGATEFC2D

                        if "AGGREGATEFC_Profile" in plotTypes:
                            print("\nGenerating AGGREGATEFC_Profile figures"+binMessage)
                            ptypePath = OSFigPath/'AGGREGATEFC_Profile'/fcDiagName
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

                            ptLoc = deepcopy(varMethodLoc)

                            #file loop 1
                            for (statName, statDiagLabel, sciTicks, signDefinite) in FCStatsMap:
                                if statName not in selectStatNames: continue

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                                iplot = 0

                                #subplot loop 1
                                for (varName, varLabel) in varMap:
                                    ptLoc['varName'] = varName

                                    # use specific x-axis limits for each varName
                                    axLimLoc['varName'] = varName
                                    if statName == 'Count':
                                        dmin = 0.
                                    else:
                                        dmin = aggCYDFW.min(axLimLoc,statName)
                                    dmax = aggCYDFW.max(axLimLoc,statName)

                                    #subplot loop 2
                                    for fcTDelta in fcTDeltas:
                                        ptLoc['fcTDelta'] = fcTDelta

                                        #Setting to avoid over-crowding
                                        if fcTDeltas.index(fcTDelta) > (MAX_FC_SUBFIGS-1): continue

                                        #collect aggregated statNames, varying across fcTDelta
                                        linesVals = []
                                        for expName in expNames:
                                            ptLoc['expName'] = expName

                                            lineVals = []
                                            for binVal in select2DBinVals:
                                                ptLoc['binVal'] = binVal
                                                lineVals.append(aggCYDFW.loc(ptLoc,statName))

                                            linesVals.append(lineVals)

                                        # define subplot title
                                        title = varLabel+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                        # perform subplot agnostic plotting (all expNames)
                                        options['profilePlotFunc'](
                                            fig,
                                            linesVals, select2DBinNumVals,
                                            expNames,
                                            title, statDiagLabel,
                                            sciTicks, signDefinite,
                                            indepLabel, invert_ind_axis,
                                            ny, nx, nsubplots, iplot,
                                            dmin = dmin, dmax = dmax,
                                            interiorLabels = interiorLabels)


                                        iplot = iplot + 1

                                # save each figure
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%s%s_Agg_%s-%sday_%s_%s_%s'%(
                                           binVar,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                            # end statName loop

                        # end AGGREGATEFC_Profile

                        if ("AGGREGATEFC_Profile_DiffCI" in plotTypes and
                            nExp > 1):
                            print("\nGenerating AGGREGATEFC_Profile_DiffCI figures"+binMessage)
                            ptypePath = OSFigPath/'AGGREGATEFC_Profile_DiffCI'/fcDiagName
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
                                if statName not in selectStatNames: continue

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                                iplot = 0

                                #subplot loop 1
                                fcLoc = deepcopy(varMethodLoc)
                                for (varName, varLabel) in varMap:
                                    fcLoc['varName'] = varName

                                    #subplot loop 2
                                    for fcTDelta in fcTDeltas:
                                        fcLoc['fcTDelta'] = fcTDelta
                                        cntrlLoc = deepcopy(fcLoc)
                                        expLoc = deepcopy(fcLoc)

                                        #Setting to avoid over-crowding
                                        if fcTDeltas.index(fcTDelta) > (MAX_FC_SUBFIGS-1): continue

                                        linesVals = {}
                                        for trait in su.ciTraits: linesVals[trait] = []
                                        for expName in noncntrlExpNames:
                                            lineVals = {}
                                            for trait in su.ciTraits: lineVals[trait] = []

                                            cntrlLoc['expName'] = cntrlName
                                            expLoc['expName'] = expName

                                            for binVal in select2DBinVals:
                                                cntrlLoc['binVal'] = binVal
                                                expLoc['binVal'] = binVal

                                                ciVals = su.bootStrapClusterFunc(
                                                             X = diagDFW.loc(expLoc),
                                                             Y = diagDFW.loc(cntrlLoc),
                                                             n_samples = 10000,
                                                             statNames = [statName])

                                                for trait in su.ciTraits:
                                                    lineVals[trait].append(ciVals[statName][trait][0])

                                            for trait in su.ciTraits:
                                                linesVals[trait].append(lineVals[trait])

                                        # define subplot title
                                        title = varLabel+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                        # use specific y-axis limits for each varName
                                        dmin = np.nanmin(linesVals[su.cimin])
                                        dmax = np.nanmax(linesVals[su.cimax])

                                        # perform subplot agnostic plotting (all expNames)
                                        options['profilePlotFunc'](
                                            fig,
                                            linesVals[su.cimean], select2DBinNumVals,
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
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%s%s_Agg_%s-%sday_%s_%s_%s'%(
                                           binVar,binMethodFile,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                           fcDiagName,statName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                            # end statName loop

                        # end AGGREGATEFC_Profile_DiffCI

                    # end binMethods loop

                # end binVars2D loop

            # end PLOT_ACROSS_BINS


            if ("AGGREGATEFC_PDF" in plotTypes):
                binVarsPDF = [
                    vu.varDictObs[vu.obsVarNormErr][1]
                ]
                for binVarPDF in binVarsPDF:
                    if binVarPDF not in diagBinVars: continue

                    print("\nGenerating AGGREGATEFC_PDF figures across "+binVarPDF)
                    ptypePath = OSFigPath/'AGGREGATEFC_PDF'/fcDiagName
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

                    #Get all float/int binVals associated with binVarPDF
                    binVarLoc = {}
                    binVarLoc['diagName'] = diagName
                    binVarLoc['binVar'] = binVarPDF
                    binVarLoc['binVal'] = binNumVals2DasStr

                    selectBinMethods = diagDFW.levels('binMethod',binVarLoc)

                    selectPDFBinVals = diagDFW.levels('binVal',binVarLoc)

                    binUnits = diagDFW.uniquevals('binUnits',binVarLoc)[0]

                    # assume all bins represent same variable/units
                    indepLabel = binVarPDF
                    if binUnits != vu.miss_s:
                        indepLabel = indepLabel+" ("+binUnits+")"

                    # bin info
                    selectPDFBinNumVals = []
                    for binVal in selectPDFBinVals:
                        ibin = allBinVals.index(binVal)
                        selectPDFBinNumVals.append(binNumVals[ibin])

                    # sort bins by numeric value
                    indices = list(range(len(selectPDFBinNumVals)))
                    indices.sort(key=selectPDFBinNumVals.__getitem__)
                    selectPDFBinNumVals = list(map(selectPDFBinNumVals.__getitem__, indices))
                    selectPDFBinVals = list(map(selectPDFBinVals.__getitem__, indices))

                    if len(selectPDFBinVals) < 2: continue

                    ptLoc = deepcopy(binVarLoc)
                    #file loop 1
                    for expName in expNames:
                        ptLoc['expName'] = expName

                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)

                        iplot = 0

                        #subplot loop 1
                        for (varName, varLabel) in varMap:
                            ptLoc['varName'] = varName

                            #subplot loop 2
                            for fcTDelta in fcTDeltas:
                                ptLoc['fcTDelta'] = fcTDelta

                                #Setting to avoid over-crowding
                                if fcTDeltas.index(fcTDelta) > (MAX_FC_SUBFIGS-1): continue

                                #collect aggregated statNames, varying across fcTDelta
                                linesVals = []
                                binMethodLabels = []
                                for binMethod in selectBinMethods:
                                    ptLoc['binMethod'] = binMethod

                                    # if binMethod != bu.identityBinMethod: do something with bu.identityBinMethod
                                    if binMethod == bu.identityBinMethod:
                                        binMethodLabels.append('CONST')
                                    else:
                                        binMethodLabels.append(binMethod)

                                    lineVals = []
                                    for binVal in selectPDFBinVals:
                                        ptLoc['binVal'] = binVal
                                        lineVals.append(aggCYDFW.loc(ptLoc,'Count'))

                                    linesVals.append(lineVals)

                                # define subplot title
                                title = varLabel+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                # perform subplot agnostic plotting (all expNames)
                                bpf.plotPDF(
                                    fig,
                                    linesVals, selectPDFBinNumVals,
                                    binMethodLabels,
                                    title,
                                    indepLabel,
                                    ny, nx, nsubplots, iplot,
                                    interiorLabels = interiorLabels)

                                iplot = iplot + 1

                        # save each figure
                        ptypePath.mkdir(parents=True, exist_ok=True)
                        filename = ptypePath/('%s_Agg_%s-%sday_%s_%s_%s'%(
                                   binVarPDF,fcTDeltas_dir[0],fcTDeltas_dir[-1],DiagSpaceName,
                                   fcDiagName,expName))

                        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                    # end expName loop

            # end AGGREGATEFC_PDF


            if "AGGREGATEFC_StatsComposite" in plotTypes:
                binVarsStats = {
                    vu.obsVarSCI: {'statsPlotFunc': bpf.plotfitRampComposite},
                }

                for selectBinVar, options in binVarsStats.items():
                    statBinVar = vu.varDictObs[selectBinVar][1]
                    if (statBinVar not in diagBinVars): continue

                    print("\nGenerating AGGREGATEFC_StatsComposite figures across "+statBinVar)
                    ptypePath = OSFigPath/'AGGREGATEFC_StatsComposite'/fcDiagName
                    # top-level 2-D plot settings
                    subplot_size = 1.9
                    aspect = 0.9

                    nsubplots = nVars
                    nx = np.int(np.ceil(np.sqrt(nsubplots)))
                    ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

                    #Get all float/int binVals associated with statBinVar
                    binVarLoc = {}
                    binVarLoc['diagName'] = diagName
                    binVarLoc['binVar'] = statBinVar
                    binVarLoc['binVal'] = binNumVals2DasStr

                    selectBinMethods = diagDFW.levels('binMethod',binVarLoc)

                    selectStatBinVals = diagDFW.levels('binVal',binVarLoc)

                    binUnits = diagDFW.uniquevals('binUnits',binVarLoc)[0]


                    # assume all bins represent same variable/units
                    indepLabel = statBinVar
                    if binUnits != vu.miss_s:
                        indepLabel = indepLabel+" ("+binUnits+")"

                    # bin info
                    selectStatBinNumVals = []
                    for binVal in selectStatBinVals:
                        ibin = allBinVals.index(binVal)
                        selectStatBinNumVals.append(binNumVals[ibin])


                    # sort bins by numeric value
                    indices = list(range(len(selectStatBinNumVals)))
                    indices.sort(key=selectStatBinNumVals.__getitem__)
                    selectStatBinNumVals = list(map(selectStatBinNumVals.__getitem__, indices))
                    selectStatBinVals = list(map(selectStatBinVals.__getitem__, indices))

                    nBins = len(selectStatBinVals)
                    if nBins < 2: continue

                    ptLoc = deepcopy(binVarLoc)

                    #file loop 1
                    for binMethod in selectBinMethods:
                        ptLoc['binMethod'] = binMethod

                        binMethodFile = ''
                        if binMethod != bu.identityBinMethod: binMethodFile = '_'+binMethod

                        #file loop 2
                        for expName in expNames:
                            ptLoc['expName'] = expName

                            #file loop 3
                            for (fcTDelta, fcTDelta_dir) in fcMap:
                                ptLoc['fcTDelta'] = fcTDelta

                                # establish a new figure
                                fig = pu.setup_fig(nx, ny, subplot_size, aspect, interiorLabels)
                                iplot = 0

                                ERRParams = {}
                                ERRParams[DiagSpaceName] = {}

                                #subplot loop 1
                                for (varName, varLabel) in varMap:
                                    ptLoc['varName'] = varName

                                    #collect aggregated statNames, varying across fcTDelta
                                    countsVals = np.full(nBins,0)
                                    meansVals  = np.full(nBins,np.NaN)
                                    rmssVals   = np.full(nBins,np.NaN)
                                    stdsVals   = np.full(nBins,np.NaN)

                                    for ibin, binVal in enumerate(selectStatBinVals):
                                        ptLoc['binVal'] = binVal
                                        countsVals[ibin] = aggCYDFW.loc(ptLoc,'Count').to_numpy()
                                        meansVals[ibin] = aggCYDFW.loc(ptLoc,'Mean').to_numpy()
                                        rmssVals[ibin] = aggCYDFW.loc(ptLoc,'RMS').to_numpy()
                                        stdsVals[ibin] = aggCYDFW.loc(ptLoc,'STD').to_numpy()

                                    # define subplot title
                                    title = varLabel

                                    # perform subplot agnostic plotting (all expNames)
                                    FitParams = options['statsPlotFunc'](
                                        fig,
                                        selectStatBinNumVals,
                                        countsVals,
                                        meansVals,
                                        rmssVals,
                                        stdsVals,
                                        title,
                                        "STATS("+fcDiagName+")",
                                        indepLabel,
                                        ny, nx, nsubplots, iplot,
                                        interiorLabels = interiorLabels)

                                    paramKey = statsDB.chlist[iplot]
                                    if paramKey == '': paramKey = varName
                                    ERRParams[DiagSpaceName][(paramKey,binMethod)] = FitParams

                                    iplot = iplot + 1

                                YAMLParams = {}
                                print("\n#For binning_params:")
                                for key in sorted(ERRParams[DiagSpaceName]):
                                    print(statBinVar+"ErrParams['"+DiagSpaceName+"'][",key,"]   = ",
                                           ERRParams[DiagSpaceName][key]['bu'])
                                    for param, val in ERRParams[DiagSpaceName][key]['YAML'].items():
                                        if param not in YAMLParams: YAMLParams[param] = []
                                        YAMLParams[param] += val
                                print("\n#For UFO YAML config:")
                                for param, val in YAMLParams.items():
                                    print('#  '+param+':',val)

                                # save each figure
                                ptypePath.mkdir(parents=True, exist_ok=True)
                                filename = ptypePath/('%s%s_Agg_%sday_%s_%s_%s'%(
                                           statBinVar,binMethodFile,fcTDelta_dir,DiagSpaceName,
                                           fcDiagName,expName))

                                pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                        # end expName loop

                    # end binMethod loop

                # end selectBinVar loop

            # end AGGREGATEFC_StatsComposite

        # end diagName loop

    # end DiagSpaceName loop

def main():
    nproc, myproc = paru.par_args(sys.argv[:])
    plot_stats_timeseries(nproc, myproc)

if __name__ == '__main__': main()
