import os, sys
from netCDF4 import Dataset
#import numpy
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import datetime as dt
import glob
import re
import plot_utils as pu
import math

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See plot_utils.par_args for more information.

#Select the type of plots
# options: 'SNAPSHOTCY', 'AGGREGATEFC','SNAPSHOTCY2D', 'AGGREGATEFC2D', 'SNAPSHOTFC'
# SNAPSHOTCY - creates a timeseries figure between firstCycleDTime and lastCycleDTime
#               for each forecast length between fcTDeltaFirst and fcTDeltaLast
#             - x-axis: cycle initial time
#             -    line: separate experiment
#             - subplot: different observed variable
#             -    file: combination of FC lead time, statistic, and bin (if applicable)
# AGGREGATEFC - creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing 
#               aggregated statistics for the period between firstCycleDTime and lastCycleDTime
#             - x-axis: forecast duration
#             -    line: separate experiment
#             - subplot: different observed variable
#             -    file: combination of FC lead time, statistic, and bin
# SNAPSHOTCY2D/AGGREGATEFC2D creates contour maps similar to above with vertical info on y-axis
#             - subplot: combination of experiment and observed variable
#             -    file: combination of FC lead time, statistic, and bin
#             - only applicable to binned observations (e.g., vertical dimension, latitude)
# SNAPSHOTFC  - similar to SNAPSHOTCY, except
#             -    line: different FC lead time
#             - subplot: different observed variable
#             -    file: combination of experiment, statistic, and bin
# SNAPSHOTCY_LAT1 - similar to SNAPSHOTCY, except
#                 -    line: named latitude bins
#                 - subplot: combination of experiment and observed variable
#                 -    file: combination of statistic and forecast length
#                 - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_LAT2 - similar to SNAPSHOTCY, except
#                 -    line: different experiments
#                 - subplot: combination of named latitude bins and observed variable
#                 -    file: combination of statistic and forecast length
#                 - by default, skipped for radiances due to large # figures (slow)
#TODO: AGGREGATEFC_LAT1/2

# NOTE: all *FC* figures require non-zero forecast length

plotTypes = ['SNAPSHOTCY','SNAPSHOTCY2D']
#plotTypes = ['AGGREGATEFC','AGGREGATEFC2D','SNAPSHOTFC']


#Select the diagnostics for plotting
# options: 'omb','oma','obs','bak','ana'
diagNames = ['omb']

#Eventually would like this capability... 
#diagGroups_to_plot = [['omb'],['omb','oma'],['obs','bak','ana']]

#Select the stats for plotting
# options: see pu.allStats
statNames = ['Count','RMS','Mean','STD']
nStat = len(statNames)

#Select the variable for 2D figures
# options: 'P','alt','lat'
# NOTE: radiances can only be binned by 'lat' so far
binVars2D = ['P','alt','lat']

user = os.getenv('USER','jban')

#ASCII statistics file example for cycling run (on cheyenne): 
#statFile = '/glade/scratch/user/pandac/DA_3dvar/2018041500/0/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.txt'
#            |                         |        |          | |                      |             |         |
#                       ^                   ^        ^      ^           ^                  ^           ^
#                  expDirectory          expName  cyDTime fcTDelta statsFilePrefix     DAMethod   ObsSpaceName

example = "SNAP"

if example == "SNAP":
    ## Example settings for SNAP (best when considering single forecast length)
    ## ------------------------------------------------------------------------
    plotTypes = ['SNAPSHOTCY','SNAPSHOTCY2D','SNAPSHOTCY_LAT2']

    user = 'guerrett'
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,4,30,0,0,0)

    #FIRST FORECAST TIME TO INCLUDE
    fcTDeltaFirst = dt.timedelta(days=0)

    #LAST FORECAST TIME TO INCLUDE
    fcTDeltaLast  = dt.timedelta(days=0)

    fcTimeInc  = dt.timedelta(hours=24)
    cyTimeInc  = dt.timedelta(hours=6)

    expLongNames = [ \
                 'VFC0day_gfsana_amsua_no-bias-correction', \
                 'VFC0day_gfsana_amsua_bias-correction' \
                ]

    expNames = [ \
                 'amsua raw obs.', \
                 'amsua bias-corrected' \
                ]

    DAMethods = [ \
                 'verify', \
                 'verify' \
                ]

if example == "AGG":
    # Example settings for AGG (best when considering multiple forecast lengths)
    # --------------------------------------------------------------------------
    plotTypes = ['AGGREGATEFC','AGGREGATEFC2D','SNAPSHOTFC']

    user = 'jban'
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,4,20,0,0,0)

    #FIRST FORECAST TIME TO INCLUDE
    fcTDeltaFirst = dt.timedelta(days=0)

    #LAST FORECAST TIME TO INCLUDE
    fcTDeltaLast  = dt.timedelta(days=10)

    fcTimeInc  = dt.timedelta(hours=24)
    cyTimeInc  = dt.timedelta(hours=24)

    expLongNames = [ \
                 'VFC10day_test20_3denvar_sst_top30km_rho', \
                 'VFC10day_test20_3denvar_sst_top30km_rho' \
                ]

    expNames = [ \
                 'junmei_test', \
                 'junmei_test_copy' \
                ]

    DAMethods = [ \
                 '3dvar_bumpcov', \
                 '3dvar_bumpcov' \
                ]

nExp  = len(expNames)

expDirectory = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac/')

statsFilePrefix = 'diagnostic_stats/stats_'

#plot settings
INTERIOR_AXES_LABELS = True

sciticks = []
for statName in statNames:
    if statName == 'Count': sciticks.append(True)
    else: sciticks.append(False)

max_lines_SNAPSHOTFC = 5

# primary script to be called by main()

def plot_stats_timeseries(nproc, myproc):
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

    if fcTDeltaFirst == fcTDeltaLast and firstCycleDTime == lastCycleDTime:
        print("\n\nERROR: Need a time difference to produce timeseries")
        os._exit(1)

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

    # Retrieve list of ObsSpaceNames for from first experiment (assume identical across expNames)
    FILEPREFIX0 = expDirectory + expLongNames[0] +'/'+cyDTimes_dir[0]+'/'+fcTDeltas_dir[0] +'/' \
                  +statsFilePrefix+DAMethods[0]+"_"
    ObsSpaceNames = []
    for File in glob.glob(FILEPREFIX0+'*.txt'):
       ObsSpaceName = re.sub(".txt","",re.sub(FILEPREFIX0,"",File))
       ObsSpaceInfo_ = ObsSpaceDict.get(ObsSpaceName,[pu.miss_s, 0, pu.nullBinKeys])
       if ObsSpaceInfo_[0] != pu.miss_s and ObsSpaceInfo_[1] > 0:
           ObsSpaceNames.append(ObsSpaceName)
    ObsSpaceNames.sort()

    print("")
    print("Processing data for these ObsSpaceNames:")
    print("---------------------------------------")
    print(ObsSpaceNames)

    # Read stats and make figures for all ObsSpaceNames
    for ObsSpaceName in ObsSpaceNames:
        print("\n==============================")
        print("ObsSpaceName = "+ObsSpaceName)
        print("==============================")

        OSFigPath = "./"+ObsSpaceName+"_figs"
        if not os.path.exists(OSFigPath):
            os.mkdir(OSFigPath)

        # create temporary csv file with statistics from all expNames, fcTDeltas, cyDTimes
        allStats = []
        tmpStatsFile = 'temp_stats_'+ObsSpaceName+'.txt'
        if os.path.exists(tmpStatsFile):
            os.remove(tmpStatsFile)
        fpw = open(tmpStatsFile, 'a')
        for iexp, expName in enumerate(expNames):
            expPrefix = expDirectory + expLongNames[iexp] +'/'
            statsFile = statsFilePrefix+DAMethods[iexp]+'_'+ObsSpaceName+'.txt'
            for ifc, fcTDelta in enumerate(fcTDeltas):
                fcstr=fcTDelta.__str__()
                for icy, cyDTime in enumerate(cyDTimes):
                    #Read all stats from ASCII file for diagName and varName
                    # note: this reading does not include binned data yet
                    cyStatsFile = expPrefix+cyDTimes_dir[icy]+'/'+fcTDeltas_dir[ifc]+'/'+statsFile

                    statLines = \
                        readStatsFile(cyStatsFile)

                    prefix=expName+pu.csvSEP+fcstr+pu.csvSEP+cyDTime.__str__()+pu.csvSEP
                    for il, line in enumerate(statLines):
                        writeLine = re.sub("\n","",line)
#                        writeLine = re.sub("\s+",pu.csvSEP,writeLine)
                        fpw.write(prefix+writeLine+"\n")
        fpw.close()

        # load csv file into pandas DataFrame
        columnNames = ['expName','fcTDelta','cyDTime','ObsSpaceGrp', \
                        'varName','varUnits','diagName','binVar','binVal','binUnits']
        for stat in pu.allFileStats:
            columnNames.append(stat)
        indexNames = ['expName','fcTDelta','cyDTime','ObsSpaceGrp', \
                      'varName','diagName','binVal']
        osDF = pd.read_csv(tmpStatsFile, sep=pu.csvSEP, parse_dates=['cyDTime'], \
                             names=columnNames, index_col=indexNames, \
                             converters={'binVar': str,'binVal': str,'binUnits': str, \
                                         'fcTDelta': pd.to_timedelta,'Count':int} \
                            )
        osDF = osDF.sort_index()

        ##  obsspace
        ObsSpaceGrp = osDF.index.levels[indexNames.index('ObsSpaceGrp')]

        ##  variables
        # get varNames and sort alphabetically
        varNames = osDF.index.levels[indexNames.index('varName')]
        nVars = len(varNames)
        indices = list(range(nVars))
        if ObsSpaceGrp[0] == pu.radiance_s:
            # sort by channel number (int) for radiances
            chlist = []
            for varName in varNames:
                chlist.append(int(''.join(varName.split("_")[-1:])))
            indices.sort(key=chlist.__getitem__)
        else:
            indices.sort(key=varNames.__getitem__)
        varNames = list(map(varNames.__getitem__, indices))

        ##  bins (may want to sort numerical values)
        binVals = osDF.index.levels[indexNames.index('binVal')].tolist()

        # extract units for all varNames from varUnits DF column
        varsLoc = (expNames[0],fcTDeltas[0],cyDTimes[0],ObsSpaceGrp,varNames,diagNames[0],binVals[0])
        varUnitss = osDF.loc[varsLoc,'varUnits'].tolist()

        # define first location for simpler extraction of additional info
        ##  bins
        # extract bin information from binVar and binUnits DF columns
        binsLoc = (expNames[0],fcTDeltas[0],cyDTimes[0],ObsSpaceGrp,varNames[0],diagNames[0],binVals)
        binVars = osDF.loc[binsLoc,'binVar'].tolist()
        binUnitss = osDF.loc[binsLoc,'binUnits'].tolist()

        # convert binVals to numeric type that can be used as axes values
        binNumVals = []
        for binVal in binVals:
            if pu.isfloat(binVal):
                binNumVals.append(float(binVal))
            elif pu.isint(binVal):
                binNumVals.append(int(binVal))
            else:
                binNumVals.append(np.NaN)

        if os.path.exists(tmpStatsFile):
            os.remove(tmpStatsFile)

        #==================================
        # Create figures for all diagNames
        #==================================

        for diagName in diagNames:
            data_labels = []
            for statName in statNames:
                if statName == 'Count': data_labels.append(statName)
                else: data_labels.append(diagName+" "+statName)

            posdef = []
            for statName in statNames:
                #This only applies to the unbounded quantities (omb, oma, ana/bak for velocity)
                if statName == 'Mean' and \
                   diagName == 'omb' or diagName == 'oma': 
                    posdef.append(False)
                else:
                    posdef.append(True)

            # reduce the index space by two dimensions
            #            expName     fcTDelta    cyDTime                    varName              binVal
            diagLoc = (slice(None),slice(None),slice(None),ObsSpaceGrp[0],slice(None),diagName,slice(None))
            diagDF = osDF.xs(diagLoc)



            #=============
            # 1-D figures
            #=============
            if "SNAPSHOTCY" in plotTypes \
               or "AGGREGATEFC" in plotTypes \
               or "SNAPSHOTFC" in plotTypes:

                ## make one figure for unbinned data
                selectBinVals1D = [pu.allBins]

                ## make one figure for each non-float and non-int bin
                #selectBinVals1D = []
                #for binVal in binVals:
                #    if not pu.isfloat(binVal) and not pu.isint(binVal):
                #        selectBinVals1D.append(binVal)

                ## make one figure for each bin
                #selectBinVals1D = []
                #for binVal in binVals:
                #    if binVal != pu.allBins:
                #        selectBinVals1D.append(binVal)

                ## make one figure for unbinned data and one for each bin
                #selectBinVals1D = binVals

                # top-level 1-D plot settings
                nsubplots = nVars
                nx = np.int(np.ceil(np.sqrt(nsubplots)))
                ny = np.int(np.ceil(np.true_divide(nsubplots,nx)))

            if "SNAPSHOTCY" in plotTypes \
                and len(selectBinVals1D) > 0:
                print("\nGenerating SNAPSHOTCY figures...")
                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal == pu.allBins:
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
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].min()
                                ymax = diagDF.loc[varLoc,statName].max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    expLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    linesVals.append(diagDF.loc[expLoc,statName].to_numpy())

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, ny, nx, nsubplots, ivar, \
                                                cyDTimes, linesVals, expNames, \
                                                title, data_labels[istat], \
                                                sciticks[istat], posdef[istat], \
                                                ymin = ymin, ymax = ymax )

                            filename = OSFigPath+'/SNAPSHOTCY_TSeries_%sday_%s_%s_%s'%( \
                                       fcTDeltas_dir[ifc],ObsSpaceName, \
                                       diagName,statName)+filebin

                            pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                # end binName loop

            # end SNAPSHOTCY

            if "AGGREGATEFC" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and len(fcTDeltas) > 1:
                print("\nGenerating AGGREGATEFC figures...")
                subplot_size = 1.9
                aspect = 0.6

                # apply aggregation over cyDTime via including all other free indices in groupby
                aggDF = diagDF.groupby(
                        ['expName','fcTDelta','varName','binVal']).apply(pu.aggStatsSeries)

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal == pu.allBins:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #file loop 2
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                        #subplot loop
                        for ivar, varName in enumerate(varNames):
                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,selectBinVals1D)
                            ymin = aggDF.loc[varLoc,statName].min()
                            ymax = aggDF.loc[varLoc,statName].max()

                            #collect aggregated statNames, varying across fcTDelta
                            linesVals = []
                            for expName in expNames:
                                #                                    fcTDelta
                                linesVals.append(aggDF.loc[(expName,slice(None),varName,binVal), \
                                                            statName].to_numpy())

                            # define subplot title
                            title = varName
                            if varUnitss[ivar] != pu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"
                            title = title+titlebin

                            # perform subplot agnostic plotting (all expNames)
                            plotTimeSeries( fig, ny, nx, nsubplots, ivar, \
                                            fcTDeltas, linesVals, expNames, \
                                            title, data_labels[istat], \
                                            sciticks[istat], posdef[istat], \
                                            ymin = ymin, ymax = ymax )

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFC_TSeries_%s-%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   diagName,statName)+filebin

                        pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                    # end statName loop

                # end binName loop

            # end AGGREGATEFC

            if "SNAPSHOTFC" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and len(fcTDeltas) > 1:
                print("\nGenerating SNAPSHOTFC figures...")
                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal == pu.allBins:
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
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].min()
                                ymax = diagDF.loc[varLoc,statName].max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                fcTDeltas_labels = []
                                for ifc, fcTDelta in enumerate(fcTDeltas):

                                    #Setting to avoid over-crowding
                                    if ifc > (max_lines_SNAPSHOTFC-1): continue

                                    fcTDelta_sec = pu.TDeltas2Seconds([fcTDelta])
                                    fcTDeltas_labels.append(pu.timeTicks(fcTDelta_sec[0],0))
                                    #                           cyDTime
                                    fcLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    linesVals.append(diagDF.loc[fcLoc,statName].to_numpy())

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, ny, nx, nsubplots, ivar, \
                                                cyDTimes, linesVals, fcTDeltas_labels, \
                                                title, data_labels[istat], \
                                                sciticks[istat], posdef[istat], \
                                                ymin = ymin, ymax = ymax )

                            expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                            filename = OSFigPath+'/SNAPSHOTFC_TSeries_%s_%s_%s_%s'%( \
                                       expName_file,ObsSpaceName, \
                                       diagName,statName)+filebin

                            pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                # end binName loop

            # end SNAPSHOTFC

            if "SNAPSHOTCY_LAT1" in plotTypes \
                and len(pu.namedLatBandsStrVals) > 0 \
                and ObsSpaceGrp[0] != pu.radiance_s:

                print("\nGenerating SNAPSHOTCY_LAT1 figures...")
                selectBinVals1D = pu.namedLatBandsStrVals
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
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                        iplot = 0

                        #subplot loop 1
                        for expName in expNames:

                            #subplot loop 2
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].min()
                                ymax = diagDF.loc[varLoc,statName].max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for binVal in selectBinVals1D:
                                    #                           cyDTime
                                    expLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    linesVals.append(diagDF.loc[expLoc,statName].to_numpy())

                                # define subplot title
                                title = expName+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, ny, nx, nsubplots, iplot, \
                                                cyDTimes, linesVals, pu.namedLatBandsStrVals, \
                                                title, data_labels[istat], \
                                                sciticks[istat], posdef[istat], \
                                                ymin = ymin, ymax = ymax )
                                iplot = iplot + 1

                        filename = OSFigPath+'/SNAPSHOTCY_LAT1_TSeries_%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[ifc],ObsSpaceName, \
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                # end binName loop

            # end SNAPSHOTCY_LAT1

            if "SNAPSHOTCY_LAT2" in plotTypes \
                and len(pu.namedLatBandsStrVals) > 0 \
                and ObsSpaceGrp[0] != pu.radiance_s:
                print("\nGenerating SNAPSHOTCY_LAT2 figures...")
                selectBinVals1D = pu.namedLatBandsStrVals
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
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                        iplot = 0
                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):


                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                            ymin = diagDF.loc[varLoc,statName].min()
                            ymax = diagDF.loc[varLoc,statName].max()

                            #subplot loop 2
                            for binVal in selectBinVals1D:

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    expLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    linesVals.append(diagDF.loc[expLoc,statName].to_numpy())

                                # define subplot title
                                title = binVal+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, ny, nx, nsubplots, iplot, \
                                                cyDTimes, linesVals, expNames, \
                                                title, data_labels[istat], \
                                                sciticks[istat], posdef[istat], \
                                                ymin = ymin, ymax = ymax )

                                iplot = iplot + 1

                        #expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                        filename = OSFigPath+'/SNAPSHOTCY_LAT2_TSeries_%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[ifc],ObsSpaceName, \
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                # end binName loop

            # end SNAPSHOTCY_LAT2



            #=============
            # 2-D figures
            #=============
            for binVar2D in binVars2D:
                if "SNAPSHOTCY2D" in plotTypes \
                   or "AGGREGATEFC2D" in plotTypes:
                    selectBinVals2D = []
                    for ibin, binVal in enumerate(binVals):
                        #Only numerical bins
                        if (pu.isfloat(binVal) or pu.isint(binVal)) \
                           and binVars[ibin] == binVar2D:
                            selectBinVals2D.append(binVal)

                    # top-level 2-D plot settings
                    nsubplots = nExp * nVars
                    nx = nExp
                    ny = nVars

                    # bin info
                    selectBinNumVals = []
                    for binVal in selectBinVals2D:
                        ibin = binVals.index(binVal)
                        selectBinNumVals.append(binNumVals[ibin])

                        # assume all bins represent same variable/units
                        if binVars[ibin] == pu.miss_s:
                            binVar = ''
                        else:
                            binVar = binVars[ibin]
                        ylabel = binVar

                        if binUnitss[ibin] != pu.miss_s:
                            ylabel = ylabel+" ("+binUnitss[ibin]+")"

                        # invert y-axis for pressure bins
                        pressure_dict = pu.varDict.get('air_pressure',['',''])
                        invert_y_axis = (pressure_dict[1] == binVar)

                    # sort bins by numeric value
                    indices = list(range(len(selectBinNumVals)))
                    indices.sort(key=selectBinNumVals.__getitem__)
                    selectBinNumVals = list(map(selectBinNumVals.__getitem__, indices))
                    selectBinVals2D = list(map(selectBinVals2D.__getitem__, indices))

                if "SNAPSHOTCY2D" in plotTypes \
                    and len(selectBinVals2D) > 0:
                    print("Generating SNAPSHOTCY2D figures across "+binVar2D+"...")
                    subplot_size = 2.4
                    aspect = 0.65

                    #file loop 1
                    for ifc, fcTDelta in enumerate(fcTDeltas):
                        #print("fcTDelta = "+fcTDeltas_dir[ifc])

                        #file loop 2
                        for istat, statName in enumerate(statNames):
                            # establish a new figure
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                            iplot = 0

                            #subplot loop 1
                            for ivar, varName in enumerate(varNames):
                                # use specific c-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals2D)
                                cmin = diagDF.loc[varLoc,statName].min()
                                cmax = diagDF.loc[varLoc,statName].max()

                                #subplot loop 2
                                # letting cyDTime and binVal vary
                                for expName in expNames:
                                    contourVals = np.empty((len(selectBinVals2D),len(cyDTimes)))
                                    for ibin, binVal in enumerate(selectBinVals2D):
                                        #                           cyDTime
                                        binLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                        contourVals[ibin,:] = diagDF.loc[binLoc,statName].to_numpy()

                                    # define subplot title
                                    title = expName+"\n"+varName
                                    if varUnitss[ivar] != pu.miss_s:
                                        title = title+" ("+varUnitss[ivar]+")"

                                    # perform subplot agnostic plotting (all expNames)
                                    plotTimeSeries2D( fig, ny, nx, nsubplots, iplot, \
                                                      cyDTimes, selectBinNumVals, contourVals, \
                                                      title, data_labels[istat], \
                                                      sciticks[istat], posdef[istat], \
                                                      ylabel, invert_y_axis, \
                                                      cmin = cmin, cmax = cmax )
                                    iplot = iplot + 1

                            filename = OSFigPath+'/SNAPSHOTCY2D_%s_TSeries_%sday_%s_%s_%s'%( \
                                       binVar2D,fcTDeltas_dir[ifc],ObsSpaceName, \
                                        diagName,statName)

                            pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                        # end statName loop

                    # end fcTDelta loop

                # end SNAPSHOTCY2D

                if "AGGREGATEFC2D" in plotTypes \
                    and len(selectBinVals2D) > 0 \
                    and len(fcTDeltas) > 1:
                    print("Generating AGGREGATEFC2D figures across "+binVar2D+"...")
                    subplot_size = 2.4
                    aspect = 0.55

                    # apply aggregation over cyDTime via including all other free indices in groupby
                    aggDF = diagDF.groupby(
                            ['expName','fcTDelta','varName','binVal']).apply(pu.aggStatsSeries)

                    #file loop 1
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, INTERIOR_AXES_LABELS)

                        iplot = 0

                        #subplot loop
                        for ivar, varName in enumerate(varNames):
                            # use specific c-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,selectBinVals2D)
                            cmin = aggDF.loc[varLoc,statName].min()
                            cmax = aggDF.loc[varLoc,statName].max()

                            #collect aggregated statNames, varying across fcTDelta
                            for expName in expNames:

                                contourVals = np.empty((len(selectBinVals2D),len(fcTDeltas)))
                                for ibin, binVal in enumerate(selectBinVals2D):
                                    #                                         fcTDelta
                                    contourVals[ibin,:] = aggDF.loc[(expName,slice(None),varName,binVal), \
                                                                    statName].to_numpy()

                                # define subplot title
                                title = expName+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries2D( fig, ny, nx, nsubplots, iplot, \
                                                  fcTDeltas, selectBinNumVals, contourVals, \
                                                  title, data_labels[istat], \
                                                  sciticks[istat], posdef[istat], \
                                                  ylabel, invert_y_axis, \
                                                  cmin = cmin, cmax = cmax )
                                iplot = iplot + 1

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFC2D_%s_TSeries_%s-%sday_%s_%s_%s'%( \
                                   binVar2D,fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, 'png', INTERIOR_AXES_LABELS)

                    # end statName loop

                # end AGGREGATEFC2D

            # end binVars2D loop

        # end diagName loop

    # end ObsSpaceName loop


def readStatsFile(statsFile):
# INPUTS
# statsFile (string) - file containing statistical information
#
# OUTPUTS
# statLines (string) - all lines in file (ascii)

    if not os.path.exists(statsFile):
        print ('\n\nERROR: statsFile is missing: '+statsFile)
        os._exit(1)

    fp = open(statsFile, 'r')
    statLines = []
    for line in fp:
        statLines.append(line)
    fp.close()
    return statLines


def plotTimeSeries(fig, ny, nx, nplots, iplot, \
                   xDates, linesVals, allLineLabels, \
                   title, ylabel, \
                   sciticks, posdef, \
                   *args, **kwargs):

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #float xVals
    xVals = pu.TDeltas2Seconds(xDates)

    #add lines
    lineLabels = []
    plotVals = []
    for iline, lineVals in enumerate(linesVals):
        if not all(np.isnan(lineVals)):
            # Plot line for each lineVals that has non-missing data
            ax.plot(xVals, lineVals, pu.plotMarkers[iline], markersize=1.5, linewidth=0.5)
            lineLabels.append(allLineLabels[iline])
            plotVals.append(lineVals)


    if len(lineLabels) == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #legend
    ax.legend(lineLabels, loc='upper right',fontsize=3,frameon=False)

    #axes settings
    pu.format_x_for_dates(ax, xDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    # add horizontal zero line for unbounded quantities
    if not posdef:
        ax.plot([xVals[0], xVals[-1]], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize y-limits
    ymin_ = kwargs.get('ymin', np.NaN)
    ymax_ = kwargs.get('ymax', np.NaN)
    if not np.isnan(ymin_) and not np.isnan(ymax_):
        ymin = ymin_
        ymax = ymax_
    else:
        ymin = np.nanmin(plotVals)
        ymax = np.nanmax(plotVals)

    ymaxabs=max(abs(ymin),abs(ymax))
    if ymaxabs == 0.0 or np.isnan(ymaxabs):
        minyval = 0.0
        maxyval = 1.0
    else:
        roundfact = np.round(1. / 10.0 ** np.floor(np.log10(ymaxabs)))
        if np.isnan(roundfact) or roundfact <= 0.0: roundfact = 1.0
        maxyval = np.ceil(  ymaxabs*roundfact ) / roundfact

        if posdef:
            minyval = np.floor(   ymin*roundfact ) / roundfact
        else:
            minyval = np.floor( - ymaxabs*roundfact ) / roundfact

    if not np.isnan(minyval) and not np.isnan(maxyval):
        ax.set_ylim(minyval,maxyval)

    ax.grid()

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not INTERIOR_AXES_LABELS \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if INTERIOR_AXES_LABELS or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)

    return

def plotTimeSeries2D(fig, ny, nx, nplots, iplot, \
                     xDates, yVals, contourVals, \
                     title, clabel, \
                     sciticks, posdef, \
                     ylabel, invert_y_axis, \
                     *args, **kwargs):

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    xVals = pu.TDeltas2Seconds(xDates)

    # standardize c-limits
    cmin_ = kwargs.get('cmin', np.NaN)
    cmax_ = kwargs.get('cmax', np.NaN)
    if not np.isnan(cmin_) and not np.isnan(cmax_):
        cmin = cmin_
        cmax = cmax_
    else:
        cmin = np.nanmin(contourVals)
        cmax = np.nanmax(contourVals)

    cmaxabs=max(abs(cmin),abs(cmax))
    if cmaxabs == 0.0:
        mincval =  0.0
        maxcval =  1.0
    else:
        roundfact = np.round(1. / 10.0 ** np.floor(np.log10(cmaxabs)))
        if np.isnan(roundfact) or roundfact <= 0.0: roundfact = 1.0
        maxcval = np.ceil(  cmaxabs*roundfact ) / roundfact

        if posdef:
            mincval = np.floor(   cmin*roundfact ) / roundfact
        else:
            mincval = np.floor( - cmaxabs*roundfact ) / roundfact

    if posdef:
        cmapName = 'BuPu'
        nlevs = 18
    else:
        cmapName = 'seismic'
        nlevs = 28

    # plot contour
    # option 1: smoothed contours
    #cp = ax.contourf(xVals, yVals, contourVals, nlevs, cmap=cmapName, extend='both', \
    #     vmin=mincval, vmax=maxcval)

    # option 2: pixel contours
    levels = MaxNLocator(nbins=nlevs).tick_values(mincval,maxcval)
    cmap = plt.get_cmap(cmapName)
    cmap.set_bad(color = 'k', alpha = 1.0)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    xVals_pcolor, yVals_pcolor = transformXY_for_pcolor(xVals,yVals)
    cp = ax.pcolormesh(xVals_pcolor, yVals_pcolor, contourVals, cmap=cmap, norm = norm)

    #title
    ax.set_title(title,fontsize=5)

    #axes settings
    pu.format_x_for_dates(ax, xDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not INTERIOR_AXES_LABELS \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if INTERIOR_AXES_LABELS or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)
    if INTERIOR_AXES_LABELS or ix == nx-1:
        #colorbar
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(contourVals)
        if not np.isnan(mincval) and not np.isnan(maxcval):
            m.set_clim(mincval,maxcval)
        cb = plt.colorbar(m, ax=ax)
        #scientific formatting
        if sciticks:
            cb.ax.set_yticklabels(['{:.2e}'.format(float(x.get_text())) \
                               for x in cb.ax.get_yticklabels()])
        cb.ax.tick_params(labelsize=3)
        cb.set_label(clabel,fontsize=5)

    if invert_y_axis:
        ax.invert_yaxis()

    # optionally add a grid
    #ax.grid()

    return

def transformXY_for_pcolor(xs,ys):
    # adjust centered x and y values to edges to work with pcolormesh 
    # note: works best for regularly spaced data
    xs_diff = xs[1] - xs[0]
    # extend xs by 2
    # fill in first endpoint
    xs_extend = [xs[0]-xs_diff]
    # fill in internal values
    for x in xs: xs_extend.append(x)
    # fill in last endpoint
    xs_extend.append(xs_extend[-1]+(xs[-1]-xs[-2]))
    # calculate the midpoints
    xs_pcolormesh_midpoints = []
    for ii, x in enumerate(xs_extend[:-1]):
        xs_pcolormesh_midpoints.append(x+0.5*(xs_extend[ii+1] - xs_extend[ii]))

    ys_diff = ys[1] - ys[0]
    # extend ys by 2
    # fill in first endpoint
    ys_extend = [ys[0]-ys_diff]
    # fill in internal values
    for y in ys: ys_extend.append(y)
    # fill in last endpoint
    ys_extend.append(ys_extend[-1]+(ys[-1]-ys[-2]))
    # calculate the midpoints
    ys_pcolormesh_midpoints = []
    for ii, y in enumerate(ys_extend[:-1]):
        ys_pcolormesh_midpoints.append(y+0.5*(ys_extend[ii+1] - ys_extend[ii]))

    return xs_pcolormesh_midpoints, ys_pcolormesh_midpoints

def main():
    nproc, myproc = pu.par_args(sys.argv[:])
    plot_stats_timeseries(nproc, myproc)

if __name__ == '__main__': main()
