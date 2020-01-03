from copy import deepcopy
import datetime as dt
import glob
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
#import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import plot_utils as pu
import re
import os
import sys

# This script can be executed normally OR with optional arguments -n and -i
# in order to run with GNU parallel. See plot_utils.par_args for more information.

#Select the type of plots
# options: 'SNAPSHOTCY', 'AGGREGATEFC','SNAPSHOTCY2D', 'AGGREGATEFC2D', 'SNAPSHOTFC'
# SNAPSHOTCY - creates a timeseries figure between firstCycleDTime and lastCycleDTime
#               for each forecast length between fcTDeltaFirst and fcTDeltaLast
#             - x-axis: cycle initial time
#             -    line: per experiment
#             - subplot: per obs. variable
#             -    file: combination of FC lead time, statistic, and bin (if applicable)
# AGGREGATEFC - creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing 
#               aggregated statistics for the period between firstCycleDTime and lastCycleDTime
#             - x-axis: forecast duration
#             -    line: per experiment
#             - subplot: per obs. variable
#             -    file: combination of statistic and bin
# SNAPSHOTCY2D/AGGREGATEFC2D creates contour maps similar to above with vertical info on y-axis
#             - only applicable to binned observations (e.g., vertical dimension, latitude)
#             - subplot: column by experiment, row by obs. variable
#             SNAPSHOTCY2D
#             -    file: combination of FC lead time and statistic
#             AGGREGATEFC2D
#             -    file: per statistic
# AGGREGATEFCPROFILE similar to AGGREGATEFC2D, except
#             - each vertical column of data is plotted as a profile on
#               a separate set of axes
#             -    line: per experiment
#             - subplot: column by lead time, row by obs. variable
#             -    file: per statistic
#             - MAX_FC_SUBFIGS determines number of FC lead times to include
MAX_FC_SUBFIGS = 3
# SNAPSHOTFC  similar to SNAPSHOTCY, except
#             -    line: per FC lead time
#             - subplot: per obs. variable
#             -    file: combination of experiment, statistic, and bin
#             - MAX_FC_LINES determines number of FC lead times to include
MAX_FC_LINES = 5
# SNAPSHOTCY_LAT1 similar to SNAPSHOTCY, except
#             -    line: named latitude bins
#             - subplot: column by experiment, row by obs. variable
#             -    file: combination of statistic and forecast length
#             - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_LAT2 similar to SNAPSHOTCY, except
#             -    line: per experiment
#             - subplot: column by LAT bin, row by obs. variable
#             -    file: combination of statistic and forecast length
#             - by default, skipped for radiances due to large # figures (slow)
# SNAPSHOTCY_QC similar to SNAPSHOTCY, except
#             - only plots 'Count' statistic
#             -    line: named QC bins
#             - subplot: column by experiment, row by obs. variable
#             -    file: forecast length
# AGGREGATEFC_DIFFCI/AGGREGATEFCPROFILE_DIFFCI similar to their per experiment counterparts
#             - shows difference between experiment(s) and control
#             - control is selected using cntrlExpInd
cntrlExpInd = 0
#             - statistics are selected with bootStrapStats
bootStrapStats = []
for x in pu.sampleableAggStats:
    if x != 'Count': bootStrapStats.append(x)


#TODO: AGGREGATEFC_LAT1/2 named latitude figures w/ FC lead time x-axis

# NOTE: all *FC* figures require non-zero forecast length


#Select the stats for plotting
# options: see pu.allFileStats
statNames = ['Count','Mean','RMS','STD']


#Select the variable for 2D figures
# options: 'P','alt','lat'
# NOTE: radiances can only be binned by 'lat' so far
binVars2D = ['P','alt','lat']

user = os.getenv('USER','jban')

plotGroup = "multiFCLen" # "singleFCLen"


# multiFCLen ASCII statistics file example for cycling run (on cheyenne):
#statFile = '/glade/scratch/user/pandac/DA_3dvar/2018041500/0/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.txt'
#            |                         |        |          | |                      |             |         |
#                       ^                   ^        ^      ^           ^                  ^           ^
#                  expDirectory          expName  cyDTime fcTDelta statsFilePrefix     DAMethod   ObsSpaceName

# singleFCLen ASCII statistics file example for cycling run (on cheyenne):
#statFile = '/glade/scratch/user/pandac/DA_3dvar/2018041500/diagnostic_stats/stats_3dvar_bumpcov_amsua_n19.txt'
#            |                         |        |          |                      |             |         |
#                       ^                   ^        ^                ^                  ^           ^
#                  expDirectory          expName  cyDTime      statsFilePrefix     DAMethod   ObsSpaceName


plotTypes = []

if plotGroup == "singleFCLen":
    #Select the diagnostics for plotting
    # options: 'omb','oma','obs','bak','ana'
    diagNames = ['omb','oma']
    #TODO (maybe): multiple diagNames on same subplot
    #diagGroups_to_plot = [['omb'],['omb','oma'],['obs','bak','ana']]

    ## plotTypes for considering single forecast length
    ## ------------------------------------------------------------------------
    plotTypes.append('SNAPSHOTCY')
    plotTypes.append('SNAPSHOTCY2D')
    plotTypes.append('SNAPSHOTCY_LAT2')
    plotTypes.append('SNAPSHOTCY_QC')

    user = 'jban'

    #First and Last CYCLE dates
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,5,14,0,0,0)
    cyTimeInc  = dt.timedelta(hours=6)

    #First and Last FORECAST durations 
    fcTDeltaFirst = dt.timedelta(days=0)
    fcTDeltaLast  = dt.timedelta(days=0)
    fcTimeInc  = dt.timedelta(hours=24)

    expLongNames = [ \
                 'conv60it-2node/DADIAG', \
                 'amua60it-2node/DADIAG' \
                ]

    expNames = [ \
                 'conv', \
                 'conv+amsua' \
                ]

    DAMethods = [ \
                 '3denvar_bumploc', \
                 '3denvar_bumploc' \
                ]

if plotGroup == "multiFCLen":
    #Select the diagnostics for plotting
    # options: 'omb','obs','bak'
    diagNames = ['omb']

    # plotTypes for considering multiple forecast lengths
    # --------------------------------------------------------------------------
    plotTypes.append('AGGREGATEFC')
    plotTypes.append('AGGREGATEFC2D')
    plotTypes.append('SNAPSHOTFC')
    plotTypes.append('AGGREGATEFC_DIFFCI')
    plotTypes.append('AGGREGATEFCPROFILE')
    plotTypes.append('AGGREGATEFCPROFILE_DIFFCI')

    user = 'jban'

    #First and Last CYCLE dates
    firstCycleDTime = dt.datetime(2018,4,15,0,0,0)
    lastCycleDTime = dt.datetime(2018,5,4,0,0,0)
    cyTimeInc  = dt.timedelta(hours=24)

    #First and Last FORECAST durations 
    fcTDeltaFirst = dt.timedelta(days=0)
    fcTDeltaLast  = dt.timedelta(days=10)
    fcTimeInc  = dt.timedelta(hours=24)
    #TODO: define FC directory names with in terms of seconds or d_hh-mm-ss
    #      in order to allow for increments < 1 day --> modify workflow scripts

    expLongNames = [ \
                 'conv60it-2node/OMF', \
                 'amua60it-2node/OMF' \
                ]

    expNames = [ \
                 'conv', \
                 'conv+amsua' \
                ]

    DAMethods = [ \
                 '3dvar', \
                 '3dvar' \
                ]

nExp  = len(expNames)


cntrlName = expNames[cntrlExpInd]
noncntrlExpNames = [x for x in expNames if x != cntrlName]
print('\nControl Experiment: '+cntrlName)
print('\nNon-control Experiment(s): ')
print(noncntrlExpNames)

expDirectory = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac/')

statsFilePrefix = 'diagnostic_stats/stats_'

#plot settings
FULL_SUBPLOT_LABELS = True

sciticks = []
for statName in statNames:
    if statName == 'Count': sciticks.append(True)
    else: sciticks.append(False)

#
# primary script to be called by main()
#

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
        print("\n\nERROR: Time difference is required for all plotTypes (either forecast or multiple cycles)")
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

    # Retrieve list of ObsSpaceNames from all experiments
    expsObsSpaceNames = []
    for expName in expLongNames:
        dateDir = cyDTimes_dir[0]
        if plotGroup == "multiFCLen":
            dateDir = dateDir+'/'+fcTDeltas_dir[0]

        FILEPREFIX0 = expDirectory + expName +'/'+dateDir+'/' \
                      +statsFilePrefix+DAMethods[0]+"_"

        ObsSpaceNames = []
        for File in glob.glob(FILEPREFIX0+'*.txt'):
           ObsSpaceName = re.sub(".txt","",re.sub(FILEPREFIX0,"",File))
           ObsSpaceInfo_ = ObsSpaceDict.get(ObsSpaceName,pu.nullObsSpaceInfo)
           if ObsSpaceInfo_['process']:
               ObsSpaceNames.append(ObsSpaceName)
        expsObsSpaceNames.append(ObsSpaceNames)

    # Remove ObsSpaceNames that are not common to all experiments
    ObsSpaceNames = deepcopy(expsObsSpaceNames[0])
    if len(expsObsSpaceNames) > 1:
        for expObsSpaceNames in expsObsSpaceNames[1:]:
            for ObsSpaceName in expObsSpaceNames:
                if not (ObsSpaceName in expObsSpaceNames):
                    ObsSpaceNames.remove(ObsSpaceName)
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

        print("\nConcatenating ASCII files...")

        # create temporary csv file with statistics from all expNames, fcTDeltas, cyDTimes
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
                    dateDir = cyDTimes_dir[icy]
                    if plotGroup == "multiFCLen":
                        dateDir = dateDir+'/'+fcTDeltas_dir[ifc]
                    cyStatsFile = expPrefix+dateDir+'/'+statsFile

                    statLines = \
                        readStatsFile(cyStatsFile)

                    prefix=expName+pu.csvSEP+fcstr+pu.csvSEP+cyDTime.__str__()+pu.csvSEP
                    for il, line in enumerate(statLines):
                        writeLine = re.sub("\n","",line)
#                        writeLine = re.sub("\s+",pu.csvSEP,writeLine)
                        fpw.write(prefix+writeLine+"\n")
        fpw.close()

        print("\nPopulating pandas DataFrame...")

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

        OSFigPath = "./"+ObsSpaceName+"_figs"
        if not os.path.exists(OSFigPath):
            os.mkdir(OSFigPath)

        for diagName in diagNames:
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
                if statName == 'Mean' and \
                   (diagName == 'omb' or diagName == 'oma'): 
                    signdef.append(False)
                else:
                    signdef.append(True)

            # reduce the index space by two dimensions
            #            expName     fcTDelta    cyDTime                    varName              binVal
            diagLoc = (slice(None),slice(None),slice(None),ObsSpaceGrp[0],slice(None),diagName,slice(None))
            diagDF = osDF.xs(diagLoc)



            #=============
            # 1-D figures
            #=============
            if "SNAPSHOTCY" in plotTypes \
               or "AGGREGATEFC" in plotTypes \
               or "AGGREGATEFC_DIFFCI" in plotTypes \
               or "SNAPSHOTFC" in plotTypes:

                ## make one figure for unbinned data
                selectBinVals1D = pu.goodBinNames

                ## make one figure for each non-float and non-int bin
                #selectBinVals1D = []
                #for binVal in binVals:
                #    if not pu.isfloat(binVal) and not pu.isint(binVal):
                #        selectBinVals1D.append(binVal)

                ## make one figure for each bin
                #selectBinVals1D = []
                #for binVal in binVals:
                #    if not (binVal in pu.goodBinNames):
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
                    elif binVal in pu.goodBinNames:
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
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].dropna().min()
                                ymax = diagDF.loc[varLoc,statName].dropna().max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                               ['expName','fcTDelta','varName','binVal'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, \
                                                cyDTimes, linesVals, expNames, \
                                                title, data_labels[istat], \
                                                sciticks[istat], signdef[istat], \
                                                ny, nx, nsubplots, ivar, \
                                                ymin = ymin, ymax = ymax )

                            filename = OSFigPath+'/SNAPSHOTCY_TSeries_%sday_%s_%s_%s'%( \
                                       fcTDeltas_dir[ifc],ObsSpaceName, \
                                       diagName,statName)+filebin

                            pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                # end binName loop

            # end SNAPSHOTCY


            if "AGGREGATEFC" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and nFC > 1:
                print("\nGenerating AGGREGATEFC figures...")
                subplot_size = 1.9
                aspect = 0.6

                # apply aggregation over cyDTime via including all other free indices in groupby
                aggDF = diagDF.groupby(
                        ['expName','fcTDelta','varName','binVal']).apply(pu.aggStatsDF)

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal in pu.goodBinNames:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #file loop 2
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        #subplot loop
                        for ivar, varName in enumerate(varNames):
                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,selectBinVals1D)
                            ymin = aggDF.loc[varLoc,statName].dropna().min()
                            ymax = aggDF.loc[varLoc,statName].dropna().max()

                            #collect aggregated statNames, varying across fcTDelta
                            linesVals = []
                            for expName in expNames:
                                #                                   fcTDelta
                                linesVals.append(aggDF.loc[(expName,slice(None),varName,binVal), \
                                                            statName].to_numpy())

                            # define subplot title
                            title = varName
                            if varUnitss[ivar] != pu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"
                            title = title+titlebin

                            # perform subplot agnostic plotting (all expNames)
                            plotTimeSeries( fig, \
                                            fcTDeltas, linesVals, expNames, \
                                            title, fcdata_labels[istat], \
                                            sciticks[istat], signdef[istat], \
                                            ny, nx, nsubplots, ivar, \
                                            ymin = ymin, ymax = ymax )

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFC_TSeries_%s-%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   fcDiagName,statName)+filebin

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                    # end statName loop

                # end binName loop

            # end AGGREGATEFC


            if "AGGREGATEFC_DIFFCI" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and nFC > 1:
                print("\nGenerating AGGREGATEFC_DIFFCI figures...")
                subplot_size = 1.9
                aspect = 0.6

                # nx = 1
                # ny = nVars
                # nsubplots = nx * ny

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal in pu.goodBinNames:
                        titlebin = ""
                        filebin = ""
                    else:
                        titlebin = " @ "+binVal
                        filebin = "_"+binVal

                    #figure loop 2
                    for statName in bootStrapStats:

                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            # define subplot title
                            title = varName
                            if varUnitss[ivar] != pu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"
                            title = title+titlebin

                            linesVals = {}
                            for trait in pu.ciTraits: linesVals[trait] = []

                            for expName in noncntrlExpNames:

                                lineVals = {}
                                for trait in pu.ciTraits: lineVals[trait] = []

                                for fcTDelta in fcTDeltas:
                                    #                              cyDTime
                                    cntrlLoc = (cntrlName,fcTDelta,slice(None),varName,binVal)

                                    #                          cyDTime
                                    expLoc = (expName,fcTDelta,slice(None),varName,binVal)

                                    ciVals = pu.bootStrapClusterFunc( \
                                                 X = diagDF.loc[expLoc,:], \
                                                 Y = diagDF.loc[cntrlLoc,:], \
                                                 n_samples = 10000, \
                                                 statNames = [statName])

                                    for trait in pu.ciTraits:
                                        lineVals[trait].append(ciVals[statName][trait][0])

                                for trait in pu.ciTraits:
                                    linesVals[trait].append( \
                                        lineVals[trait] )

                            # use specific y-axis limits for each varName
                            ymin = np.nanmin(linesVals[pu.cimin])
                            ymax = np.nanmax(linesVals[pu.cimax])

                            # perform subplot agnostic plotting (all expNames)
                            plotTimeSeries( fig, \
                                            fcTDeltas, linesVals[pu.cimean], \
                                            noncntrlExpNames, \
                                            title, \
                                            statName+"("+fcDiagName+"): [EXP - CTRL]", \
                                            False, False, \
                                            ny, nx, nsubplots, ivar, \
                                            linesValsMinCI = linesVals[pu.cimin], \
                                            linesValsMaxCI = linesVals[pu.cimax], \
                                            ymin = ymin, ymax = ymax, \
                                            lineAttribOffset = 1)

                        # end varName loop

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFC_DIFFCI_TSeries_%s-%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   fcDiagName,statName)+filebin

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)
                    # end statName loop

                # end binName loop

            # end AGGREGATEFC

            if "SNAPSHOTFC" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and nFC > 1:
                print("\nGenerating SNAPSHOTFC figures...")
                subplot_size = 1.9
                aspect = 0.75

                #file loop 1
                for binVal in selectBinVals1D:
                    if pu.isfloat(binVal) or pu.isint(binVal):
                        ibin = binVals.index(binVal)
                        titlebin = " @ "+binVars[ibin]+"="+binVal+" "+binUnitss[ibin]
                        filebin = "_"+binVal+binUnitss[ibin]
                    elif binVal in pu.goodBinNames:
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
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                            #subplot loop
                            for ivar, varName in enumerate(varNames):
                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].dropna().min()
                                ymax = diagDF.loc[varLoc,statName].dropna().max()

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
                                    dataLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                               ['expName','fcTDelta','varName','binVal'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+titlebin

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, \
                                                xsVals, linesVals, fcTDeltas_labels, \
                                                title, data_labels[istat], \
                                                sciticks[istat], signdef[istat], \
                                                ny, nx, nsubplots, ivar, \
                                                ymin = ymin, ymax = ymax )

                            expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                            filename = OSFigPath+'/SNAPSHOTFC_TSeries_%s_%s_%s_%s'%( \
                                       expName_file,ObsSpaceName, \
                                       fcDiagName,statName)+filebin

                            pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                # end binName loop

            # end SNAPSHOTFC


            selectBinVals1D = pu.namedLatBandsStrVals
            if "SNAPSHOTCY_LAT1" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and ObsSpaceGrp[0] != pu.radiance_s:

                print("\nGenerating SNAPSHOTCY_LAT1 figures...")
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
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        iplot = 0

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            #subplot loop 2
                            for expName in expNames:

                                # use specific y-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                                ymin = diagDF.loc[varLoc,statName].dropna().min()
                                ymax = diagDF.loc[varLoc,statName].dropna().max()

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for binVal in selectBinVals1D:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                               ['expName','fcTDelta','varName','binVal'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = expName+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, \
                                                cyDTimes, linesVals, selectBinVals1D, \
                                                title, data_labels[istat], \
                                                sciticks[istat], signdef[istat], \
                                                ny, nx, nsubplots, iplot, \
                                                ymin = ymin, ymax = ymax )
                                iplot = iplot + 1

                        filename = OSFigPath+'/SNAPSHOTCY_LAT1_TSeries_%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[ifc],ObsSpaceName, \
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                # end binName loop

            # end SNAPSHOTCY_LAT1


            selectBinVals1D = pu.namedLatBandsStrVals
            if "SNAPSHOTCY_LAT2" in plotTypes \
                and len(selectBinVals1D) > 0 \
                and ObsSpaceGrp[0] != pu.radiance_s:
                print("\nGenerating SNAPSHOTCY_LAT2 figures...")
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
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        iplot = 0
                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):


                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                            ymin = diagDF.loc[varLoc,statName].dropna().min()
                            ymax = diagDF.loc[varLoc,statName].dropna().max()

                            #subplot loop 2
                            for binVal in selectBinVals1D:

                                # collect statName for all lines on this subplot, letting cyDTime vary
                                linesVals = []
                                for expName in expNames:
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                               ['expName','fcTDelta','varName','binVal'])
                                    dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                    lineVals = np.empty(nCY)
                                    lineVals[:] = np.NaN
                                    for cyDTime in dataCYDTimes:
                                        icy = cyDTimes.index(cyDTime)
                                        lineVals[icy] = dataDF.loc[cyDTime]
                                    linesVals.append(lineVals)

                                # define subplot title
                                title = binVal+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries( fig, \
                                                cyDTimes, linesVals, expNames, \
                                                title, data_labels[istat], \
                                                sciticks[istat], signdef[istat], \
                                                ny, nx, nsubplots, iplot, \
                                                ymin = ymin, ymax = ymax )

                                iplot = iplot + 1

                        #expName_file = re.sub("\.","",re.sub("\s+","-",expName))
                        filename = OSFigPath+'/SNAPSHOTCY_LAT2_TSeries_%sday_%s_%s_%s'%( \
                                   fcTDeltas_dir[ifc],ObsSpaceName, \
                                   diagName,statName)

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                # end binName loop

            # end SNAPSHOTCY_LAT2


            selectBinVals1D = pu.goodFlagNames + pu.badFlagNames
            if "SNAPSHOTCY_QC" in plotTypes \
                and len(selectBinVals1D) > 0:
                print("\nGenerating SNAPSHOTCY_QC figures...")
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
                    fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                    iplot = 0

                    #subplot loop 1
                    for ivar, varName in enumerate(varNames):
                        #subplot loop 2
                        for expName in expNames:

                            # use specific y-axis limits for each varName
                            varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals1D)
                            ymin = diagDF.loc[varLoc,statName].dropna().min()
                            ymax = diagDF.loc[varLoc,statName].dropna().max()

                            # collect statName for all lines on this subplot, letting cyDTime vary
                            linesVals = []
                            for binVal in selectBinVals1D:
                                #                           cyDTime
                                dataLoc = (expName,fcTDelta,slice(None),varName,binVal)
                                dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                           ['expName','fcTDelta','varName','binVal'])
                                dataCYDTimes = dataDF.index.get_level_values('cyDTime')

                                lineVals = np.empty(nCY)
                                lineVals[:] = np.NaN
                                for cyDTime in dataCYDTimes:
                                    icy = cyDTimes.index(cyDTime)
                                    lineVals[icy] = dataDF.loc[cyDTime]
                                linesVals.append(lineVals)

                            # define subplot title
                            title = expName+"\n"+varName
                            if varUnitss[ivar] != pu.miss_s:
                                title = title+" ("+varUnitss[ivar]+")"

                            # perform subplot agnostic plotting (all expNames)
                            plotTimeSeries( fig, \
                                            cyDTimes, linesVals, selectBinVals1D, \
                                            title, data_labels[istat], \
                                            sciticks[istat], signdef[istat], \
                                            ny, nx, nsubplots, iplot, \
                                            ymin = ymin, ymax = ymax, \
                                            legend_inside = False )
                            iplot = iplot + 1

                    filename = OSFigPath+'/SNAPSHOTCY_QC_TSeries_%sday_%s_%s_%s'%( \
                               fcTDeltas_dir[ifc],ObsSpaceName, \
                               diagName,statName)

                    pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                # end binName loop

            # end SNAPSHOTCY_QC


            #=============
            # 2-D figures
            #=============
            for binVar2D in binVars2D:
                if "SNAPSHOTCY2D" in plotTypes \
                   or "AGGREGATEFC2D" in plotTypes \
                   or "AGGREGATEFCPROFILE" in plotTypes \
                   or "AGGREGATEFCPROFILE_DIFFCI" in plotTypes:
                    selectBinVals2D = []
                    for ibin, binVal in enumerate(binVals):
                        #Only numerical bins
                        if (pu.isfloat(binVal) or pu.isint(binVal)) \
                           and binVars[ibin] == binVar2D:
                            selectBinVals2D.append(binVal)

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
                    print("\nGenerating SNAPSHOTCY2D figures across "+binVar2D+"...")
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
                            fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                            iplot = 0

                            #subplot loop 1
                            for ivar, varName in enumerate(varNames):
                                # use specific c-axis limits for each varName
                                varLoc = (slice(None),slice(None),slice(None),varName,selectBinVals2D)
                                cmin = diagDF.loc[varLoc,statName].dropna().min()
                                cmax = diagDF.loc[varLoc,statName].dropna().max()

                                #subplot loop 2
                                # letting cyDTime and binVal vary
                                for expName in expNames:
                                    contourVals = np.empty((len(selectBinVals2D),nCY))
                                    contourVals[:,:] = np.NaN
                                    #                           cyDTime
                                    dataLoc = (expName,fcTDelta,slice(None),varName,selectBinVals2D)
                                    dataDF = diagDF.loc[dataLoc,statName].droplevel( \
                                               ['expName','fcTDelta','varName'])
                                    dataCYDTimes = dataDF.loc[(slice(None),selectBinVals2D[0]) \
                                                             ].index.get_level_values('cyDTime')
                                    for ibin, binVal in enumerate(selectBinVals2D):
                                        tmp = dataDF.loc[(slice(None),binVal)].to_numpy()
                                        for jcy, cyDTime in enumerate(dataCYDTimes):
                                            icy = cyDTimes.index(cyDTime)
                                            contourVals[ibin,icy] = tmp[jcy]

                                    # define subplot title
                                    title = expName+"\n"+varName
                                    if varUnitss[ivar] != pu.miss_s:
                                        title = title+" ("+varUnitss[ivar]+")"

                                    # perform subplot agnostic plotting (all expNames)
                                    plotTimeSeries2D( fig, \
                                                      cyDTimes, selectBinNumVals, contourVals, \
                                                      title, data_labels[istat], \
                                                      sciticks[istat], signdef[istat], \
                                                      ylabel, invert_y_axis, \
                                                      ny, nx, nsubplots, iplot, \
                                                      cmin = cmin, cmax = cmax )
                                    iplot = iplot + 1

                            filename = OSFigPath+'/SNAPSHOTCY2D_%s_TSeries_%sday_%s_%s_%s'%( \
                                       binVar2D,fcTDeltas_dir[ifc],ObsSpaceName, \
                                       diagName,statName)

                            pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                        # end statName loop

                    # end fcTDelta loop

                # end SNAPSHOTCY2D


                if "AGGREGATEFC2D" in plotTypes \
                    and len(selectBinVals2D) > 0 \
                    and nFC > 1:
                    print("\nGenerating AGGREGATEFC2D figures across "+binVar2D+"...")

                    # top-level 2-D plot settings
                    subplot_size = 2.4
                    aspect = 0.55
                    nx = nExp
                    ny = nVars
                    nsubplots = nx * ny

                    # apply aggregation over cyDTime via including all other free indices in groupby
                    aggDF = diagDF.groupby(
                            ['expName','fcTDelta','varName','binVal']).apply(pu.aggStatsDF)

                    #file loop 1
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        iplot = 0
                        #subplot loop
                        for ivar, varName in enumerate(varNames):
                            # use specific c-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,selectBinVals2D)
                            cmin = aggDF.loc[varLoc,statName].dropna().min()
                            cmax = aggDF.loc[varLoc,statName].dropna().max()

                            #collect aggregated statName, varying across fcTDelta+binVal
                            for expName in expNames:
                                contourVals = np.empty((len(selectBinVals2D),nFC))
                                contourVals[:,:] = np.NaN
                                #                  fcTDelta
                                dataLoc = (expName,slice(None),varName,selectBinVals2D)
                                dataDF = aggDF.loc[dataLoc,statName].droplevel( \
                                           ['expName','varName'])
                                dataFCTDeltas = dataDF.loc[(slice(None),selectBinVals2D[0]) \
                                                          ].index.get_level_values('fcTDelta')
                                for ibin, binVal in enumerate(selectBinVals2D):
                                    tmp = dataDF.loc[(slice(None),binVal)].to_numpy()
                                    for jfc, fcTDelta in enumerate(dataFCTDeltas):
                                        ifc = fcTDeltas.index(fcTDelta)
                                        contourVals[ibin,ifc] = tmp[jfc]

                                # define subplot title
                                title = expName+"\n"+varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"

                                # perform subplot agnostic plotting (all expNames)
                                plotTimeSeries2D( fig, \
                                                  fcTDeltas, selectBinNumVals, contourVals, \
                                                  title, fcdata_labels[istat], \
                                                  sciticks[istat], signdef[istat], \
                                                  ylabel, invert_y_axis, \
                                                  ny, nx, nsubplots, iplot, \
                                                  cmin = cmin, cmax = cmax )
                                iplot = iplot + 1

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFC2D_%s_TSeries_%s-%sday_%s_%s_%s'%( \
                                   binVar2D,fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   fcDiagName,statName)

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS)

                    # end statName loop

                # end AGGREGATEFC2D


                if "AGGREGATEFCPROFILE" in plotTypes \
                    and len(selectBinVals2D) > 0:
                    print("\nGenerating AGGREGATEFCPROFILE figures across "+binVar2D+"...")

                    # top-level 2-D plot settings
                    subplot_size = 1.2
                    aspect = 1.3
                    nx = min([nFC,MAX_FC_SUBFIGS])
                    ny = nVars
                    nsubplots = nx * ny

                    # apply aggregation over cyDTime via including all other free indices in groupby
                    aggDF = diagDF.groupby(
                            ['expName','fcTDelta','varName','binVal']).apply(pu.aggStatsDF)

                    #file loop 1
                    for istat, statName in enumerate(statNames):
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        iplot = 0

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            # use specific x-axis limits for each varName
                            varLoc = (slice(None),slice(None),varName,selectBinVals2D)
                            xmin = aggDF.loc[varLoc,statName].dropna().min()
                            xmax = aggDF.loc[varLoc,statName].dropna().max()

                            #subplot loop 2
                            for ifc, fcTDelta in enumerate(fcTDeltas):

                                #Setting to avoid over-crowding
                                if ifc > (MAX_FC_SUBFIGS-1): continue

                                #collect aggregated statNames, varying across fcTDelta
                                linesVals = []
                                for expName in expNames:
                                    linesVals.append(aggDF.loc[(expName,fcTDelta,varName,selectBinVals2D), \
                                                                    statName].to_numpy())

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                # perform subplot agnostic plotting (all expNames)
                                plotProfile( fig, \
                                             linesVals, selectBinVals2D, \
                                             expNames, \
                                             title, fcdata_labels[istat], \
                                             sciticks[istat], signdef[istat], \
                                             ylabel, invert_y_axis, \
                                             ny, nx, nsubplots, iplot, \
                                             xmin = xmin, xmax = xmax )

                                iplot = iplot + 1

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFCPROFILE_%s_TSeries_%s-%sday_%s_%s_%s'%( \
                                   binVar2D,fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   fcDiagName,statName)

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS, True)

                    # end statName loop

                # end AGGREGATEFCPROFILE


                if "AGGREGATEFCPROFILE_DIFFCI" in plotTypes \
                    and len(selectBinVals2D) > 0:
                    print("\nGenerating AGGREGATEFCPROFILE_DIFFCI figures across "+binVar2D+"...")

                    # top-level 2-D plot settings
                    subplot_size = 1.2
                    aspect = 1.3

                    nx = min([nFC,MAX_FC_SUBFIGS])
                    ny = nVars
                    nsubplots = nx * ny

                    #file loop 1
                    for statName in bootStrapStats:
                        # establish a new figure
                        fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                        iplot = 0

                        #subplot loop 1
                        for ivar, varName in enumerate(varNames):
                            #subplot loop 2
                            for ifc, fcTDelta in enumerate(fcTDeltas):

                                #Setting to avoid over-crowding
                                if ifc > (MAX_FC_SUBFIGS-1): continue

                                linesVals = {}
                                for trait in pu.ciTraits: linesVals[trait] = []
                                for expName in noncntrlExpNames:
                                    lineVals = {}
                                    for trait in pu.ciTraits: lineVals[trait] = []

                                    for binVal in selectBinVals2D:
                                        #                              cyDTime
                                        cntrlLoc = (cntrlName,fcTDelta,slice(None),varName,binVal)

                                        #                          cyDTime
                                        expLoc = (expName,fcTDelta,slice(None),varName,binVal)

                                        ciVals = pu.bootStrapClusterFunc( \
                                                     X = diagDF.loc[expLoc,:], \
                                                     Y = diagDF.loc[cntrlLoc,:], \
                                                     n_samples = 10000, \
                                                     statNames = [statName])

                                        for trait in pu.ciTraits:
                                            lineVals[trait].append(ciVals[statName][trait][0])

                                    for trait in pu.ciTraits:
                                        linesVals[trait].append(lineVals[trait])

                                # define subplot title
                                title = varName
                                if varUnitss[ivar] != pu.miss_s:
                                    title = title+" ("+varUnitss[ivar]+")"
                                title = title+" @ "+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+"days"

                                # use specific y-axis limits for each varName
                                xmin = np.nanmin(linesVals[pu.cimin])
                                xmax = np.nanmax(linesVals[pu.cimax])

                                # perform subplot agnostic plotting (all expNames)
                                plotProfile( fig, \
                                             linesVals[pu.cimean], selectBinVals2D, \
                                             noncntrlExpNames, \
                                             title, \
                                             statName+"("+fcDiagName+"): [EXP - CTRL]", \
                                             False, False, \
                                             ylabel, invert_y_axis, \
                                             ny, nx, nsubplots, iplot, \
                                             linesValsMinCI = linesVals[pu.cimin], \
                                             linesValsMaxCI = linesVals[pu.cimax], \
                                             xmin = xmin, xmax = xmax, \
                                             lineAttribOffset = 1)


                                iplot = iplot + 1

                        # save each figure
                        filename = OSFigPath+'/AGGREGATEFCPROFILE_DIFFCI_%s_TSeries_%s-%sday_%s_%s_%s'%( \
                                   binVar2D,fcTDeltas_dir[0],fcTDeltas_dir[-1],ObsSpaceName, \
                                   fcDiagName,statName)

                        pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS, True)

                    # end statName loop

                # end AGGREGATEFCPROFILE_DIFFCI

            # end binVars2D loop

        # end diagName loop

    # end ObsSpaceName loop


###############################################################################
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


###############################################################################
lenWarnProf = 0
nanWarnProf = 0
def plotProfile(fig, \
                linesVals, yVals, \
                linesLabel, \
                title="", xlabel="x", \
                sciticks=False, signdef=False, \
                ylabel="y", invert_y_axis=False, \
                ny=1, nx=1, nplots=1, iplot=0, \
                linesValsMinCI=None, linesValsMaxCI=None, \
                xmin=np.NaN, xmax=np.NaN, \
                lineAttribOffset=0, \
                legend_inside=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# yVals            - independent variable on y-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# xlabel           - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# signdef          - whether linesVals is positive/negative definite, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# xmin, xmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if all(np.isnan(lineVals)):
            global nanWarnProf
            if nanWarnProf==0:
                print("\nWARNING: skipping all-NaN data")
                print(title,xlabel,linesLabel[iline])
            nanWarnProf=nanWarnProf+1
            continue
        if len(lineVals)!=len(yVals):
            global lenWarnProf
            if lenWarnProf==0:
                print("\nWARNING: skipping data where len(x)!=len(y)")
                print(title,xlabel,linesLabel[iline])
            lenWarnProf=lenWarnProf+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColors[iline+lineAttribOffset]

        ax.plot(lineVals, yVals, \
                color=pColor, \
                label=linesLabel[iline], \
                ls=pu.plotLineStyles[iline+lineAttribOffset], \
                linewidth=0.5)
        nLines = nLines + 1
        plotVals.append(lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if not np.isnan(x) else -1.0 for x in significant])

            siginds = np.array([i for i,x in enumerate(significant) if x > 0.0],dtype=int)
            if len(siginds) > 0:
                ax.plot(np.array(lineVals)[siginds], np.array(yVals)[siginds], \
                        color=pColor, \
                        ls='', \
                        marker=pu.plotMarkers[iline+lineAttribOffset], \
                        markersize=1.5)
            ax.plot(linesValsMinCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(linesValsMaxCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax.xaxis.get_offset_text().set_fontsize(3)

    # add vertical zero line for unbounded quantities
    if not signdef:
        ax.plot([0., 0.], [yVals[0], yVals[-1]], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize x-limits
    minxval, maxxval = pu.get_clean_ax_limits(xmin,xmax,plotVals,signdef)
    if not np.isnan(minxval) and not np.isnan(maxxval):
        ax.set_xlim(minxval,maxxval)
        if maxxval-minxval < 1.0 or \
           maxxval-minxval > 100.0:
            ax.tick_params(axis='x',rotation=-35)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not FULL_SUBPLOT_LABELS \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if FULL_SUBPLOT_LABELS or ix == 0:
        ax.set_xlabel(xlabel,fontsize=4)
        ax.set_ylabel(ylabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if invert_y_axis:
        ax.invert_yaxis()

    ax.grid()

    return


###############################################################################
lenWarnTS=0
nanWarnTS=0
def plotTimeSeries(fig, \
                   xsDates, linesVals, \
                   linesLabel, \
                   title="", ylabel="", \
                   sciticks=False, signdef=False, \
                   ny=1, nx=1, nplots=1, iplot=0, \
                   linesValsMinCI=None, linesValsMaxCI=None, \
                   ymin=np.NaN, ymax=np.NaN, \
                   lineAttribOffset=0, \
                   legend_inside=True):

# ARGUMENTS
# fig              - matplotlib figure object
# xsDates          - date x-values (list/array or list of lists/arrays
#                                   of float seconds, dt.timedelta, dt.datetime)
# linesVals        - dependent variable (list of arrays)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# ylabel           - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# signdef          - whether linesVals is positive/negative definite, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# ymin, ymax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if all(np.isnan(lineVals)):
            global nanWarnTS
            if nanWarnTS==0:
                print("\nWARNING: skipping all-NaN data")
                print(title,ylabel,linesLabel[iline])
            nanWarnTS=nanWarnTS+1
            continue

        #float xVals
        if isinstance(xsDates[0],(list,np.ndarray)):
            xVals = pu.TDeltas2Seconds(xsDates[min([iline,len(xsDates)-1])])
        else:
            xVals = pu.TDeltas2Seconds(xsDates)

        if len(lineVals)!=len(xVals):
            global lenWarnTS
            if lenWarnTS==0:
                print("\nWARNING: skipping data where len(x)!=len(y)")
                print(title,ylabel,linesLabel[iline])
            lenWarnTS=lenWarnTS+1
            continue

        if iline == 0:
            minX = xVals[0]
            maxX = xVals[-1]
        else:
            minX = min([xVals[0], minX])
            maxX = max([xVals[-1], maxX])

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColors[iline+lineAttribOffset]

        ax.plot(xVals, lineVals, \
                label=linesLabel[iline], \
                color=pColor, \
                ls=pu.plotLineStyles[iline+lineAttribOffset], \
                linewidth=0.5)
        nLines = nLines + 1
        plotVals.append(lineVals)

        # Add shaded CI regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if not np.isnan(x) else -1.0 for x in significant])
            siginds = np.array([i for i,x in enumerate(significant) if x > 0.0],dtype=int)
            if len(siginds) > 0:
                ax.plot(np.array(xVals)[siginds], np.array(lineVals)[siginds], \
                        color=pColor, \
                        ls='', \
                        marker=pu.plotMarkers[iline+lineAttribOffset], \
                        markersize=1.5)
            ax.plot(xVals, linesValsMinCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(xVals, linesValsMaxCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    #axes settings
    if isinstance(xsDates[0],(list,np.ndarray)):
        pu.format_x_for_dates(ax, xsDates[0])
    else:
        pu.format_x_for_dates(ax, xsDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    if sciticks:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    # add horizontal zero line for unbounded quantities
    if not signdef:
        ax.plot([minX, maxX], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize y-limits
    minyval, maxyval = pu.get_clean_ax_limits(ymin,ymax,plotVals,signdef)
    if not np.isnan(minyval) and not np.isnan(maxyval):
        ax.set_ylim(minyval,maxyval)

    ax.grid()

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not FULL_SUBPLOT_LABELS \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if FULL_SUBPLOT_LABELS or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        nlcol = np.int(np.ceil(np.sqrt(nLines)))
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=nlcol)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)


    return


###############################################################################
def plotTimeSeries2D(fig, \
                     xDates, yVals, contourVals, \
                     title="", clabel="", \
                     sciticks=False, signdef=False, \
                     ylabel="y", invert_y_axis=False, \
                     ny=1, nx=1, nplots=1, iplot=0, \
                     cmin=np.NaN, cmax=np.NaN):

# ARGUMENTS
# fig           - matplotlib figure object
# xDates        - date x-values (array of float seconds, dt.timedelta, dt.datetime)
# yVals         - second independent variable
# contourVals   - dependent variable (2d array)

# title         - subplot title, optional
# clabel        - label for dependent variable, optional
# sciticks      - whether contourVals needs scientific formatting for ticks, optional
# signdef       - whether contourVals is positive/negative definite, optional
# ylabel        - label for yVals, optional
# invert_y_axis - whether to invert y-axis orientation, optional

# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# cmin, cmax    - min/max values of contourVals, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    xVals = pu.TDeltas2Seconds(xDates)

    # standardize c-limits
    mincval, maxcval = pu.get_clean_ax_limits(cmin,cmax,contourVals,signdef)
    if signdef:
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
    if not FULL_SUBPLOT_LABELS \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if FULL_SUBPLOT_LABELS or ix == 0:
        ax.set_ylabel(ylabel,fontsize=4)
    if FULL_SUBPLOT_LABELS or ix == nx-1:
        #colorbar
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(contourVals)
        if not np.isnan(mincval) and not np.isnan(maxcval):
            m.set_clim(mincval,maxcval)
        cb = plt.colorbar(m, ax=ax)
        #scientific formatting
        if sciticks:
            cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            cb.ax.yaxis.get_offset_text().set_fontsize(3)

        cb.ax.tick_params(labelsize=3)
        cb.set_label(clabel,fontsize=5)

    if invert_y_axis:
        ax.invert_yaxis()

    # optionally add a grid
    #ax.grid()

    return


###############################################################################
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
