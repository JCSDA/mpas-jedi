#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import binning_configs as bcs
from collections.abc import Iterable
from collections import defaultdict
from copy import deepcopy
import collections
import datetime as dt
import inspect
import logging
import multiprocessing as mp
import numpy as np
from pathlib import Path
import plot_utils as pu
import re
import os
import stat_utils as su
import StatisticsDatabase as sdb
import var_utils as vu

bootStrapStats = []
for x in su.sampleableAggStats:
    if x != 'Count': bootStrapStats.append(x)

#Select the stats for plotting
# options: see su.allFileStats
statNames = ['Count','Mean','RMS','STD']

## plot settings
figureFileType = 'pdf' #['pdf','png']

interiorLabels = True


def anWorkingDir(DiagSpace):
  return DiagSpace+'_analyses'

###################################
## Base class for all analysisTypes
###################################
class AnalyzeStatisticsBase():
    def __init__(self, db, analysisType):
        self.db = db
        self.DiagSpaceName = db.DiagSpaceName
        self.analysisType = analysisType

        self.logger = logging.getLogger(__name__+'.'+self.DiagSpaceName+'.'+self.analysisType)

        ## Extract useful variables from the database
        self.diagNames = db.diagNames
        self.allDiagNames = db.allDiagNames

        self.expNames = db.expNames
        self.nExp = len(self.expNames)
        self.cntrlExpName = db.cntrlExpName
        self.noncntrlExpNames = db.noncntrlExpNames

        self.fcTDeltas = db.fcTDeltas
        self.fcTDeltas_totmin = db.fcTDeltas_totmin
        self.fcMap = list(zip(self.fcTDeltas, self.fcTDeltas_totmin))
        self.nFC = len(self.fcTDeltas)

        self.cyDTimes = db.cyDTimes
        self.nCY = len(self.cyDTimes)

        self.varNames = db.varNames
        self.nVars = len(self.varNames)

        varLabels = []
        for (varName, varUnits) in zip(self.varNames, db.varUnitss):
            label = varName
            if varUnits != vu.miss_s:
                label = label+' ('+varUnits+')'
            varLabels.append(label)
        self.varMap = list(zip(self.varNames, varLabels))
        self.chlist = db.chlist

        self.allBinVals = db.allBinVals
        self.binNumVals = db.binNumVals
        self.binNumVals2DasStr = db.binNumVals2DasStr

        ## Establish default configuration
        self.blocking = False
        self.parallelism = False

        # TODO(JJG): decide if nproc is needed
        # nproc could be used to initialize workers within a single
        # analysisType to be distributed however they would be most useful.
        # That would allow for more granularity of computational effort
        # but requires that analysisType to be "blocking" due to the nature
        # of multiprocessing Pool's.
        # self.nproc = nproc

        self.requestAggDFW = False
        self.blankBinMethodFile = bu.identityBinMethod

        self.subWidth = 2.5
        self.subAspect = 1.0

        self.MAX_FC_SUBFIGS = 6
        self.MAX_FC_LINES = 6

        ## Setup paths
        CWD = os.getcwd()
        wd = CWD.split('/')[-1]
        DSSubDir = anWorkingDir(self.DiagSpaceName)
        self.DiagSpacePath = Path('./')
        if wd != DSSubDir:
            self.DiagSpacePath = self.DiagSpacePath/DSSubDir
        self.myFigPath = self.DiagSpacePath/self.analysisType
        self.myFigPath.mkdir(parents=True, exist_ok=True)

    def binMethodFile(self, binMethod, before = True):
        '''
        Format the binMethod for file naming
        '''
        binMethodFile = ''
        if binMethod != self.blankBinMethodFile:
            if before:
                binMethodFile = '_'+binMethod
            else:
                binMethodFile = binMethod+'_'

        return binMethodFile

    def fcName(self, diagName):
        '''
        Format the diagName for forecast analysisType's
        '''
        fcDiagName = diagName
        if self.fcTDeltas[-1] > dt.timedelta(0):
            fcDiagName = fcDiagName.replace('omb','omf')
            fcDiagName = fcDiagName.replace('bmo','fmo')
            fcDiagName = fcDiagName.replace('bak','fc')
        return fcDiagName

    def statPlotAttributes(self, diagName, statName):
        '''
        Define plotting attributes for the combination of diagName and statName
        '''
        fcDiagName = self.fcName(diagName)
        if statName == 'Count':
            statDiagLabel = statName
            fcstatDiagLabel = statName
        else:
            statDiagLabel = statName+'('+diagName+')'
            fcstatDiagLabel = statName+'('+fcDiagName+')'

        #This only applies to the unbounded quantities (omb, oma, ana/bak for velocity)
        signDefinite = True
        fcstatDiagDiffLabel = statName+'('+fcDiagName+'): [EXP - '+self.cntrlExpName+']'
        if statName == 'Mean':
            if ('omb' in diagName or 'oma' in diagName):
                signDefinite = False
                fcstatDiagDiffLabel = statName+': ['+self.cntrlExpName+' - EXP]'
            if ('bmo' in diagName or 'amo' in diagName):
                signDefinite = False
                fcstatDiagDiffLabel = statName+': [EXP - '+self.cntrlExpName+']'

        sciTicks = (statName == 'Count')

        return statDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite

    def UNIONcntrlANDexpCYDTimes(self, dfw, myLoc = {}):
        '''
        Determine the union of cyDTimes available between
        the control and each other experiment at each fcTDelta
        '''
        cntrlLoc = deepcopy(myLoc)
        cntrlLoc['expName'] = self.cntrlExpName

        expsCYDTimes = {}
        for fcTDelta in self.fcTDeltas:
            cntrlLoc['fcTDelta'] = fcTDelta
            cntrlCYDTimes = set(dfw.levels('cyDTime', cntrlLoc))

            expLoc = deepcopy(cntrlLoc)
            for expName in self.noncntrlExpNames:
                expLoc['expName'] = expName
                expCYDTimes = set(dfw.levels('cyDTime', expLoc))
                expsCYDTimes[(expName, fcTDelta)] = list(cntrlCYDTimes & expCYDTimes)
                if len(cntrlCYDTimes) != len(expCYDTimes):
                    self.logger.warning(self.cntrlExpName+' and '+expName+' have different number of CYDTimes at forecast length ', fcTDelta, ' Only using common CYDTimes for CI calculation.')
        return expsCYDTimes

    def analyze(self, workers = None):
        self.logger.info('analyze()')
        if self.blocking:
            # analyses with internal blocking
            self.analyze_()
        elif self.parallelism:
            # each analysis determines how to split external workers
            self.analyze_(workers)
        else:
            # divide workers acros analyses without internal parallelism/blocking
            workers.apply_async(self.analyze_)

    def analyze_(self, workers = None):
        '''
        stub method
        '''
        pass


def categoryBinValsAttributes(dfw, fullBinVar, binMethod, options):
    '''
    Utility function for providing an ordered list of 
    pairs of binVals and associated labels for
    category binMethods in the context of a DFWrapper
    '''

    binVar = vu.varDictObs[fullBinVar][1]

    dbSelect1DBinVals = dfw.levels('binVal')
    binUnitss = dfw.uniquevals('binUnits')
    #if (len(binUnitss) == 0 or
    #    len(dbSelect1DBinVals) == 1): return None, None
    assert (len(binUnitss) != 0 and len(dbSelect1DBinVals) > 0), 'ERROR: categoryBinValsAttributes received invalid binVar/binMethod'

    binUnits = binUnitss[0]

    # reorder select1DBinVals to match binMethod definition
    # TODO(JJG): clean up for readability
    tmp = deepcopy(bcs.binVarConfigs.get(
        fullBinVar,{}).get(
        binMethod,{}).get(
        'values', dbSelect1DBinVals))
    select1DBinVals = []
    if (not isinstance(tmp, Iterable) or
        isinstance(tmp, str)):
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
            t = ' @ '+binVar+'='+binVal
            if binUnits != vu.miss_s:
                t = t+' '+binUnits
        elif binVal in [bcs.goodFlagName]:
            t = ''
        else:
            t = ' @ '+binVal
        binTitles.append(t)

    binValsMap = list(zip(select1DBinVals, binTitles))

    return binValsMap


#=============
# 1-D figures
#=============
class CategoryBinMethodBase(AnalyzeStatisticsBase):
    '''
    Base class used to analyze statistics across binMethods with zero-dimensioned or
      category binValues, e.g., QC flag, named latitude band, cloudiness regime, surface type
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)
        self.parallelism = True

        # default binVar/binMethod combinations
        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarPrs, bu.PjetMethod): {},
            (vu.obsVarAlt, bu.altjetMethod): {},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {},
        }

    def subplotArrangement(self, binValsMap):
        # subplot configuration
        if len(binValsMap) > 1:
            nxplots = len(binValsMap)
            nyplots = self.nVars
            nsubplots = nxplots * nyplots
        else:
            nsubplots = self.nVars
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        return nxplots, nyplots, nsubplots

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)

        # TODO(JJG): another category base class for multiple diagNames on same subplot
        #      i.e., binVarDict has some attribute like 
        #            diagNameGroups = [['omb','oma'],['obs','bak','ana']]
        for diagName in self.db.diagNames:
            if diagName not in self.db.allDiagNames: continue

            diagLoc = {'diagName': diagName}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)

            for (fullBinVar, binMethod), options in self.binVarDict.items():
                binVar = vu.varDictObs[fullBinVar][1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue
                for statName in statNames:
                    if statName not in options.get('onlyStatNames', statNames): continue

                    self.logger.info(binVar+', '+binMethod+', '+statName)

                    if useWorkers:
                        workers.apply_async(self.innerloopsWrapper,
                            args = (diagName, fullBinVar, binMethod, statName, options))
                    else:
                        self.innerloopsWrapper(
                            diagName, fullBinVar, binMethod, statName, options)

    def innerloopsWrapper(self,
        diagName, fullBinVar, binMethod, statName, options):

        binVar = vu.varDictObs[fullBinVar][1]

        myLoc = {}
        myLoc['diagName'] = diagName
        myLoc['binVar'] = binVar
        myLoc['binMethod'] = binMethod

        # reducing to mydfwDict speeds extractions in innerloops
        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # include aggregated statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'],['cyDTime'])

        binValsMap = categoryBinValsAttributes(
            mydfwDict['dfw'], fullBinVar, binMethod, options)

        nxplots, nyplots, nsubplots = self.subplotArrangement(binValsMap)

        self.innerloops(
            mydfwDict, myLoc, statName, binValsMap, options,
            nsubplots, nxplots, nyplots)

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):
        pass


class CYAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between firstCycleDTime and lastCycleDTime
      for each forecast length between fcTDeltaFirst and fcTDeltaLast
      -  x-axis: cycle initial time
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and FC lead time (if applicable)
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.subWidth = 1.9
        self.subAspect = 0.75

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        myPath = self.myFigPath/myLoc['diagName']
        myPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            lineLoc['fcTDelta'] = fcTDelta
            axisLimitsLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                #subplot loop 2
                for binVal, binTitle in binValsMap:
                    lineLoc['binVal'] = binVal
                    axisLimitsLoc['binVal'] = binVal

                    # use common y-axis limits across axisLimitsLoc database locations
                    if statName == 'Count':
                        dmin = 0.
                    else:
                        dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                    dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                    # collect statName for all lines on this subplot
                    linesVals = []
                    for expName in self.expNames:
                        lineLoc['expName'] = expName
                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc(cyLoc, statName)
                        linesVals.append(lineVals)

                    # define subplot title
                    title = varLabel+binTitle

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        self.cyDTimes, linesVals, self.expNames,
                        title, bgstatDiagLabel,
                        sciTicks, signDefinite,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

                # end binVal loop

            # end varMap loop

            # save each figure
            filename = myPath/('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), fcTDelta_totmin,
                       self.DiagSpaceName, myLoc['diagName'], statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end fcMap loop


class FCAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing
      aggregated statistics for the period between firstCycleDTime and lastCycleDTime
      -  x-axis: forecast duration
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.requestAggDFW = True

        self.subWidth = 1.9
        self.subAspect = 0.9

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
        iplot = 0

        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            lineLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for binVal, binTitle in binValsMap:
                lineLoc['binVal'] = binVal
                axisLimitsLoc['binVal'] = binVal

                # use common y-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

                #collect aggregated statNames, varying across fcTDelta
                linesVals = []
                for expName in self.expNames:
                    lineLoc['expName'] = expName
                    lineFCTDeltas = dfwDict['agg'].levels('fcTDelta', lineLoc)

                    lineVals = np.full(self.nFC, np.NaN)
                    fcLoc = deepcopy(lineLoc)
                    for fcTDelta in lineFCTDeltas:
                        ifc = self.fcTDeltas.index(fcTDelta)
                        fcLoc['fcTDelta'] = fcTDelta
                        lineVals[ifc] = dfwDict['agg'].loc(fcLoc, statName)
                    linesVals.append(lineVals)

                # define subplot title
                title = varLabel+binTitle

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals, self.expNames,
                    title, fcstatDiagLabel,
                    sciTicks, signDefinite,
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

            # end statMap loop

        # end varMap loop

        # save each figure
        filename = myPath/('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)


class FCAxisExpLinesDiffCI(CategoryBinMethodBase):
    '''
    Similar to FCAxisExpLines, except
      - shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by bootStrapStats
      - confidence intervals (CI) are shown at each lead time
      -    line+shaded region: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        # OPTIONAL: implement fine-grained parallelism for bootStrapping
        #self.blocking = True

        self.subWidth = 1.9
        self.subAspect = 0.9

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = bootStrapStats

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2 or self.nExp < 2: return

        #if statName not in bootStrapStats: return
        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
        iplot = 0

        binValLoc = {}
        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            binValLoc['varName'] = varName

            #subplot loop 2
            for binVal, binTitle in binValsMap:
                binValLoc['binVal'] = binVal

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], binValLoc)

                cntrlLoc = deepcopy(binValLoc)
                cntrlLoc['expName'] = self.cntrlExpName

                # define subplot title
                title = varLabel+binTitle

                linesVals = defaultdict(list)
                for expName in self.noncntrlExpNames:
                    lineVals = defaultdict(list)

                    for fcTDelta in self.fcTDeltas:
                        cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                        expLoc = deepcopy(cntrlLoc)
                        expLoc['expName'] = expName

                        cntrlLoc['fcTDelta'] = fcTDelta
                        expLoc['fcTDelta'] = fcTDelta

                        X = tempdfw.loc(expLoc)
                        Y = tempdfw.loc(cntrlLoc)

                        ciVals = su.bootStrapClusterFunc(
                                     X, Y,
                                     n_samples = 10000,
                                     statNames = [statName])

                        for trait in su.ciTraits:
                            lineVals[trait] += [ciVals[statName][trait][0]]

                    for trait in su.ciTraits:
                        linesVals[trait].append(lineVals[trait])

                # use specific y-axis limits for each varName
                dmin = np.nanmin(linesVals[su.cimin])
                dmax = np.nanmax(linesVals[su.cimax])

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals[su.cimean],
                    self.noncntrlExpNames,
                    title,
                    fcstatDiagDiffLabel,
                    False, False,
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = interiorLabels)
                iplot = iplot + 1

            # end binValsMap loop

        # end varName loop

        # save each figure
        filename = myPath/('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)


class CYAxisFCLines(CategoryBinMethodBase):
    '''
    Similar to CYAxisExpLines, except
      each line is for a different forecast lead time and
      each experiment is in a different file
      -  x-axis: valid time of forecast
      -    line: per FC lead time
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and experiment
      - self.MAX_FC_LINES determines number of FC lead time lines to include
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.subWidth = 1.9
        self.subAspect = 0.75

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2 or self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for expName in self.expNames:
            lineLoc['expName'] = expName

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                #subplot loop 2
                for binVal, binTitle in binValsMap:
                    lineLoc['binVal'] = binVal
                    axisLimitsLoc['binVal'] = binVal

                    # use common y-axis limits across axisLimitsLoc database locations
                    if statName == 'Count':
                        dmin = 0.
                    else:
                        dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                    dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                    # collect statName for all lines on this subplot, letting cyDTime vary
                    xsVals = []
                    linesVals = []
                    self.fcTDeltas_labels = []
                    for fcTDelta in self.fcTDeltas:
                        lineLoc['fcTDelta'] = fcTDelta

                        # calculate valid time for x-axis
                        xVals = []
                        for cyDTime in self.cyDTimes:
                            xVals.append(cyDTime+fcTDelta)
                        xsVals.append(xVals)

                        #Setting to avoid over-crowding
                        if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_LINES-1): continue

                        self.fcTDeltas_labels.append(
                            pu.timeDeltaTicks(fcTDelta.total_seconds(),0))

                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc(cyLoc, statName)

                        linesVals.append(lineVals)

                    # define subplot title
                    title = varLabel+binTitle

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xsVals, linesVals, self.fcTDeltas_labels,
                        title, bgstatDiagLabel,
                        sciTicks, signDefinite,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)
                    iplot = iplot + 1

                # end binValsMap loop

            # end varMap loop

            expFileName = re.sub('\.', '', re.sub('\s+', '-', expName))
            filename = myPath/('%s%s_TSeries_%s_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), expFileName,
                       self.DiagSpaceName, fcDiagName, statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end expName loop


###########################################
## Figures with individual lines per binVal
###########################################

class BinValLinesAnalysisType(CategoryBinMethodBase):
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        # TODO(JJG): allow for multiple binMethods in one figure, such as 
        #self.binVarDict = {
        #    (vu.obsVarQC, [bu.goodQCMethod, bu.badQCMethod]): {
        #        'onlyStatNames': ['Count'],
        #    },
        #}
        # OR if this is the only case for which it's needed
        # TODO: replace badQCMethod with all-encompassing QC Method, e.g., allQCMethod
        self.binVarDict = {
            (vu.obsVarQC, bu.badQCMethod): {
                'onlyStatNames': ['Count'],
            },
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {},
        }

    def subplotArrangement(self, dummy):
        # subplot configuration
        return self.nExp, self.nVars, self.nExp * self.nVars


class CYAxisBinValLines(BinValLinesAnalysisType):
    '''
    Similar to CYAxisExpLines, except
      each line is for a different binVal (e.g., latitude band, cloudiness, etc.)
      -    line: binVals for named bins (e.g., NXTro, Tro, SXTro for latitude)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of statistic and forecast length
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.subWidth = 1.9
        self.subAspect = 0.75

    def innerloops(self,
        dfwDict, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)


        myPath = self.myFigPath/myLoc['diagName']
        myPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        binVals = []
        for binVal, binTitle in binValsMap: binVals.append(binVal)
        lineLoc['binVal'] = binVals

        axisLimitsLoc = deepcopy(lineLoc)

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            lineLoc['fcTDelta'] = fcTDelta
            axisLimitsLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                lineLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common y-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                #subplot loop 2
                for expName in self.expNames:
                    lineLoc['expName'] = expName

                    # collect statName for all lines on this subplot, letting cyDTime vary
                    linesVals = []
                    for binVal in binVals:
                        lineLoc['binVal'] = binVal
                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc(cyLoc, statName)

                        linesVals.append(lineVals)

                    # end binVal loop

                    # define subplot title
                    title = expName+'\n'+varLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        self.cyDTimes, linesVals, binVals,
                        title, bgstatDiagLabel,
                        sciTicks, signDefinite,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

                # end expName Loop

            # end varMap Loop

            filename = myPath/('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       myLoc['diagName'], statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end fcMap loop


# TODO(JJG): implement FCAxisBinValLines similar to FCAxisExpLines


#########################################################
## Figures with binVal on one axis, i.e., 2D and profiles
#########################################################
class OneDimBinMethodBase(AnalyzeStatisticsBase):
    '''
    Base class used to analyze statistics across binMethods with one-dimensional binValues
      that are assigned numerical values, e.g., altitude, pressure, latitude, cloud fraction
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)
        self.parallelism = True

        # default binVars
        self.binVarDict = {
            vu.obsVarAlt: {'profilefunc': bpf.plotProfile},
            vu.obsVarACI: {'profilefunc': bpf.plotSeries},
            vu.obsVarCldFrac: {'profilefunc': bpf.plotSeries},
            vu.obsVarGlint: {'profilefunc': bpf.plotSeries},
            vu.obsVarLandFrac: {'profilefunc': bpf.plotSeries},
            vu.obsVarLat: {'profilefunc': bpf.plotProfile},
            vu.obsVarLT: {'profilefunc': bpf.plotSeries},
            vu.obsVarPrs: {'profilefunc': bpf.plotProfile},
            vu.obsVarSCI: {'profilefunc': bpf.plotSeries},
            vu.obsVarSenZen: {'profilefunc': bpf.plotSeries},
        }

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)

        for diagName in self.db.diagNames:
            if diagName not in self.db.allDiagNames: continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagName})

            for fullBinVar, options in self.binVarDict.items():
                binVar = vu.varDictObs[fullBinVar][1]
                if (binVar not in diagBinVars): continue

                binVarLoc = {}
                binVarLoc['diagName'] = diagName
                binVarLoc['binVar'] = binVar
                binVarLoc['binVal'] = self.binNumVals2DasStr

                #Make figures for all binMethods
                binMethods = self.db.dfw.levels('binMethod', binVarLoc)
                for binMethod in binMethods:
                    for statName in statNames:
                        if statName not in options.get('onlyStatNames', statNames): continue

                        self.logger.info(binVar+', '+binMethod+', '+statName)

                        if useWorkers:
                            workers.apply_async(self.innerloopsWrapper,
                                args = (diagName, binVar, binMethod, statName, options))
                        else:
                            self.innerloopsWrapper(
                                diagName, binVar, binMethod, statName, options)

    def innerloopsWrapper(self,
        diagName, binVar, binMethod, statName, options):

        myLoc = {}
        myLoc['diagName'] = diagName
        myLoc['binVar'] = binVar
        myLoc['binVal'] = self.binNumVals2DasStr
        myLoc['binMethod'] = binMethod

        # reducing to mydfwDict speeds extractions in innerloops
        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # include aggregated statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'],['cyDTime'])

        ## Get all float/int binVals associated with binVar
        binVals = mydfwDict['dfw'].levels('binVal')
        binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

        # assume all bins represent same variable/units
        indepLabel = binVar
        if binUnits != vu.miss_s:
            indepLabel = indepLabel+' ('+binUnits+')'

        # bin info
        binNumVals = []
        for binVal in binVals:
            ibin = self.allBinVals.index(binVal)
            binNumVals.append(self.binNumVals[ibin])

            # invert independent variable axis for pressure bins
            pressure_dict = vu.varDictObs.get(vu.obsVarPrs,['',''])
            invert_ind_axis = (pressure_dict[1] == binVar)

        # sort bins by numeric value
        indices = list(range(len(binNumVals)))
        indices.sort(key=binNumVals.__getitem__)
        binNumVals = list(map(binNumVals.__getitem__, indices))
        binVals = list(map(binVals.__getitem__, indices))

        myBinConfigs = {
            'str': binVals,
            'values': binNumVals,
            'indepLabel': indepLabel,
            'invert_ind_axis': invert_ind_axis,
        }
        if len(binVals) < 2: return

        self.innerloops(
            mydfwDict, myLoc, statName, myBinConfigs, options)

    def innerloops(self,
        dfwDict, myLoc, statName, myBinConfigs, options):
        pass


class CYandBinValAxes2D(OneDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, statistic, and FC lead time
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.subWidth = 2.4
        self.subAspect = 0.65

    def innerloops(self,
        dfwDict, myLoc, statName, myBinConfigs, options):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        myPath = self.myFigPath/myLoc['diagName']
        myPath.mkdir(parents=True, exist_ok=True)

        nxplots = self.nExp
        nyplots = self.nVars
        nsubplots = nxplots * nyplots

        planeLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            planeLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in self.varMap:
                planeLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common c-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['dfw'].min(axisLimitsLoc, statName)
                dmax = dfwDict['dfw'].max(axisLimitsLoc, statName)

                #subplot loop 2
                # letting cyDTime and binVal vary
                for expName in self.expNames:
                    planeLoc['expName'] = expName
                    planeCYDTimes = dfwDict['dfw'].levels('cyDTime', planeLoc)

                    planeVals = np.full((len(myBinConfigs['str']), self.nCY), np.NaN)
                    binLoc = deepcopy(planeLoc)
                    for ibin, binVal in enumerate(myBinConfigs['str']):
                        binLoc['binVal'] = binVal
                        tmp = dfwDict['dfw'].loc(binLoc, statName).to_numpy()
                        for jcy, cyDTime in enumerate(planeCYDTimes):
                            icy = self.cyDTimes.index(cyDTime)
                            planeVals[ibin, icy] = tmp[jcy]

                    # define subplot title
                    title = expName+'\n'+varLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries2D(
                        fig,
                        self.cyDTimes, myBinConfigs['values'], planeVals,
                        title, bgstatDiagLabel,
                        sciTicks, signDefinite,
                        myBinConfigs['indepLabel'], myBinConfigs['invert_ind_axis'],
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = interiorLabels)

                    iplot = iplot + 1

            filename = myPath/('%s%s_BinValAxisTSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       myLoc['diagName'], statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end fcTDelta loop


class FCandBinValAxes2D(OneDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''

    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.requestAggDFW = True

        self.subWidth = 2.4
        self.subAspect = 0.55

    def innerloops(self,
        dfwDict, myLoc, statName, myBinConfigs, options):

        if self.nFC < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        nxplots = self.nExp
        nyplots = self.nVars
        nsubplots = nxplots * nyplots

        planeLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)

        iplot = 0
        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            planeLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            # use common c-axis limits across axisLimitsLoc database locations
            if statName == 'Count':
                dmin = 0.
            else:
                dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
            dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

            #subplot loop 2
            #collect aggregated statName, varying across fcTDelta+binVal
            for expName in self.expNames:
                planeLoc['expName'] = expName
                planeFCTDeltas = dfwDict['agg'].levels('fcTDelta', planeLoc)

                planeVals = np.full((len(myBinConfigs['str']), self.nFC), np.NaN)
                binLoc = deepcopy(planeLoc)
                for ibin, binVal in enumerate(myBinConfigs['str']):
                    binLoc['binVal'] = binVal
                    tmp = dfwDict['agg'].loc(binLoc, statName).to_numpy()
                    for jfc, fcTDelta in enumerate(planeFCTDeltas):
                        ifc = self.fcTDeltas.index(fcTDelta)
                        planeVals[ibin, ifc] = tmp[jfc]

                # define subplot title
                title = expName+'\n'+varLabel

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries2D(
                    fig,
                    self.fcTDeltas, myBinConfigs['values'], planeVals,
                    title, fcstatDiagLabel,
                    sciTicks, signDefinite,
                    myBinConfigs['indepLabel'], myBinConfigs['invert_ind_axis'],
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

        # save each figure
        filename = myPath/('%s%s_BinValAxisTSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)


class BinValAxisProfile(OneDimBinMethodBase):
    '''
    Similar to FCandBinValAxes2D, except
      - each vertical column of raster points is plotted as a profile on
        a separate set of axes instead of in 2-D color
      - therefore this is a valid plot even for a single forecast length (omb)
      -    line: per experiment
      - subplot: column by lead time, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
      - self.MAX_FC_SUBFIGS determines number of FC lead times to include
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.requestAggDFW = True

        self.subWidth = 1.2
        self.subAspect = 1.3

    def innerloops(self,
        dfwDict, myLoc, statName, myBinConfigs, options):

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = self.nVars
            nsubplots = nxplots * nyplots
        else:
            nsubplots = self.nVars
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        ptLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
        iplot = 0

        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            ptLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                ptLoc['fcTDelta'] = fcTDelta
                axisLimitsLoc['fcTDelta'] = fcTDelta

                # use common x-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                #collect aggregated statNames, varying across fcTDelta
                linesVals = []
                for expName in self.expNames:
                    ptLoc['expName'] = expName

                    lineVals = []
                    for binVal in myBinConfigs['str']:
                        ptLoc['binVal'] = binVal
                        lineVals.append(dfwDict['agg'].loc(ptLoc, statName))

                    linesVals.append(lineVals)

                # define subplot title
                title = varLabel+' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+'days'

                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals, myBinConfigs['values'],
                    self.expNames,
                    title, fcstatDiagLabel,
                    sciTicks, signDefinite,
                    myBinConfigs['indepLabel'], myBinConfigs['invert_ind_axis'],
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels)


                iplot = iplot + 1

        # save each figure
        filename = myPath/('%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)


class BinValAxisProfileDiffCI(OneDimBinMethodBase):
    '''
    Similar to BinValAxisProfile, except
      shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by bootStrapStats
      - confidence intervals (CI) are shown at each lead time and binVal
      -    line+shaded region: per experiment
      - subplot: column by lead time, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
      - self.MAX_FC_SUBFIGS determines number of FC lead times to include
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        # OPTIONAL: implement fine-grained parallelism for bootStrapping
        #self.blocking = True

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = bootStrapStats

        self.subWidth = 1.2
        self.subAspect = 1.3

    def innerloops(self,
        dfwDict, myLoc, statName, myBinConfigs, options):

        if self.nExp < 2: return

        if statName not in bootStrapStats: return
        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(myLoc['diagName'], statName)

        fcDiagName = self.fcName(myLoc['diagName'])
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = self.nVars
            nsubplots = nxplots * nyplots
        else:
            nsubplots = self.nVars
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
        iplot = 0

        fcLoc = {}
        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            fcLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                fcLoc['fcTDelta'] = fcTDelta

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], fcLoc)

                cntrlLoc = deepcopy(fcLoc)
                cntrlLoc['expName'] = self.cntrlExpName

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                linesVals = defaultdict(list)
                for expName in self.noncntrlExpNames:
                    lineVals = defaultdict(list)

                    cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                    for binVal in myBinConfigs['str']:
                        cntrlLoc['binVal'] = binVal
                        expLoc = deepcopy(cntrlLoc)
                        expLoc['expName'] = expName

                        X = tempdfw.loc(expLoc)
                        Y = tempdfw.loc(cntrlLoc)

                        ciVals = su.bootStrapClusterFunc(
                                     X, Y,
                                     n_samples = 10000,
                                     statNames = [statName])

                        for trait in su.ciTraits:
                            lineVals[trait] += [ciVals[statName][trait][0]]

                    for trait in su.ciTraits:
                        linesVals[trait].append(lineVals[trait])

                # define subplot title
                title = varLabel+' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+'days'

                # use specific y-axis limits for each varName
                dmin = np.nanmin(linesVals[su.cimin])
                dmax = np.nanmax(linesVals[su.cimax])

                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals[su.cimean], myBinConfigs['values'],
                    self.noncntrlExpNames,
                    title,
                    fcstatDiagDiffLabel,
                    False, False,
                    myBinConfigs['indepLabel'], myBinConfigs['invert_ind_axis'],
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = interiorLabels)

                iplot = iplot + 1

        # save each figure
        filename = myPath/('%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)


class BinValAxisPDF(AnalyzeStatisticsBase):
    '''
    Similar to BinValAxisProfile, except
      uses Count statistic to analyze a PDF across binVals
      -  x-axis: binVal
      -    line: per binMethod
      - subplot: combination of FC lead time and DiagSpace variable
      -    file: per experiment (if applicable)
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)
        # TODO(JJG): Make a generic version of bpf.plotPDF, which 
        # currently overlays a standard Gaussian model. That should
        # be a special case only for vu.obsVarNormErr.
        self.binVarDict = {
            vu.obsVarNormErr: {'pdffunc': bpf.plotPDF},
        }

        self.requestAggDFW = True

        self.subWidth = 1.2
        self.subAspect = 1.3

    def analyze_(self, workers = None):
        for diagName in self.db.diagNames:
            if diagName not in self.db.allDiagNames: continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagName})

            for fullBinVar, options in self.binVarDict:
                binVar = vu.varDictObs[fullBinVar][1]
                if binVar not in diagBinVars: continue

                myLoc = {}
                myLoc['diagName'] = diagName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.binNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                # include aggregated statistics when requested
                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'],['cyDTime'])

                ## Get all float/int binVals associated with binVar
                binMethods = mydfwDict['dfw'].levels('binMethod')
                binVals = mydfwDict['dfw'].levels('binVal')
                binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

                # assume all bins represent same variable/units
                indepLabel = binVar
                if binUnits != vu.miss_s:
                    indepLabel = indepLabel+' ('+binUnits+')'

                # bin info
                binNumVals = []
                for binVal in binVals:
                    ibin = self.allBinVals.index(binVal)
                    binNumVals.append(self.binNumVals[ibin])

                # sort bins by numeric value
                indices = list(range(len(binNumVals)))
                indices.sort(key=binNumVals.__getitem__)
                binNumVals = list(map(binNumVals.__getitem__, indices))
                binVals = list(map(binVals.__getitem__, indices))

                if len(binVals) < 2: continue

                self.logger.info('binVar=>'+binVar)

                fcDiagName = self.fcName(diagName)
                myPath = self.myFigPath/fcDiagName
                myPath.mkdir(parents=True, exist_ok=True)

                if self.nFC > 1:
                    nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
                    nyplots = self.nVars
                    nsubplots = nxplots * nyplots
                else:
                    nsubplots = self.nVars
                    nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
                    nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

                ptLoc = deepcopy(myLoc)

                #file loop 1
                for expName in self.expNames:
                    ptLoc['expName'] = expName

                    # establish a new figure
                    fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)

                    iplot = 0

                    #subplot loop 1
                    for (varName, varLabel) in self.varMap:
                        ptLoc['varName'] = varName

                        #subplot loop 2
                        for fcTDelta in self.fcTDeltas:
                            ptLoc['fcTDelta'] = fcTDelta

                            #Setting to avoid over-crowding
                            if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                            #collect aggregated statNames, varying across fcTDelta
                            linesVals = []
                            binMethodLabels = []
                            for binMethod in binMethods:
                                ptLoc['binMethod'] = binMethod

                                # if binMethod != bu.identityBinMethod: do something with bu.identityBinMethod
                                if binMethod == bu.identityBinMethod:
                                    binMethodLabels.append('ObsSpace')
                                else:
                                    binMethodLabels.append(binMethod)

                                lineVals = []
                                for binVal in binVals:
                                    ptLoc['binVal'] = binVal
                                    lineVals.append(dfwDict['agg'].loc(ptLoc,'Count'))

                                linesVals.append(lineVals)

                            # define subplot title
                            title = varLabel+' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+'days'

                            # perform subplot agnostic plotting (all expNames)
                            options['pdffunc'](
                                fig,
                                linesVals, binNumVals,
                                binMethodLabels,
                                title,
                                indepLabel,
                                nyplots, nxplots, nsubplots, iplot,
                                interiorLabels = interiorLabels)

                            iplot = iplot + 1

                        # end fcTDelta loop

                    # end varMap loop

                    # save each figure
                    filename = myPath/('%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                               binVar, self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                               self.DiagSpaceName, fcDiagName, expName))

                    pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                # end expName loop


# TODO: generalize as a sub-class of OneDimBinMethodBase
class BinValAxisStatsComposite(AnalyzeStatisticsBase):
    '''
    Similar to BinValAxisProfile, except
      all statistics (Count, Mean, RMS, STD) are placed on the same axis
      -  x-axis: binVal
      -    line: per statistic
      - subplot: per DiagSpace variable
      -    file: combination of FC lead time, experiment, and binMethod (if applicable)
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)
        self.binVarDict = {
            # TODO(JJG): Make a generic version of bpf.plotComposite, because 
            # bpf.plotfitRampComposite also provides parameters for a ramp fitting
            # function that may not be useful for binVars besides vu.obsVarSCI.
            vu.obsVarSCI: {'statsfunc': bpf.plotfitRampComposite},
        }

        self.requestAggDFW = True

        self.subWidth = 1.9
        self.subAspect = 0.9

    def analyze_(self, workers = None):
        for diagName in self.db.diagNames:
            if diagName not in self.db.allDiagNames: continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagName})

            for fullBinVar, options in self.binVarDict.items():
                binVar = vu.varDictObs[fullBinVar][1]
                if (binVar not in diagBinVars): continue

                myLoc = {}
                myLoc['diagName'] = diagName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.binNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                # include aggregated statistics when requested
                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'],['cyDTime'])

                ## Get all float/int binVals associated with binVar
                binMethods = mydfwDict['dfw'].levels('binMethod')
                binVals = mydfwDict['dfw'].levels('binVal')
                binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

                # assume all bins represent same variable/units
                indepLabel = binVar
                if binUnits != vu.miss_s:
                    indepLabel = indepLabel+' ('+binUnits+')'

                # bin info
                binNumVals = []
                for binVal in binVals:
                    ibin = self.allBinVals.index(binVal)
                    binNumVals.append(self.binNumVals[ibin])

                # sort bins by numeric value
                indices = list(range(len(binNumVals)))
                indices.sort(key=binNumVals.__getitem__)
                binNumVals = list(map(binNumVals.__getitem__, indices))
                binVals = list(map(binVals.__getitem__, indices))

                nBins = len(binVals)
                if nBins < 2: continue

                fcDiagName = self.fcName(diagName)
                myPath = self.myFigPath/fcDiagName
                myPath.mkdir(parents=True, exist_ok=True)

                nsubplots = self.nVars
                nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
                nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

                ptLoc = {}

                #file loop 1
                for binMethod in binMethods:
                    ptLoc['binMethod'] = binMethod

                    self.logger.info('binVar=>'+binVar+', binMethod=>'+binMethod)

                    #file loop 2
                    for expName in self.expNames:
                        ptLoc['expName'] = expName

                        #file loop 3
                        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
                            ptLoc['fcTDelta'] = fcTDelta

                            # establish a new figure
                            fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
                            iplot = 0

                            ERRParams = {}
                            ERRParams[self.DiagSpaceName] = {}

                            #subplot loop 1
                            for (varName, varLabel) in self.varMap:
                                ptLoc['varName'] = varName

                                #collect aggregated statNames, varying across fcTDelta
                                countsVals = np.full(nBins,0)
                                meansVals  = np.full(nBins, np.NaN)
                                rmssVals   = np.full(nBins, np.NaN)
                                stdsVals   = np.full(nBins, np.NaN)

                                for ibin, binVal in enumerate(binVals):
                                    ptLoc['binVal'] = binVal
                                    countsVals[ibin] = dfwDict['agg'].loc(ptLoc,'Count').to_numpy()
                                    meansVals[ibin] = dfwDict['agg'].loc(ptLoc,'Mean').to_numpy()
                                    rmssVals[ibin] = dfwDict['agg'].loc(ptLoc,'RMS').to_numpy()
                                    stdsVals[ibin] = dfwDict['agg'].loc(ptLoc,'STD').to_numpy()

                                # define subplot title
                                title = varLabel

                                # perform subplot agnostic plotting (all expNames)
                                FitParams = options['statsfunc'](
                                    fig,
                                    binNumVals,
                                    countsVals,
                                    meansVals,
                                    rmssVals,
                                    stdsVals,
                                    title,
                                    'STATS('+fcDiagName+')',
                                    indepLabel,
                                    nyplots, nxplots, nsubplots, iplot,
                                    interiorLabels = interiorLabels)

                                paramKey = self.chlist[iplot]
                                if paramKey == '': paramKey = varName
                                ERRParams[self.DiagSpaceName][(paramKey, binMethod)] = FitParams

                                iplot = iplot + 1

                            YAMLParams = {}
                            print('\n#For binning_params:')
                            for key in sorted(ERRParams[self.DiagSpaceName]):
                                print(binVar+"ErrParams['"+self.DiagSpaceName+"'][", key, "]   = ",
                                       ERRParams[self.DiagSpaceName][key]['bu'])
                                for param, val in ERRParams[self.DiagSpaceName][key]['YAML'].items():
                                    if param not in YAMLParams: YAMLParams[param] = []
                                    YAMLParams[param] += val
                            print('\n#For UFO YAML config:')
                            for param, val in YAMLParams.items():
                                print('#  '+param+':', val)

                            # save each figure
                            filename = myPath/('%s%s_BinValAxis_%smin_%s_%s_%s'%(
                                       binVar, self.binMethodFile(binMethod), fcTDelta_totmin,
                                       self.DiagSpaceName, fcDiagName, expName))

                            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels, True)

                    # end expName loop

                # end binMethod loop

            # end fullBinVar loop


#===========================
# Calculate gross statistics
#===========================
class GrossValues(AnalyzeStatisticsBase):
    '''
    Calculate gross statistics for specified category binMethods at first forecast length 
      NOTE: currently only calculates statistics at self.fcTDeltas[0]
            adjust minimum forecast length in order to calculate
            for non-zero forecast lengths, assuming those lengths
            are present in db
    '''
    def __init__(self, db, analysisType):
        super().__init__(db, analysisType)

        self.requestAggDFW = True

        # Force serial processing so that console output is contiguous
        # TODO(JJG): output to an ascii file and remove this line
        self.blocking = True

        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {},
        }

    def analyze_(self, workers = None):
        for diagName in self.db.diagNames:
            if diagName not in self.db.allDiagNames: continue

            diagLoc = {'diagName': diagName}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)

            for (fullBinVar, binMethod), options in self.binVarDict.items():
                binVar = vu.varDictObs[fullBinVar][1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue

                myLoc = {}
                myLoc['diagName'] = diagName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.binNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}
                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'],['cyDTime'])

                print(' Calculate gross statistics: binVar=>'+binVar+', binMethod=>'+binMethod)

                binValsMap = categoryBinValsAttributes(
                    mydfwDict['dfw'], fullBinVar, binMethod, options)

                print(' at FC length ', self.fcTDeltas[0])
                # Calculate gross statistics for this binVal
                statsLoc = {}
                statsLoc['fcTDelta'] = self.fcTDeltas[0]
                for binVal, binTitle in binValsMap:
                    statsLoc['binVal'] = binVal
                    GrossValues = {}
                    for varName in self.varNames:
                        statsLoc['varName'] = varName
                        for expName in self.expNames:
                            statsLoc['expName'] = expName
                            statsDFW = sdb.DFWrapper.fromLoc(mydfwDict['agg'], statsLoc)
                            for statName in statNames:
                                GrossValues[(statName, expName, varName)] = statsDFW.var(statName).to_numpy()
                    for expName in self.expNames:
                        print('Gross statistics for')
                        print('experiment=>'+expName)
                        if len(binValsMap) > 1:
                            print('binVal=>'+binVal)
                        print(' variables: ', self.varNames)

                        for statName in statNames:
                            print(statName)
                            tmp = np.asarray([])
                            for varName in self.varNames:
                                tmp = np.append(tmp, GrossValues[(statName, expName, varName)])
                            print(tmp)

AnalysisTypeDict = {
    #Derived from CategoryBinMethodBase(AnalyzeStatisticsBase)
    'CYAxisExpLines': CYAxisExpLines,
    'FCAxisExpLines': FCAxisExpLines,
    'FCAxisExpLinesDiffCI': FCAxisExpLinesDiffCI,
    'CYAxisFCLines': CYAxisFCLines,
    'CYAxisBinValLines': CYAxisBinValLines,
    #Derived from OneDimBinMethodBase(AnalyzeStatisticsBase)
    'CYandBinValAxes2D': CYandBinValAxes2D,
    'FCandBinValAxes2D': FCandBinValAxes2D,
    'BinValAxisProfile': BinValAxisProfile,
    'BinValAxisProfileDiffCI': BinValAxisProfileDiffCI,
    # TODO(JJG): TwoDimBinMethodBase(AnalyzeStatisticsBase)
    #'BinValAxes2D': BinValAxes2D,
    #Derived from AnalyzeStatisticsBase
    'BinValAxisPDF': BinValAxisPDF,
    'BinValAxisStatsComposite': BinValAxisStatsComposite,
    'GrossValues': GrossValues,
}

# NOTES:
# (1) FCAxis* types require non-zero forecast length
# (2) CYAxis* types require > 1 analysis cycle
# (3) CYAxisFCLines requires (1) and (2)
# (4) *DiffCI types require more than one experiment

def AnalysisFactory(db, analysisType):
    myClass = AnalysisTypeDict.get(analysisType, None)
    assert (myClass is not None and inspect.isclass(myClass)), \
        ('\n\nERROR: AnalysisFactory cannot construct ', analysisType, ' without instructions in AnalysisTypeDict')
    return myClass(db, analysisType)


class Analyses():
    def __init__(self, db, analysisTypes, nproc = 1):
        self.nproc = nproc
        self.analyses = []
        for anType in analysisTypes:
            self.analyses.append(AnalysisFactory(db, anType))
        self.logger = logging.getLogger(__name__+'.'+db.DiagSpaceName)
        self.logger.info('Analyses Constructed')

    def analyze(self):
        self.logger.info("Entering Analyses.analyze()")

        if self.nproc > 1:
            workers = mp.Pool(self.nproc)
        else:
            workers = None

        for an in self.analyses:
            an.analyze(workers)

        if workers is not None:
            workers.close()
            workers.join()

        self.logger.info("Exiting Analyses.analyze()")
