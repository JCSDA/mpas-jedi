#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import predefined_configs as pconf
from collections.abc import Iterable
from collections import defaultdict
from copy import deepcopy
import collections
import datetime as dt
import diag_utils as du
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

## plot settings
figureFileType = 'pdf' #['pdf','png']

interiorLabels = True


def anWorkingDir(DiagSpace):
  return DiagSpace+'_analyses'

###################################
## Base class for all analysisTypes
###################################
class AnalysisBase():
    def __init__(self, db, analysisType, diagnosticGroupings = {}):
        self.analysisType = analysisType
        self.DiagSpaceName = db.DiagSpaceName
        self.diagnosticGroupings = diagnosticGroupings

        self.logger = logging.getLogger(__name__+'.'+self.DiagSpaceName+'.'+self.analysisType)

        ## Extract useful variables from the database
        self.db = db
        self.diagnosticConfigs = db.diagnosticConfigs
        self.availableDiagnostics = list(self.diagnosticConfigs.keys())

        self.expNames = self.db.expNames
        self.nExp = len(self.expNames)
        self.cntrlExpName = self.db.cntrlExpName
        self.noncntrlExpNames = self.db.noncntrlExpNames

        self.fcTDeltas = self.db.fcTDeltas
        self.fcTDeltas_totmin = self.db.fcTDeltas_totmin
        self.fcMap = list(zip(self.fcTDeltas, self.fcTDeltas_totmin))
        self.nFC = len(self.fcTDeltas)

        self.cyDTimes = self.db.cyDTimes
        self.nCY = len(self.cyDTimes)

        self.varNames = self.db.varNames
        self.nVars = len(self.varNames)

        varLabels = []
        for (varName, varUnits) in zip(self.varNames, self.db.varUnitss):
            label = varName
            if varUnits != vu.miss_s:
                label = label+' ('+varUnits+')'
            varLabels.append(label)
        self.varMap = list(zip(self.varNames, varLabels))
        self.chlist = self.db.chlist

        self.allBinVals = self.db.allBinVals
        self.binNumVals = self.db.binNumVals
        self.binNumVals2DasStr = self.db.binNumVals2DasStr

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

        self.requiredStatistics = []
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

    def fcName(self, diagnosticGroup):
        '''
        Format the diagnosticGroup for forecast analysisType's
        '''
        fcDiagName = diagnosticGroup
        if self.fcTDeltas[-1] > dt.timedelta(0):
            fcDiagName = fcDiagName.replace('omm','omf')
            fcDiagName = fcDiagName.replace('omb','omf')
            fcDiagName = fcDiagName.replace('bmo','fmo')
            fcDiagName = fcDiagName.replace('mmo','fmo')
            fcDiagName = fcDiagName.replace('hmo','fmo')
            fcDiagName = fcDiagName.replace('bak','fc')
        return fcDiagName

    def statPlotAttributes(self, diagnosticGroup, statName,
                           allDiagnosticNames = None):
        '''
        Define plotting attributes for the combination of diagnosticGroup and statName
        '''
        ommDiagnostics = ['omb', 'oma', 'omm', 'omf']
        mmoDiagnostics = ['bmo', 'amo', 'mmo', 'fmo', 'mmgfsan']
        truncateDiagnostics = ommDiagnostics+mmoDiagnostics
        diagnosticGroup_ = diagnosticGroup
        for diag in truncateDiagnostics:
            if diag in diagnosticGroup_:
                diagnosticGroup_ = diag

        fcDiagName = self.fcName(diagnosticGroup_)
        if statName in ['Count']+du.CRStatNames:
            statDiagLabel = statName
            fcstatDiagLabel = statName
        else:
            statDiagLabel = statName+'('+diagnosticGroup_+')'
            fcstatDiagLabel = statName+'('+fcDiagName+')'

        #These only apply to unbounded quantities (omb, oma, ana/bak for velocity, differences)
        signDefinite = True
        allDiagnosticNames_ = deepcopy(allDiagnosticNames)
        if allDiagnosticNames_ is None:
            cntrlDiagnosticName = diagnosticGroup_
            allDiagnosticNames_ = [diagnosticGroup_]
        else:
            cntrlDiagnosticName = allDiagnosticNames_[0]

        for diag in truncateDiagnostics:
            if diag in cntrlDiagnosticName:
                cntrlDiagnosticName = diag
            for idiag, adiag in enumerate(allDiagnosticNames_):
                if diag in adiag:
                    allDiagnosticNames_[idiag] = diag

        cntrlExpDiagnosticLabel = expDiagnosticLabel(self.cntrlExpName, cntrlDiagnosticName, allDiagnosticNames_)
        fcstatDiagDiffLabel = statName+'('+fcDiagName+'): [EXP - '+cntrlExpDiagnosticLabel+']'
        if statName == 'Mean':
            if diagnosticGroup_ in ommDiagnostics:
                signDefinite = False
                fcstatDiagDiffLabel = statName+': ['+cntrlExpDiagnosticLabel+' - EXP]'
            if diagnosticGroup_ in mmoDiagnostics:
                signDefinite = False
                fcstatDiagDiffLabel = statName+': [EXP - '+cntrlExpDiagnosticLabel+']'

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
            for expName in self.expNames:
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
        virtual method
        '''
        raise NotImplementedError()


def expDiagnosticLabel(expName, diagnosticName, allDiagnosticNames):
    if len(allDiagnosticNames) > 1:
        return expName+'-'+diagnosticName
    else:
        return expName


def categoryBinValsAttributes(dfw, fullBinVar, binMethod, options):
    '''
    Utility function for providing an ordered list of
    pairs of binVals and associated labels for
    category binMethods in the context of a DFWrapper
    '''

    binVar = vu.varDictAll[fullBinVar][1]

    dbSelect1DBinVals = dfw.levels('binVal')
    binUnitss = dfw.uniquevals('binUnits')
    #if (len(binUnitss) == 0 or
    #    len(dbSelect1DBinVals) == 1): return None, None
    assert (len(binUnitss) != 0 and len(dbSelect1DBinVals) > 0), 'ERROR: categoryBinValsAttributes received invalid binVar/binMethod'

    binUnits = binUnitss[0]

    # reorder select1DBinVals to match binMethod definition
    # TODO(JJG): clean up for readability
    tmp = deepcopy(pconf.binVarConfigs.get(
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
        else:
            t = ' @ '+binVal
        binTitles.append(t)

    binValsMap = list(zip(select1DBinVals, binTitles))

    return binValsMap


#=============
# 1-D figures
#=============
class CategoryBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with zero-dimensioned or
      category binValues, e.g., QC flag, named latitude band, cloudiness regime, surface type
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        self.parallelism = True
        self.maxBinVarTier = 2

        # default binVar/binMethod combinations
        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {'binVarTier': 1},
            (vu.obsVarLat, bu.latbandsMethod): {'binVarTier': 1},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {'binVarTier': 1},
#            (vu.modVarLat, bu.latbandsMethod): {'binVarTier': 1},
            (vu.noBinVar, bu.noBinMethod): {'binVarTier': 1},
            (vu.obsRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.modelRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.obsVarPrs, bu.PjetMethod): {'binVarTier': 3},
            (vu.obsVarAlt, bu.altjetMethod): {'binVarTier': 3},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {'binVarTier': 3},
        }
        self.maxDiagnosticsPerAnalysis = 10 // self.nExp

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

        # TODO(JJG): construct member Diagnostic objects (create new class) from
        #            diagnosticConfigs instead of referencing dictionary
        #            entries below.
        # TODO(JJG): use same color, vary line style across diagnosticGroupings
        diagnosticGrouped = {}
        for diag in self.availableDiagnostics:
            diagnosticGrouped[diag] = False

        diagnosticGroupings = deepcopy(self.diagnosticGroupings)
        for group in list(diagnosticGroupings.keys()):
            diags = diagnosticGroupings[group]
            if (len(diags) > self.maxDiagnosticsPerAnalysis or
               not set(diags).issubset(set(list(self.availableDiagnostics)))):
                del diagnosticGroupings[group]
                continue
            for diag in diags: diagnosticGrouped[diag] = True

        for diag in self.availableDiagnostics:
            if not diagnosticGrouped[diag]:
                diagnosticGroupings[diag] = [diag]

        for diagnosticGroup, diagnosticNames in diagnosticGroupings.items():
            if len(diagnosticNames) > self.maxDiagnosticsPerAnalysis: continue
            if len(set(diagnosticNames) & set(self.availableDiagnostics)) == 0: continue
            diagnosticConfigs = {}
            analysisStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                analysisStatistics = set(list(analysisStatistics) +
                                         diagnosticConfigs[diagnosticName]['analysisStatistics'])
            if not set(self.requiredStatistics).issubset(analysisStatistics): continue

            diagLoc = {'diagName': diagnosticNames}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)
            for (fullBinVar, binMethod), options in self.binVarDict.items():
                if options.get('binVarTier', 10) > self.maxBinVarTier: continue
                binVar = vu.varDictAll[fullBinVar][1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue

                self.logger.info(diagnosticGroup+', '+binVar+', '+binMethod)

                if useWorkers:
                    workers.apply_async(self.innerloopsWrapper,
                        args = (diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, analysisStatistics, options))
                else:
                    self.innerloopsWrapper(
                        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, analysisStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, analysisStatistics, options):

        binVar = vu.varDictAll[fullBinVar][1]

        # narrow mydfwDict by binVar and binMethod to reduce run-time and memory
        myLoc = {}
        myLoc['binVar'] = binVar
        myLoc['binMethod'] = binMethod

        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # aggregate statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
            sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], diagnosticConfigs)

        # further narrow mydfwDict by diagName
        # NOTE: derived diagnostics may require multiple diagName values;
        # can only narrow by diagName after aggregation
        myLoc['diagName'] = list(diagnosticConfigs.keys())
        for key in mydfwDict.keys():
            mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

        binValsMap = categoryBinValsAttributes(
            mydfwDict['dfw'], fullBinVar, binMethod, options)

        nxplots, nyplots, nsubplots = self.subplotArrangement(binValsMap)

        for statName in analysisStatistics:
            if statName not in options.get('onlyStatNames', analysisStatistics): continue

            self.innerloops(
                mydfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
                nsubplots, nxplots, nyplots)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):
        '''
        virtual method
        '''
        raise NotImplementedError()


class CYAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between firstCycleDTime and lastCycleDTime
      for each forecast length between fcTDeltaFirst and fcTDeltaLast
      -  x-axis: cycle initial time
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and FC lead time (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subWidth = 1.9
        self.subAspect = 0.75

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        myPath = self.myFigPath/diagnosticGroup
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
                    linesLabel = []
                    linesGroup = []
                    for expName in self.expNames:
                        lineLoc['expName'] = expName
                        for diagnosticName in myLoc['diagName']:
                            linesGroup.append(expName)
                            linesLabel.append(expDiagnosticLabel(
                                expName, diagnosticName, myLoc['diagName']))

                            lineLoc['diagName'] = diagnosticName

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
                        self.cyDTimes, linesVals, linesLabel,
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
                       self.DiagSpaceName, diagnosticGroup, statName))

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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subWidth = 1.9
        self.subAspect = 0.9

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
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
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    lineLoc['expName'] = expName
                    for diagnosticName in myLoc['diagName']:
                        linesGroup.append(expName)
                        linesLabel.append(expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))
                        lineLoc['diagName'] = diagnosticName

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
                    self.fcTDeltas, linesVals, linesLabel,
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

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
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return
        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return

        #if statName not in bootStrapStats: return
        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
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
                linesLabel = []
                linesGroup = []
                cntrlLoc['diagName'] = myLoc['diagName'][0]
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        linesGroup.append(expName)
                        linesLabel.append(expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        lineVals = defaultdict(list)
                        for fcTDelta in self.fcTDeltas:
                            cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                            cntrlLoc['fcTDelta'] = fcTDelta

                            expLoc = deepcopy(cntrlLoc)
                            expLoc['diagName'] = diagnosticName
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

                # use specific y-axis limits for each varName
                dmin = np.nanmin(linesVals[su.cimin])
                dmax = np.nanmax(linesVals[su.cimax])

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals[su.cimean],
                    linesLabel,
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subWidth = 1.9
        self.subAspect = 0.75

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2 or self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

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
                'binVarTier': 1,
            },
            (vu.obsVarQC, bu.allQCMethod): {
                'onlyStatNames': ['Count'],
                'binVarTier': 1,
            },
            (vu.obsVarLat, bu.latbandsMethod): {'binVarTier': 1},
#            (vu.modVarLat, bu.latbandsMethod): {'binVarTier': 1},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {'binVarTier': 1},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {'binVarTier': 3},
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subWidth = 1.9
        self.subAspect = 0.75

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName)


        myPath = self.myFigPath/diagnosticGroup
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
#                            print(dfwDict['dfw'].loc(cyLoc, statName))
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
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end fcMap loop


# TODO(JJG): implement FCAxisBinValLines similar to FCAxisExpLines


#########################################################
## Figures with binVal on one axis, i.e., 2D and profiles
#########################################################
class OneDimBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with one-dimensional binValues
      that are assigned numerical values, e.g., altitude, pressure, latitude, cloud fraction
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        self.parallelism = True
        self.maxBinVarTier = 1

        # default binVars
        self.binVarDict = {
            vu.obsVarAlt: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.obsVarACI: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
            vu.obsVarCldFrac: {'profilefunc': bpf.plotSeries, 'binVarTier': 1},
            vu.obsVarLat: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.obsVarPrs: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.obsVarSCI: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
#            vu.modVarLat: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.modVarLev: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.obsVarGlint: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLandFrac: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLT: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarSenZen: {'profilefunc': bpf.plotSeries, 'binbinVarTier': 3},
        }
        self.maxDiagnosticsPerAnalysis = 10 // self.nExp

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)
        diagnosticGrouped = {}
        for diag in self.availableDiagnostics:
            diagnosticGrouped[diag] = False

        diagnosticGroupings = deepcopy(self.diagnosticGroupings)
        for group in list(diagnosticGroupings.keys()):
            diags = diagnosticGroupings[group]
            if (len(diags) > self.maxDiagnosticsPerAnalysis or
               not set(diags).issubset(set(list(self.availableDiagnostics)))):
                del diagnosticGroupings[group]
                continue
            for diag in diags: diagnosticGrouped[diag] = True

        for diag in self.availableDiagnostics:
            if not diagnosticGrouped[diag]:
                diagnosticGroupings[diag] = [diag]

        for diagnosticGroup, diagnosticNames in diagnosticGroupings.items():
            if len(diagnosticNames) > self.maxDiagnosticsPerAnalysis: continue
            if len(set(diagnosticNames) & set(self.availableDiagnostics)) == 0: continue
            diagnosticConfigs = {}
            analysisStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                analysisStatistics = set(list(analysisStatistics) +
                                         diagnosticConfigs[diagnosticName]['analysisStatistics'])
            if not set(self.requiredStatistics).issubset(analysisStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticNames})

            for fullBinVar, options in self.binVarDict.items():
                if options.get('binVarTier', 10) > self.maxBinVarTier: continue
                binVar = vu.varDictAll[fullBinVar][1]
                if (binVar not in diagBinVars): continue

                binVarLoc = {}
                binVarLoc['diagName'] = diagnosticNames
                binVarLoc['binVar'] = binVar
                binVarLoc['binVal'] = self.binNumVals2DasStr

                #Make figures for all binMethods
                binMethods = self.db.dfw.levels('binMethod', binVarLoc)
                for binMethod in binMethods:

                    self.logger.info(diagnosticGroup+', '+binVar+', '+binMethod)

                    if useWorkers:
                        workers.apply_async(self.innerloopsWrapper,
                            args = (diagnosticGroup, diagnosticConfigs, binVar, binMethod, analysisStatistics, options))
                    else:
                        self.innerloopsWrapper(
                            diagnosticGroup, diagnosticConfigs, binVar, binMethod, analysisStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, binVar, binMethod, analysisStatistics, options):

        myLoc = {}
        myLoc['binVar'] = binVar
        myLoc['binVal'] = self.binNumVals2DasStr
        myLoc['binMethod'] = binMethod

        # narrow mydfwDict by binVar, binVal, and binMethod to reduce run-time and memory
        mydfwDict = {'dfw': self.db.loc(myLoc)}

        # aggregate statistics when requested
        if self.requestAggDFW:
            mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
            sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], diagnosticConfigs)

        # further narrow mydfwDict by diagName
        # NOTE: derived diagnostics may require multiple diagName values;
        # can only narrow by diagName after aggregation
        myLoc['diagName'] = list(diagnosticConfigs.keys())
        for key in mydfwDict.keys():
            mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

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
            pressure_dict = vu.varDictAll.get(vu.obsVarPrs,['',''])
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

        # only analyze variables that have non-zero Count when sliced by myLoc
        nVarsLoc = 0
        varMapLoc = []
        for (varName, varLabel) in self.varMap:
          countDF = mydfwDict['dfw'].loc({'varName': varName}, 'Count')
          if countDF.shape[0] > 0:
            counts = countDF.to_numpy()
            if np.nansum(counts) > 0:
              nVarsLoc += 1
              varMapLoc.append((varName, varLabel))

        for statName in analysisStatistics:
            if statName not in options.get('onlyStatNames', analysisStatistics): continue

            self.innerloops(
                mydfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):
        '''
        virtual method
        '''
        raise NotImplementedError()


class CYandBinValAxes2D(OneDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, statistic, and FC lead time
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.subWidth = 2.4
        self.subAspect = 0.65

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName)

        myPath = self.myFigPath/diagnosticGroup
        myPath.mkdir(parents=True, exist_ok=True)

        nxplots = self.nExp
        nyplots = nVarsLoc
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
            for (varName, varLabel) in varMapLoc:
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
                            if jcy > len(tmp)-1: continue
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
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(filename), figureFileType, interiorLabels)

        # end fcTDelta loop


class FCandBinValAxes2D(OneDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''

    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subWidth = 2.4
        self.subAspect = 0.55

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nFC < 2: return

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        planeLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)

        iplot = 0
        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
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
                        if jfc > len(tmp)-1: continue
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        self.subWidth = 1.2
        self.subAspect = 1.3

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = nVarsLoc
            nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        ptLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subWidth, self.subAspect, interiorLabels)
        iplot = 0

        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            ptLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                ptLoc['fcTDelta'] = fcTDelta
                axisLimitsLoc['fcTDelta'] = fcTDelta

#                if len(dfwDict['agg'].loc(axisLimitsLoc, statName)) == 0:
#                  iplot += 1
#                  continue

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
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    ptLoc['expName'] = expName
                    for diagnosticName in myLoc['diagName']:
                        linesGroup.append(expName)
                        linesLabel.append(expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        ptLoc['diagName'] = diagnosticName

                        lineVals = []
                        for binVal in myBinConfigs['str']:
                            ptLoc['binVal'] = binVal
                            pt = dfwDict['agg'].loc(ptLoc, statName).to_numpy()
                            if len(pt) == 1:
                              lineVals.append(pt[0])
                            else:
                              lineVals.append(np.NaN)

                        linesVals.append(lineVals)

                # define subplot title
                title = varLabel+' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+'days'

                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals, myBinConfigs['values'],
                    linesLabel,
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
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

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
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return

        if statName not in bootStrapStats: return
        bgstatDiagLabel, fcstatDiagLabel, fcstatDiagDiffLabel, sciTicks, signDefinite = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        myPath = self.myFigPath/fcDiagName
        myPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = nVarsLoc
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
        for (varName, varLabel) in varMapLoc:
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
                linesLabel = []
                linesGroup = []
                cntrlLoc['diagName'] = myLoc['diagName'][0]
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                        linesGroup.append(expName)
                        linesLabel.append(expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        lineVals = defaultdict(list)

                        for binVal in myBinConfigs['str']:
                            cntrlLoc['binVal'] = binVal
                            expLoc = deepcopy(cntrlLoc)
                            expLoc['diagName'] = diagnosticName
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
                    linesLabel,
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


class BinValAxisPDF(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      uses Count statistic to analyze a PDF across binVals
      -  x-axis: binVal
      -    line: per binMethod
      - subplot: combination of FC lead time and DiagSpace variable
      -    file: per experiment (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        # TODO(JJG): Make a generic version of bpf.plotPDF, which
        # currently overlays a standard Gaussian model. That should
        # be a special case only for vu.obsVarNormErr.
        self.binVarDict = {
            vu.obsVarNormErr: {'pdffunc': bpf.plotPDF},
        }

        self.requestAggDFW = True

        self.subWidth = 1.2
        self.subAspect = 1.3

        self.requiredStatistics = ['Count']

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            analysisStatistics = diagnosticConfig['analysisStatistics']
            if not set(self.requiredStatistics).issubset(analysisStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

            for fullBinVar, options in self.binVarDict:
                binVar = vu.varDictAll[fullBinVar][1]
                if binVar not in diagBinVars: continue

                myLoc = {}
                myLoc['diagName'] = diagnosticName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.binNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                # include aggregated statistics when requested
                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])

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

                fcDiagName = self.fcName(diagnosticName)
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
class BinValAxisStatsComposite(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      all statistics (Count, Mean, RMS, STD) are placed on the same axis
      -  x-axis: binVal
      -    line: per statistic
      - subplot: per DiagSpace variable
      -    file: combination of FC lead time, experiment, and binMethod (if applicable)
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)
        self.binVarDict = {
            # TODO(JJG): Make a generic version of bpf.plotComposite, because
            # bpf.plotfitRampComposite also provides parameters for a ramp fitting
            # function that may not be useful for binVars besides vu.obsVarSCI.
            vu.obsVarSCI: {'statsfunc': bpf.plotfitRampComposite},
        }

        self.requestAggDFW = True

        # Force serial processing so that console output is contiguous
        # TODO(JJG): output to an ascii file and remove this line
        self.blocking = True

        self.subWidth = 1.9
        self.subAspect = 0.9

        self.requiredStatistics = ['Count', 'Mean', 'RMS', 'STD']

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            analysisStatistics = diagnosticConfig['analysisStatistics']
            if not set(self.requiredStatistics).issubset(set(analysisStatistics)): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

            for fullBinVar, options in self.binVarDict.items():
                binVar = vu.varDictAll[fullBinVar][1]
                if (binVar not in diagBinVars): continue

                myLoc = {}
                myLoc['diagName'] = diagnosticName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.binNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                # include aggregated statistics when requested
                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])

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

                fcDiagName = self.fcName(diagnosticName)
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
                                    countsVals[ibin] = mydfwDict['agg'].loc(ptLoc,'Count').to_numpy()
                                    meansVals[ibin] = mydfwDict['agg'].loc(ptLoc,'Mean').to_numpy()
                                    rmssVals[ibin] = mydfwDict['agg'].loc(ptLoc,'RMS').to_numpy()
                                    stdsVals[ibin] = mydfwDict['agg'].loc(ptLoc,'STD').to_numpy()

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
class GrossValues(AnalysisBase):
    '''
    Calculate gross statistics for specified category binMethods at first forecast length
      NOTE: currently only calculates statistics at self.fcTDeltas[0]
            adjust minimum forecast length in order to calculate
            for non-zero forecast lengths, assuming those lengths
            are present in db
    '''
    def __init__(self, db, analysisType, diagnosticGroupings):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.requestAggDFW = True

        # Force serial processing so that console output is contiguous
        # TODO(JJG): output to an ascii file and remove this line
        self.blocking = True

        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarCldFrac, bu.cloudbandsMethod): {},
        }

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            analysisStatistics = diagnosticConfig['analysisStatistics']
            if not set(self.requiredStatistics).issubset(set(analysisStatistics)): continue

            diagLoc = {'diagName': diagnosticName}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)

            for (fullBinVar, binMethod), options in self.binVarDict.items():
                binVar = vu.varDictAll[fullBinVar][1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue

                # narrow mydfwDict by binVar and binMethod to reduce run-time and memory
                myLoc = {}
                myLoc['binVar'] = binVar
                myLoc['binMethod'] = binMethod

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
                    sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], {diagnosticName: diagnosticConfig})

                # further narrow mydfwDict by diagName
                # NOTE: derived diagnostics may require multiple diagName values;
                # can only narrow by diagName after aggregation
                myLoc['diagName'] = diagnosticName
                for key in mydfwDict.keys():
                    mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

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
                            for statName in analysisStatistics:
                                GrossValues[(statName, expName, varName)] = statsDFW.var(statName).to_numpy()
                    for expName in self.expNames:
                        print('Gross statistics for')
                        print('experiment=>'+expName)
                        if len(binValsMap) > 1:
                            print('binVal=>'+binVal)
                        print(' variables: ', self.varNames)

                        for statName in analysisStatistics:
                            print(statName)
                            tmp = np.asarray([])
                            for varName in self.varNames:
                                tmp = np.append(tmp, GrossValues[(statName, expName, varName)])
                            print(tmp)

AnalysisTypeDict = {
    #Derived from CategoryBinMethodBase(AnalysisBase)
    'CYAxisExpLines': CYAxisExpLines,
    'FCAxisExpLines': FCAxisExpLines,
    'FCAxisExpLinesDiffCI': FCAxisExpLinesDiffCI,
    'CYAxisFCLines': CYAxisFCLines,
    'CYAxisBinValLines': CYAxisBinValLines,
    #Derived from OneDimBinMethodBase(AnalysisBase)
    'CYandBinValAxes2D': CYandBinValAxes2D,
    'FCandBinValAxes2D': FCandBinValAxes2D,
    'BinValAxisProfile': BinValAxisProfile,
    'BinValAxisProfileDiffCI': BinValAxisProfileDiffCI,
    # TODO(JJG): TwoDimBinMethodBase(AnalysisBase)
    #'BinValAxes2D': BinValAxes2D,
    #Derived from AnalysisBase
    'BinValAxisPDF': BinValAxisPDF,
    'BinValAxisStatsComposite': BinValAxisStatsComposite,
    'GrossValues': GrossValues,
}

# NOTES:
# (1) FCAxis* types require non-zero forecast length
# (2) CYAxis* types require > 1 analysis cycle
# (3) CYAxisFCLines requires (1) and (2)
# (4) *DiffCI types require more than one experiment

def AnalysisFactory(db, analysisType, diagnosticGroupings):
    myClass = AnalysisTypeDict.get(analysisType, None)
    assert (myClass is not None and inspect.isclass(myClass)), \
        ('\n\nERROR: AnalysisFactory cannot construct ', analysisType, ' without instructions in AnalysisTypeDict')
    return myClass(db, analysisType, diagnosticGroupings)


class Analyses():
    def __init__(self, db, analysisTypes, diagnosticGroupings = {}, nproc = 1):
        self.nproc = nproc
        self.analyses = []
        for anType in analysisTypes:
            self.analyses.append(AnalysisFactory(db, anType, diagnosticGroupings))
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
