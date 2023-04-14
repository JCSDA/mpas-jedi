#!/usr/bin/env python3

import basic_plot_functions as bpf
from collections import defaultdict
from copy import deepcopy
import diag_utils as du
import numpy as np
import plot_utils as pu
import stat_utils as su

from analysis.category.CategoryBinMethodBase import CategoryBinMethodBase
import analysis.StatisticsDatabase as sdb

class FCAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between fcTDeltaFirst and fcTDeltaLast containing
      aggregated statistics for the period between firstCycleDTime and lastCycleDTime
      -  x-axis: forecast duration
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''

    requestAggDFW = True

    subplotWidth = 2.4
    subplotAspect = 0.75

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return
        assert myLoc['binMethod'] is not None, 'FCAxisExpLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
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
                        linesLabel.append(self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))
                        lineLoc['diagName'] = diagnosticName

                        lineFCTDeltas = dfwDict['agg'].levels('fcTDelta', lineLoc)

                        lineVals = np.full(self.nFC, np.NaN)
                        fcLoc = deepcopy(lineLoc)
                        for fcTDelta in lineFCTDeltas:
                            ifc = self.fcTDeltas.index(fcTDelta)
                            fcLoc['fcTDelta'] = fcTDelta
                            lineVals[ifc] = dfwDict['agg'].loc1(fcLoc, statName)
                        linesVals.append(lineVals)

                # define subplot title
                title = varLabel+binTitle

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals, linesLabel,
                    title, 'Lead Time', fcstatDiagLabel,
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = self.interiorLabels)

                iplot = iplot + 1

            # end statMap loop

        # end varMap loop

        # save figure
        filename = ('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)


class FCAxisExpLinesDiffCI(CategoryBinMethodBase):
    '''
    Similar to FCAxisExpLines, except
      - shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by su.bootStrapStats
      - confidence intervals (CI) are shown at each lead time
      -    line+shaded region: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar and statistic
    '''
    # used for percent ratio plots
    requestAggDFW = True

    subplotWidth = 2.4
    subplotAspect = 0.75

    requiredStatistics = ['Count']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2: return
        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return
        assert myLoc['binMethod'] is not None, 'FCAxisExpLinesDiffCI.innerloops: binMethod cannot be None'

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'], isDifferencePlot=True)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
        iplot = 0

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        binValLoc = {}
        #subplot loop 1
        for (varName, varLabel) in self.varMap:
            binValLoc['varName'] = varName

            #subplot loop 2
            for binVal, binTitle in binValsMap:
                binValLoc['binVal'] = binVal

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], binValLoc)
                normdfw = sdb.DFWrapper.fromLoc(dfwDict['agg'], binValLoc)

                cntrlLoc = deepcopy(binValLoc)
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlLoc['diagName'] = myLoc['diagName'][0]

                normLoc = deepcopy(cntrlLoc)

                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel
                title += binTitle

                linesVals = defaultdict(list)
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        linesGroup.append(expName)
                        linesLabel.append(self.expDiagnosticLabel(
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

                            # normalizing value for ratio
                            normLoc['fcTDelta'] = fcTDelta
                            normalizingStat = normdfw.loc1(normLoc, statName)

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                lineVals[trait] += [t]

                        for trait in su.ciTraits:
                            linesVals[trait].append(lineVals[trait])

                # use specific y-axis limits for each varName
                #dmin = np.nanmin(linesVals[su.cimin])
                #dmax = np.nanmax(linesVals[su.cimax])
                dmin = np.nanmin(linesVals[su.cimean])
                dmax = np.nanmax(linesVals[su.cimean])

                if not (np.isfinite(dmin) or np.isfinite(dmax)):
                  iplot = iplot + 1
                  continue

                fcstatDiagLabel = fcstatDiagLabel_abs
                sciTicks = sciTicks_abs
                logScale = logScale_abs
                centralValue = 0.

                if useRelativeDifference:

                  fcstatDiagLabel = self.relativeErrorLabeler()
                  fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+fcstatDiagLabel

                  #dmin, dmax = self.relativeErrorLimiter(dmin, dmax)

                  centralValue = self.relativeErrorCenter
                  sciTicks = False
                  logScale = False

                # perform subplot agnostic plotting (all expNames)
                bpf.plotTimeSeries(
                    fig,
                    self.fcTDeltas, linesVals[su.cimean],
                    linesLabel,
                    title, 'Lead Time', fcstatDiagLabel,
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = self.interiorLabels)

                iplot = iplot + 1

            # end binValsMap loop

        # end varName loop

        # save figure
        filename = ('%s%s_TSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)
