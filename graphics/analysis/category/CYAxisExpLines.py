#!/usr/bin/env python3

import basic_plot_functions as bpf
from copy import deepcopy
import numpy as np
import plot_utils as pu

from analysis.category.CategoryBinMethodBase import CategoryBinMethodBase
import analysis.StatisticsDatabase as sdb

class CYAxisExpLines(CategoryBinMethodBase):
    '''
    Creates a timeseries figure between firstCycleDTime and lastCycleDTime
      for each forecast length between fcTDeltaFirst and fcTDeltaLast
      -  x-axis: cycle initial time
      -    line: per experiment
      - subplot: combination of DiagSpace variable and binVal
      -    file: combination of binVar, statistic, and FC lead time (if applicable)
    '''

    subplotWidth = 2.4
    subplotAspect = 0.75

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return
        assert myLoc['binMethod'] is not None, 'CYAxisExpLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            # calculate valid time for x-axis
            xVals = []
            for cyDTime in self.cyDTimes:
                xVals.append(cyDTime+fcTDelta)

            lineLoc['fcTDelta'] = fcTDelta
            axisLimitsLoc['fcTDelta'] = fcTDelta

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
                            linesLabel.append(self.expDiagnosticLabel(
                                expName, diagnosticName, myLoc['diagName']))

                            lineLoc['diagName'] = diagnosticName

                            lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                            lineVals = np.full(self.nCY, np.NaN)
                            cyLoc = deepcopy(lineLoc)
                            for cyDTime in lineCYDTimes:
                                icy = self.cyDTimes.index(cyDTime)
                                cyLoc['cyDTime'] = cyDTime
                                lineVals[icy] = dfwDict['dfw'].loc1(cyLoc, statName)
                            linesVals.append(lineVals)

                    # define subplot title
                    title = varLabel+binTitle

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xVals, linesVals, linesLabel,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = self.interiorLabels)

                    iplot = iplot + 1

                # end binVal loop

            # end varMap loop

            # save each figure
            filename = ('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), fcTDelta_totmin,
                       self.DiagSpaceName, diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)

        # end fcMap loop
