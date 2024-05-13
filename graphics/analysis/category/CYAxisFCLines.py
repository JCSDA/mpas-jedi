#!/usr/bin/env python3

import basic_plot_functions as bpf
from copy import deepcopy
import numpy as np
import plot_utils as pu
import re

from analysis.category.CategoryBinMethodBase import CategoryBinMethodBase
import analysis.StatisticsDatabase as sdb

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

    subplotWidth = 2.4
    subplotAspect = 0.75

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nFC < 2 or self.nCY < 2: return
        assert myLoc['binMethod'] is not None, 'CYAxisFCLines.innerloops: binMethod cannot be None'

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        axisLimitsLoc = {}

        #file loop 1
        for expName in self.expNames:
            lineLoc['expName'] = expName

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
                            pu.timeDeltaTicks([fcTDelta.total_seconds()]))

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
                        xsVals, linesVals, self.fcTDeltas_labels,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = self.interiorLabels)

                    iplot = iplot + 1

                # end binValsMap loop

            # end varMap loop

            expFileName = re.sub('\.', '', re.sub('\s+', '-', expName))
            filename = ('%s%s_TSeries_%s_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']), expFileName,
                       self.DiagSpaceName, fcDiagName, statName))

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)

        # end expName loop
