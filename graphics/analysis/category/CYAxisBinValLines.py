#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
from copy import deepcopy
import numpy as np
import plot_utils as pu
import var_utils as vu

from analysis.category.CategoryBinMethodBase import CategoryBinMethodBase
import analysis.StatisticsDatabase as sdb

class BinValLines(CategoryBinMethodBase):
    '''
    Figures with individual lines per binVal
    '''
    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.binVarDict = {
            (vu.obsVarQC, bu.badQCMethod): {
                'onlyStatNames': ['Count'],
            },
            (vu.obsVarQC, bu.allQCMethod): {
                'onlyStatNames': ['Count'],
            },
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarLat, bu.troplatbandsMethod): {},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #(vu.modVarLat, bu.latbandsMethod): {'binVarTier': 1},
            (vu.obsVarCldFracX, bu.cloudbandsMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
            (vu.obsVarLandFrac, bu.surfbandsMethod): {'binVarTier': 3},
        }

    def subplotArrangement(self, dummy):
        # subplot configuration
        return self.nExp, self.nVars, self.nExp * self.nVars


class CYAxisBinValLines(BinValLines):
    '''
    Similar to CYAxisExpLines, except
      each line is for a different binVal (e.g., latitude band, cloudiness, etc.)
      -    line: binVals for named bins (e.g., NXTro, Tro, SXTro for latitude)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of statistic and forecast length
    '''
    subplotWidth = 2.4
    subplotAspect = 0.75

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.maxDiagnosticsPerAnalysis = 1

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, binValsMap, options,
        nsubplots, nxplots, nyplots):

        if self.nCY < 2: return

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName)


        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        lineLoc = {}
        binStrVals = [binVal for binVal, binTitle in binValsMap]
        lineLoc['binVal'] = binStrVals

        axisLimitsLoc = deepcopy(lineLoc)

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
                    for binVal in binStrVals:
                        lineLoc['binVal'] = binVal
                        lineCYDTimes = dfwDict['dfw'].levels('cyDTime', lineLoc)

                        lineVals = np.full(self.nCY, np.NaN)
                        cyLoc = deepcopy(lineLoc)
                        for cyDTime in lineCYDTimes:
                            icy = self.cyDTimes.index(cyDTime)
                            cyLoc['cyDTime'] = cyDTime
                            lineVals[icy] = dfwDict['dfw'].loc1(cyLoc, statName)

                        linesVals.append(lineVals)

                    # end binVal loop

                    # define subplot title
                    title = expName+'\n'+varLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plotTimeSeries(
                        fig,
                        xVals, linesVals, binStrVals,
                        title, None, bgstatDiagLabel,
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = self.interiorLabels)

                    iplot = iplot + 1

                # end expName Loop

            # end varMap Loop

            filename = ('%s%s_TSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)

        # end fcMap loop
