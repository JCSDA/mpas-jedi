#!/usr/bin/env python3

import basic_plot_functions as bpf
from copy import deepcopy
import diag_utils as du
import numpy as np
import plot_utils as pu
import stat_utils as su

from analysis.multidim.MultiDimBinMethodBase import MultiDimBinMethodBase
import analysis.StatisticsDatabase as sdb

class CYandBinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, statistic, and FC lead time
    '''

    subplotWidth = 2.5
    subplotAspect = 0.70
    maxDiagnosticsPerAnalysis = 1

    cldFracTransform = 'logit'

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nCY < 2: return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName)

        figPath = self.myFigPath/diagnosticGroup
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # axes settings
        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        xLabel = None

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        nBinVals = len(binStrVals)
        binLabel = myBinConfigs['binLabel']

        planeLoc = {}
        axisLimitsLoc = {}

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        #file loop 1
        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
            # calculate valid time for x-axis
            xVals = []
            for cyDTime in self.cyDTimes:
                xVals.append(cyDTime+fcTDelta)

            planeLoc['fcTDelta'] = fcTDelta

            # establish a new figure
            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)

            iplot = 0

            #subplot loop 1
            for (varName, varLabel) in varMapLoc:
                planeLoc['varName'] = varName
                axisLimitsLoc['varName'] = varName

                # use common c-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin_abs = 0.
                else:
                    dmin_abs = dfwDict['dfw'].min(axisLimitsLoc, statName)
                dmax_abs = dfwDict['dfw'].max(axisLimitsLoc, statName)

                # letting cyDTime and binVal vary
                # extract control experiment
                cntrlLoc = deepcopy(planeLoc)
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlPlaneCYDTimes = dfwDict['dfw'].levels('cyDTime', cntrlLoc)
                cntrlPlaneVals = np.full((nBinVals, self.nCY), np.NaN)

                for ibin, binVal in enumerate(binStrVals):
                    cntrlLoc['binVal'] = binVal
                    tmp = dfwDict['dfw'].loc(cntrlLoc, statName).to_numpy()
                    for jcy, cyDTime in enumerate(cntrlPlaneCYDTimes):
                        if jcy > len(tmp)-1: continue
                        icy = self.cyDTimes.index(cyDTime)
                        cntrlPlaneVals[ibin, icy] = tmp[jcy]

                #subplot loop 2
                dmin_relative = np.NaN
                dmax_relative = np.NaN
                for expName in self.expNames:
                    expLoc = deepcopy(planeLoc)
                    expLoc['expName'] = expName

                    # define subplot title
                    if useRelativeDifference:
                        title = varName
                    else:
                        title = varLabel
                    title = expName+'\n'+title

                    bgstatDiagLabel = bgstatDiagLabel_abs
                    sciTicks = sciTicks_abs
                    logScale = logScale_abs
                    centralValue = centralValue_abs
                    dmin = dmin_abs
                    dmax = dmax_abs
                    if (expName == cntrlLoc['expName']):
                        planeVals = deepcopy(cntrlPlaneVals)
                    else:
                        # letting cyDTime and binVal vary
                        # extract this experiment
                        expPlaneCYDTimes = dfwDict['dfw'].levels('cyDTime', expLoc)
                        expPlaneVals = np.full_like(cntrlPlaneVals, np.NaN)

                        for ibin, binVal in enumerate(binStrVals):
                            expLoc['binVal'] = binVal
                            tmp = dfwDict['dfw'].loc(expLoc, statName).to_numpy()
                            for jcy, cyDTime in enumerate(expPlaneCYDTimes):
                                if jcy > len(tmp)-1: continue
                                icy = self.cyDTimes.index(cyDTime)
                                expPlaneVals[ibin, icy] = tmp[jcy]

                        # automatically generate relative difference plots for positive-semi-definite statistics
                        if useRelativeDifference:
                            sciTicks = False
                            logScale = False

                            planeVals, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                              expPlaneVals,
                              cntrlPlaneVals,
                              dmin_relative,
                              dmax_relative,
                            )
                            dmin = dmin_relative
                            dmax = dmax_relative
                            bgstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label
                        else:
                            planeVals = deepcopy(expPlaneVals)

                    cLabel = bgstatDiagLabel

                    # perform subplot agnostic plotting (all expNames)
                    bpf.plot2D(
                        fig,
                        xVals, binNumVals, planeVals,
                        title, xLabel, binLabel, cLabel,
                        bpf.defaultIndepConfig,
                        myBinConfigs['binConfig'],
                        sciTicks, logScale, centralValue,
                        nyplots, nxplots, nsubplots, iplot,
                        dmin = dmin, dmax = dmax,
                        interiorLabels = self.interiorLabels)

                    iplot = iplot + 1

            filename = ('%s%s_BinValAxisTSeries_%smin_%s_%s_%s'%(
                       myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                       fcTDelta_totmin, self.DiagSpaceName,
                       diagnosticGroup, statName))

            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)

        # end fcTDelta loop
