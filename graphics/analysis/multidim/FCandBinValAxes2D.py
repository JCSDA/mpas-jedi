#!/usr/bin/env python3

import basic_plot_functions as bpf
from copy import deepcopy
import diag_utils as du
import numpy as np
import plot_utils as pu
import stat_utils as su

from analysis.multidim.MultiDimBinMethodBase import MultiDimBinMethodBase
import analysis.StatisticsDatabase as sdb

class FCandBinValAxes2D(MultiDimBinMethodBase):
    '''
    Creates raster maps with binVar binVals on y-axis
      - only applicable to binned diagnostics (e.g., vertical dimension, latitude, zenith angle)
      - subplot: column by experiment, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
    '''

    requestAggDFW = True

    subplotWidth = 2.6
    subplotAspect = 0.75

    maxDiagnosticsPerAnalysis = 1

    cldFracTransform = None
 
    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nFC < 2: return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        # axes settings
        nxplots = self.nExp
        nyplots = nVarsLoc
        nsubplots = nxplots * nyplots

        xVals = self.fcTDeltas
        xLabel = 'Lead Time'

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

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

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
                dmin_abs = dfwDict['agg'].min(axisLimitsLoc, statName)
            dmax_abs = dfwDict['agg'].max(axisLimitsLoc, statName)

            # letting fcTDelta and binVal vary
            # extract control experiment
            cntrlAggLoc = deepcopy(planeLoc)
            cntrlAggLoc['expName'] = self.cntrlExpName
            cntrlPlaneFCTDeltas = dfwDict['agg'].levels('fcTDelta', cntrlAggLoc)
            cntrlAggPlaneVals = np.full((nBinVals, self.nFC), np.NaN)

            for ibin, binVal in enumerate(binStrVals):
                cntrlAggLoc['binVal'] = binVal
                aggFCStats = dfwDict['agg'].loc(cntrlAggLoc, statName).to_numpy()
                for jfc, fcTDelta in enumerate(cntrlPlaneFCTDeltas):
                    if jfc > len(aggFCStats)-1: continue
                    ifc = self.fcTDeltas.index(fcTDelta)
                    cntrlAggPlaneVals[ibin, ifc] = aggFCStats[jfc]

            # intermediate tempdfw reduces extraction time in inner loops
            tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], planeLoc)

            #subplot loop 2
            dmin_relative = np.NaN
            dmax_relative = np.NaN
            for expName in self.expNames:
                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel
                title = expName+'\n'+title

                fcstatDiagLabel = fcstatDiagLabel_abs
                sciTicks = sciTicks_abs
                logScale = logScale_abs
                centralValue = centralValue_abs
                dmin = dmin_abs
                dmax = dmax_abs

                planeVals = {}
                for trait in su.ciTraits:
                    planeVals[trait] = np.full_like(cntrlAggPlaneVals, np.NaN)

                if (expName == cntrlAggLoc['expName']):
                    planeVals[su.cimean] = deepcopy(cntrlAggPlaneVals)
                else:
                    expAggLoc = deepcopy(cntrlAggLoc)
                    expAggLoc['expName'] = expName

                    cntrlLoc = deepcopy(cntrlAggLoc)
                    expLoc = deepcopy(cntrlLoc)
                    expLoc['expName'] = expName

                    # letting fcTDelta and binVal vary
                    # extract this experiment
                    expPlaneFCTDeltas = dfwDict['agg'].levels('fcTDelta', expAggLoc)
                    expAggPlaneVals = np.full_like(cntrlAggPlaneVals, np.NaN)
                    for ibin, binVal in enumerate(binStrVals):
                        iy = ibin

                        cntrlLoc['binVal'] = binVal
                        expLoc['binVal'] = binVal

                        expAggLoc['binVal'] = binVal
                        aggFCStats = dfwDict['agg'].loc(expAggLoc, statName).to_numpy()

                        for jfc, fcTDelta in enumerate(expPlaneFCTDeltas):
                            if jfc > len(aggFCStats)-1: continue
                            ifc = self.fcTDeltas.index(fcTDelta)
                            ix = ifc
                            expAggPlaneVals[iy, ix] = aggFCStats[jfc]

                            # normalizing value for ratio
                            normalizingStat = cntrlAggPlaneVals[iy, ix]

                            cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                            expLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]

                            if useRelativeDifference:
                                if statName in su.bootStrapStats:
                                    expLoc['fcTDelta'] = fcTDelta
                                    cntrlLoc['fcTDelta'] = fcTDelta

                                    X = tempdfw.loc(expLoc)
                                    Y = tempdfw.loc(cntrlLoc)

                                    # confidence intervals
                                    ciVals = su.bootStrapClusterFunc(
                                                 X, Y,
                                                 n_samples = 10000,
                                                 statNames = [statName])

                                else:
                                    ciVals = {statName: {
                                        su.cimean: expAggPlaneVals[iy, ix] - normalizingStat,
                                        su.cimin: np.NaN,
                                        su.cimax: np.NaN,
                                    }}
                            else:
                                ciVals = {statName: {
                                    su.cimean: expAggPlaneVals[iy, ix],
                                    su.cimin: np.NaN,
                                    su.cimax: np.NaN,
                                }}

                            for trait in su.ciTraits:
                                t = float(ciVals[statName][trait])
                                # automatically generate relative difference plots for positive-semi-definite statistics
                                if useRelativeDifference:
                                  # divide by cntrlLoc aggregated statName
                                  t /= normalizingStat
                                  if self.relativeErrorType == 'one hundred centered':
                                    t += 1.0
                                  t *= 100.0
                                planeVals[trait][iy, ix] = t


                    # automatically generate relative difference plots for positive-semi-definite statistics
                    if useRelativeDifference:
                        sciTicks = False
                        logScale = False

                        notused, dmin_relative, dmax_relative, centralValue, label = self.relativeErrorFunction(
                          expAggPlaneVals,
                          cntrlAggPlaneVals,
                          dmin_relative,
                          dmax_relative,
                        )
                        dmin = dmin_relative
                        dmax = dmax_relative
                        fcstatDiagLabel = statName.replace('RMS','rms').replace('Mean','mean')+': '+label

                cLabel = fcstatDiagLabel

                # perform subplot agnostic plotting (all expNames)
                bpf.plot2D(
                    fig,
                    xVals, binNumVals, planeVals[su.cimean],
                    title, xLabel, binLabel, cLabel,
                    bpf.defaultIndepConfig,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    contourValsMinCI = planeVals[su.cimin],
                    contourValsMaxCI = planeVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    interiorLabels = self.interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = ('%s%s_BinValAxisTSeries_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.35, 0.55)
