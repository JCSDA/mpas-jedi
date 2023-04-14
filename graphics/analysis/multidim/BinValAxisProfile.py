#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy
import diag_utils as du
import numpy as np
import plot_utils as pu
import stat_utils as su
import yaml

from analysis.multidim.MultiDimBinMethodBase import MultiDimBinMethodBase
import analysis.StatisticsDatabase as sdb

class BinValAxisProfile(MultiDimBinMethodBase):
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

    requestAggDFW = True

    subplotWidth = 1.6
    subplotAspect = 1.0

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        bgstatDiagLabel, fcstatDiagLabel, sciTicks, logScale, centralValue = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'])

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = np.nanmax([nVarsLoc, 1])
            nxplots = np.nanmax([np.int(np.ceil(np.sqrt(nsubplots))), 1])
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        ptLoc = {}
        axisLimitsLoc = {}

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
        iplot = 0

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        binLabel = myBinConfigs['binLabel']

        figureData = {}
        figureData['binNumVals'] = [float(f) for f in binNumVals]
        figureData['binLabel'] = binLabel
        figureData['dataLabel'] = fcstatDiagLabel
        figureData['subplots'] = []

        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            subplotData = {}
            subplotData['varName'] = str(varName)

            ptLoc['varName'] = varName
            axisLimitsLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                tDeltaDays = float(fcTDelta.total_seconds()) / 3600.0 / 24.0
                subplotData['tDeltaDays'] = tDeltaDays

                ptLoc['fcTDelta'] = fcTDelta
                axisLimitsLoc['fcTDelta'] = fcTDelta

                # define subplot title
                title = varLabel
                if len(self.fcTDeltas) > 1:
                  title += ' @ '+str(tDeltaDays)+' days'

                # use common x-axis limits across axisLimitsLoc database locations
                if statName == 'Count':
                    dmin = 0.
                else:
                    dmin = dfwDict['agg'].min(axisLimitsLoc, statName)
                dmax = dfwDict['agg'].max(axisLimitsLoc, statName)

                subplotData['title'] = title
                #subplotData['iplot'] = iplot
                subplotData['dmin'] = self.dataYAMLFmtFloat(dmin)
                subplotData['dmax'] = self.dataYAMLFmtFloat(dmax)

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                #collect aggregated statNames, varying across fcTDelta
                linesVals = []
                linesLabel = []
                linesGroup = []
                linesData = {}
                for expName in self.expNames:
                    linesData[str(expName)] = {}
                    ptLoc['expName'] = expName
                    for diagnosticName in myLoc['diagName']:
                        linesGroup.append(expName)
                        lineLabel = self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName'])
                        linesLabel.append(lineLabel)

                        ptLoc['diagName'] = diagnosticName

                        lineVals = []
                        for binVal in binStrVals:
                            ptLoc['binVal'] = binVal
                            pt = dfwDict['agg'].loc1(ptLoc, statName)
                            lineVals.append(pt)
                            #if len(pt) == 1:
                              #lineVals.append(pt[0])
                            #else:
                              #lineVals.append(np.NaN)
                        linesVals.append(lineVals)
                        key = str(diagnosticName)+' '+self.DiagSpaceName+str(myLoc['binMethod'])+'_'+str(varName)
                        linesData[str(expName)][key] = []
                        for d in lineVals: linesData[str(expName)][key].append(self.dataYAMLFmtFloat(d))

                subplotData['linesData'] = linesData

                figureData['subplots'].append(deepcopy(subplotData))

                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals, binNumVals,
                    linesLabel,
                    title, binLabel, fcstatDiagLabel,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = self.interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = '%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName)

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.50, 0.55)

        # save figure data as yaml
        figureYAML = yaml.safe_dump(
          figureData,
          indent=2,
          width=2147483647,
          allow_unicode=False,
          default_flow_style=None,
        )
        with open(str(dataPath/filename)+'.yaml', 'w') as file:
           file.write(figureYAML)


class BinValAxisProfileDiffCI(MultiDimBinMethodBase):
    '''
    Similar to BinValAxisProfile, except
      shows difference between experiment(s) and control
      - control is selected using cntrlExpIndex
      - statistics are narrowed down by su.bootStrapStats
      - confidence intervals (CI) are shown at each lead time and binVal
      -    line+shaded region: per experiment
      - subplot: column by lead time, row by DiagSpace variable
      -    file: combination of binVar, binMethod, and statistic
      - self.MAX_FC_SUBFIGS determines number of FC lead times to include
    '''

    # used for percent ratio plots
    requestAggDFW = True

    subplotWidth = 1.6
    subplotAspect = 1.0

    requiredStatistics = ['Count']

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        for key in self.binVarDict:
            if 'onlyStatNames' in self.binVarDict[key]:
                self.binVarDict[key]['onlyStatNames'] += su.bootStrapStats
            else:
                self.binVarDict[key]['onlyStatNames'] = su.bootStrapStats

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):

        if self.nExp * len(myLoc['diagName']) < 2: return
        if self.cntrlExpName not in dfwDict['dfw'].levels('expName'): return

        bgstatDiagLabel_abs, fcstatDiagLabel_abs, sciTicks_abs, logScale_abs, centralValue_abs = \
            self.statPlotAttributes(diagnosticGroup, statName, myLoc['diagName'], isDifferencePlot=True)

        fcDiagName = self.fcName(diagnosticGroup)
        figPath = self.myFigPath/fcDiagName
        figPath.mkdir(parents=True, exist_ok=True)
        dataPath = figPath/'data'
        dataPath.mkdir(parents=True, exist_ok=True)

        if self.nFC > 1:
            nxplots = min([self.nFC, self.MAX_FC_SUBFIGS])
            nyplots = nVarsLoc
            nsubplots = nxplots * nyplots
        else:
            nsubplots = np.nanmax([nVarsLoc, 1])
            nxplots = np.nanmax([np.int(np.ceil(np.sqrt(nsubplots))), 1])
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        # Only bootstrap over the union of cyDTimes available
        # from both experiments at each fcTDelta
        myExpsCYDTimes = self.UNIONcntrlANDexpCYDTimes(dfwDict['dfw'])

        # establish a new figure
        fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
        iplot = 0

        binNumVals = myBinConfigs['num']
        binStrVals = myBinConfigs['str']
        binLabel = myBinConfigs['binLabel']

        useRelativeDifference = (
          statName in su.posSemiDefiniteStats and
          self.relativeErrorType != 'disable' and
          not set(myLoc['diagName']).issubset(du.absoluteOnlyDiagnostics))

        fcLoc = {}
        #subplot loop 1
        for (varName, varLabel) in varMapLoc:
            fcLoc['varName'] = varName

            #subplot loop 2
            for fcTDelta in self.fcTDeltas:
                fcLoc['fcTDelta'] = fcTDelta

                # intermediate tempdfw reduces extraction time in inner loops
                tempdfw = sdb.DFWrapper.fromLoc(dfwDict['dfw'], fcLoc)
                normdfw = sdb.DFWrapper.fromLoc(dfwDict['agg'], fcLoc)

                cntrlLoc = deepcopy(fcLoc)
                #cntrlLoc = {}
                cntrlLoc['expName'] = self.cntrlExpName
                cntrlLoc['diagName'] = myLoc['diagName'][0]
                normLoc = deepcopy(cntrlLoc)

                #Setting to avoid over-crowding
                if self.fcTDeltas.index(fcTDelta) > (self.MAX_FC_SUBFIGS-1): continue

                linesVals = defaultdict(list)
                linesLabel = []
                linesGroup = []
                for expName in self.expNames:
                    for diagnosticName in myLoc['diagName']:
                        if (expName == cntrlLoc['expName'] and
                            diagnosticName == cntrlLoc['diagName']): continue
                        cntrlLoc['cyDTime'] = myExpsCYDTimes[(expName, fcTDelta)]
                        linesGroup.append(expName)
                        linesLabel.append(self.expDiagnosticLabel(
                            expName, diagnosticName, myLoc['diagName']))

                        lineVals = defaultdict(list)

                        for binVal in binStrVals:
                            cntrlLoc['binVal'] = binVal
                            expLoc = deepcopy(cntrlLoc)
                            expLoc['diagName'] = diagnosticName
                            expLoc['expName'] = expName

                            X = tempdfw.loc(expLoc)
                            Y = tempdfw.loc(cntrlLoc)

                            # confidence intervals
                            ciVals = su.bootStrapClusterFunc(
                                         X, Y,
                                         n_samples = 10000,
                                         statNames = [statName])

                            # normalizing value for ratio
                            normLoc['binVal'] = binVal
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

                # define subplot title
                if useRelativeDifference:
                    title = varName
                else:
                    title = varLabel

                if len(self.fcTDeltas) > 1:
                  title += ' @ '+str(float(fcTDelta.total_seconds()) / 3600.0 / 24.0)+' days'


                # perform subplot agnostic plotting (all expNames)
                options['profilefunc'](
                    fig,
                    linesVals[su.cimean], binNumVals,
                    linesLabel,
                    title, binLabel, fcstatDiagLabel,
                    myBinConfigs['binConfig'],
                    sciTicks, logScale, centralValue,
                    nyplots, nxplots, nsubplots, iplot,
                    linesValsMinCI = linesVals[su.cimin],
                    linesValsMaxCI = linesVals[su.cimax],
                    dmin = dmin, dmax = dmax,
                    lineAttribOffset = 1,
                    interiorLabels = self.interiorLabels)

                iplot = iplot + 1

        # save figure
        filename = ('%s%s_BinValAxis_%s-%smin_%s_%s_%s'%(
                   myLoc['binVar'], self.binMethodFile(myLoc['binMethod']),
                   self.fcTDeltas_totmin[0], self.fcTDeltas_totmin[-1],
                   self.DiagSpaceName, fcDiagName, statName))

        pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.50, 0.55)
