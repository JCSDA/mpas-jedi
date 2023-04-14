#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
import numpy as np
import plot_utils as pu
import var_utils as vu

from analysis.AnalysisBase import AnalysisBase
import analysis.StatisticsDatabase as sdb

# TODO: generalize as a sub-class of MultiDimBinMethodBase
class BinValAxisStatsComposite(AnalysisBase):
    '''
    Similar to BinValAxisProfile, except
      all statistics (Count, Mean, RMS, STD) are placed on the same axis
      -  x-axis: binVal
      -    line: per statistic
      - subplot: per DiagSpace variable
      -    file: combination of FC lead time, experiment, and binMethod (if applicable)
    '''
    requestAggDFW = True

    maxBinVarTier = 2

    subplotWidth = 1.9
    subplotAspect = 0.9

    requiredStatistics = ['Count', 'Mean', 'RMS', 'STD']

    # Force serial processing so that console output is contiguous
    # TODO(JJG): output to an ascii file and remove this line
    blocking = True

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.binVarDict = {
            vu.obsVarAlt: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarImpact: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarACI: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarCI: {'statsfunc': bpf.plotfitRampComposite},
            vu.obsVarCldFracX: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 2},
            vu.obsVarCldFracY: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 1},
            vu.obsVarLat: {'statsfunc': bpf.plotCompositeProfile},
            #vu.obsVarLogCI: {'statsfunc': bpf.plotfitRampComposite},
            vu.obsVarPrs: {'statsfunc': bpf.plotCompositeProfile},
            vu.modVarDiagPrs: {'profilefunc': bpf.plotCompositeProfile},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #vu.modVarLat: {'statsfunc': bpf.plotCompositeProfile, 'binVarTier': 1},
            vu.modVarLev: {
              'statsfunc': bpf.plotCompositeProfile,
              'binFilter': {
                # maximum model level to show on all figures
                #'maxvalue': 40,
              },
            },
            vu.obsVarGlint: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarLandFrac: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarLT: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
            vu.obsVarSenZen: {'statsfunc': bpf.plotCompositeSeries, 'binVarTier': 3},
        }

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            selectedStatistics = self.requiredStatistics
            availableStatistics = diagnosticConfig['availableStatistics']
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticName})

            for fullBinVar, options in self.binVarDict.items():
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars): continue
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue

                binFilter = options.get('binFilter', None)

                myLoc = {}
                myLoc['diagName'] = diagnosticName
                myLoc['binVar'] = binVar
                myLoc['binVal'] = self.allBinNumVals2DasStr

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                ## Get all float/int binVals associated with binVar
                binMethods = mydfwDict['dfw'].levels('binMethod')
                binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

                # assume all bins represent same variable/units
                binLabel = binVar
                if binUnits != vu.miss_s:
                    binLabel += ' ('+binUnits+')'
                for orig, sub in self.labelReplacements.items():
                    binLabel = binLabel.replace(orig, sub)

                fcDiagName = self.fcName(diagnosticName)
                figPath = self.myFigPath/fcDiagName
                figPath.mkdir(parents=True, exist_ok=True)
                dataPath = figPath/'data'
                dataPath.mkdir(parents=True, exist_ok=True)

                # only analyze variables that have non-zero Count when sliced by myLoc
                nVarsLoc = 0
                varMapLoc = []
                for (varName, varLabel) in self.varMap:
                  if 'Count' in selectedStatistics:
                      countDF = mydfwDict['dfw'].loc({'varName': varName}, 'Count')
                      if countDF.shape[0] > 0:
                        if np.nansum(countDF.to_numpy()) > 0:
                          nVarsLoc += 1
                          varMapLoc.append((varName, varLabel))
                  else:
                      statDF = mydfwDict['dfw'].loc({'varName': varName}, list(selectedStatistics)[0])
                      if statDF.shape[0] > 0:
                        if np.isfinite(statDF.to_numpy()).sum() > 0:
                          nVarsLoc += 1
                          varMapLoc.append((varName, varLabel))

                nsubplots = nVarsLoc
                nxplots = np.int(np.ceil(np.sqrt(nsubplots)))
                while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
                nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

                ptLoc = {}

                #file loop 1
                for binMethod in binMethods:
                    ptLoc['binMethod'] = binMethod

                    methodDFW = sdb.DFWrapper.fromLoc(mydfwDict['dfw'], {'binMethod': binMethod})

                    # parse binVals, which may be different for each binMethod
                    binStrVals = methodDFW.levels('binVal')

                    # bin info
                    binNumVals = []
                    for binVal in binStrVals:
                        ibin = self.allBinStrVals.index(binVal)
                        binNumVals.append(self.allBinNumVals[ibin])

                    # sort bins by numeric value
                    indices = list(range(len(binNumVals)))
                    indices.sort(key=binNumVals.__getitem__)
                    binNumVals = list(map(binNumVals.__getitem__, indices))
                    binStrVals = list(map(binStrVals.__getitem__, indices))

                    if binFilter is not None:
                      remove = np.full_like(binNumVals, False, bool)
                      minvalue = binFilter.get('minvalue', None)
                      if minvalue is not None:
                        less = bu.lessBound(np.asarray(binNumVals), minvalue)
                        remove[less] = True
                      maxvalue = binFilter.get('maxvalue', None)
                      if maxvalue is not None:
                        great = bu.greatBound(np.asarray(binNumVals), maxvalue)
                        remove[great] = True

                      binStrVals = list(np.asarray(binStrVals)[~remove])
                      binNumVals = list(np.asarray(binNumVals)[~remove])

                    nBins = len(binStrVals)
                    if nBins < 2: continue

                    self.logger.info('binVar=>'+binVar+', binMethod=>'+binMethod)

                    methodDFWAgg = sdb.DFWrapper.fromAggStats(methodDFW, ['cyDTime'])

                    #file loop 2
                    for expName in self.expNames:
                        ptLoc['expName'] = expName

                        #file loop 3
                        for (fcTDelta, fcTDelta_totmin) in self.fcMap:
                            ptLoc['fcTDelta'] = fcTDelta

                            # establish a new figure
                            fig = pu.setup_fig(nxplots, nyplots, self.subplotWidth, self.subplotAspect, self.interiorLabels)
                            iplot = 0

                            ERRParams = {}
                            ERRParams[self.DiagSpaceName] = {}

                            #subplot loop 1
                            for (varName, varLabel) in varMapLoc:
                                ptLoc['varName'] = varName

                                #collect aggregated statNames, varying across fcTDelta
                                countsVals = np.full(nBins, 0)
                                meansVals  = np.full(nBins, np.NaN)
                                rmssVals   = np.full(nBins, np.NaN)
                                stdsVals   = np.full(nBins, np.NaN)

                                for ibin, binVal in enumerate(binStrVals):
                                    ptLoc['binVal'] = binVal
                                    c = methodDFWAgg.loc1(ptLoc,'Count')
                                    if np.isfinite(c):
                                      countsVals[ibin] = c
                                    meansVals[ibin] = methodDFWAgg.loc1(ptLoc,'Mean')
                                    rmssVals[ibin] = methodDFWAgg.loc1(ptLoc,'RMS')
                                    stdsVals[ibin] = methodDFWAgg.loc1(ptLoc,'STD')

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
                                    binLabel,
                                    nyplots, nxplots, nsubplots, iplot,
                                    interiorLabels = self.interiorLabels)

                                paramKey = self.chlist[iplot]
                                if paramKey < 0: paramKey = varName
                                if FitParams is not None:
                                    ERRParams[self.DiagSpaceName][(paramKey, binMethod)] = FitParams

                                iplot = iplot + 1

                            YAMLParams = {}
                            self.logger.info('#For binning_params('+expName+'):')
                            for key in sorted(ERRParams[self.DiagSpaceName]):
                                self.logger.info(binVar+"ErrParams['"+self.DiagSpaceName+"']["+str(key)+"]   = "+
                                       str(ERRParams[self.DiagSpaceName][key]['bin_utils']))
                                for param, val in ERRParams[self.DiagSpaceName][key]['YAML'].items():
                                    if param not in YAMLParams: YAMLParams[param] = []
                                    YAMLParams[param] += val
                            self.logger.info('#For UFO YAML config('+expName+'):')
                            for param, val in YAMLParams.items():
                                self.logger.info('#  '+param+': '+str(val))

                            # save each figure
                            filename = ('%s%s_BinValAxis_%smin_%s_%s_%s'%(
                                       binVar, self.binMethodFile(binMethod), fcTDelta_totmin,
                                       self.DiagSpaceName, fcDiagName, expName))

                            pu.finalize_fig(fig, str(figPath/filename), self.figureFileType, self.interiorLabels, 0.6)

                    # end expName loop

                # end binMethod loop

            # end fullBinVar loop
