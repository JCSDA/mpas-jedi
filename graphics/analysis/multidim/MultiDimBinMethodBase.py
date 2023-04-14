#!/usr/bin/env python3

import basic_plot_functions as bpf
import binning_utils as bu
from copy import deepcopy
import numpy as np
import var_utils as vu

from analysis.AnalysisBase import AnalysisBase
import analysis.StatisticsDatabase as sdb

class MultiDimBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with numerical binValues
      that are assigned numerical values, e.g., altitude, pressure, latitude, cloud fraction
    '''

    parallelism = True
    maxBinVarTier = 2

    cldFracTransform = 'logit'

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)
        # default 1D binVars
        self.binVarDict = {
            vu.obsVarAlt: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarImpact: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum altitude to show on all figures
                'maxvalue': 30000.,
              },
            },
            vu.obsVarACI: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarCldFracX: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
            vu.obsVarCldFracY: {'profilefunc': bpf.plotSeries, 'binVarTier': 1},
            vu.obsVarLat: {'profilefunc': bpf.plotProfile},
            vu.obsVarPrs: {'profilefunc': bpf.plotProfile},
            vu.obsVarCI: {'profilefunc': bpf.plotSeries, 'binVarTier': 2},
            vu.obsVarLogCI: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.modVarDiagPrs: {'profilefunc': bpf.plotProfile},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #vu.modVarLat: {'profilefunc': bpf.plotProfile, 'binVarTier': 1},
            vu.modVarLev: {
              'profilefunc': bpf.plotProfile,
              'binFilter': {
                # maximum model level to show on all figures
                #'maxvalue': 40,
              },
            },
            vu.obsVarGlint: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLandFrac: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarLT: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
            vu.obsVarSenZen: {'profilefunc': bpf.plotSeries, 'binVarTier': 3},
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
            selectedStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                selectedStatistics = set(list(selectedStatistics) +
                                         diagnosticConfigs[diagnosticName]['selectedStatistics'])
            availableStatistics = set([])
            for diagnosticName in diagnosticNames:
                diagnosticConfigs[diagnosticName] = deepcopy(self.diagnosticConfigs[diagnosticName])
                availableStatistics = set(list(availableStatistics) +
                                         diagnosticConfigs[diagnosticName]['availableStatistics'])
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagBinVars = self.db.dfw.levels('binVar', {'diagName': diagnosticNames})

            for fullBinVar, options in self.binVarDict.items():
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars): continue
                binVarLoc = {}
                binVarLoc['diagName'] = diagnosticNames
                binVarLoc['binVar'] = binVar
                binVarLoc['binVal'] = self.allBinNumVals2DasStr

                #Make figures for all binMethods
                binMethods = self.db.dfw.levels('binMethod', binVarLoc)
                for binMethod in binMethods:

                    #TODO: REMOVE, for testing only
                    #if binMethod != bu.identityBinMethod: continue

                    self.logger.info(diagnosticGroup+', '+binVar+', '+binMethod)

                    if useWorkers:
                        workers.apply_async(self.innerloopsWrapper,
                            args = (diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options))
                    else:
                        self.innerloopsWrapper(
                            diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, binVar, binMethod, selectedStatistics, options):

        myLoc = {}
        myLoc['binVar'] = binVar
        myLoc['binVal'] = self.allBinNumVals2DasStr
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
        binStrVals = mydfwDict['dfw'].levels('binVal')
        binUnits = mydfwDict['dfw'].uniquevals('binUnits')[0]

        # assume all bins represent same variable/units
        binLabel = binVar
        if binUnits != vu.miss_s:
            binLabel += ' ('+binUnits+')'
        for orig, sub in self.labelReplacements.items():
            binLabel = binLabel.replace(orig, sub)

        # bin info
        binNumVals = []
        for binVal in binStrVals:
            ibin = self.allBinStrVals.index(binVal)
            binNumVals.append(self.allBinNumVals[ibin])

        # filter out binVals less than (greater than) minvalue (maxvalue)
        binFilter = options.get('binFilter', None)
        # Example to include only values from 3 to 40:
        # 'binFilter': {
        #   'minvalue': 3,
        #   'maxvalue': 40,
        # }
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

        # special independent variable axis configs
        binVarIs = {}
        specialBinVars = [
          vu.obsVarPrs,
          vu.obsVarMCI,
          vu.obsVarOCI,
          vu.obsVarLogCI,
          vu.obsVarCldFracX,
          vu.obsVarCldFracY,
          vu.modVarDiagPrs,
        ]
        for var in specialBinVars:
            var_dict = vu.varDictAll.get(var,['',''])
            binVarIs[var] = (var_dict[1] == binVar)

        binConfig = deepcopy(bpf.defaultIndepConfig)
        pCoord = binVarIs[vu.obsVarPrs] or binVarIs[vu.modVarDiagPrs]
        binConfig['invert'] = pCoord
        if pCoord: binConfig['transform'] = 'Pressure'
        if binVarIs[vu.obsVarMCI] or binVarIs[vu.obsVarOCI] or binVarIs[vu.obsVarLogCI]:
            binConfig['transform'] = 'CloudImpact'
        if binVarIs[vu.obsVarCldFracX] or binVarIs[vu.obsVarCldFracY]:
            binConfig['transform'] = self.cldFracTransform

        # sort bins by numeric value
        indices = list(range(len(binNumVals)))
        indices.sort(key=binNumVals.__getitem__)
        binNumVals = list(map(binNumVals.__getitem__, indices))
        binStrVals = list(map(binStrVals.__getitem__, indices))

        myBinConfigs = {
            'str': binStrVals,
            'num': binNumVals,
            'binLabel': binLabel,
            'binConfig': binConfig,
        }
        if len(binStrVals) < 2: return

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

        for statName in selectedStatistics:
            if statName not in options.get('onlyStatNames', selectedStatistics): continue

            self.innerloops(
                mydfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options)

    def innerloops(self,
        dfwDict, diagnosticGroup, myLoc, statName, nVarsLoc, varMapLoc, myBinConfigs, options):
        '''
        virtual method
        '''
        raise NotImplementedError()
