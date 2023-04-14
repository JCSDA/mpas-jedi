#!/usr/bin/env python3

import binning_utils as bu
import predefined_configs as pconf
from collections.abc import Iterable
from copy import deepcopy
import numpy as np
import plot_utils as pu
import analysis.StatisticsDatabase as sdb
import var_utils as vu

from analysis.AnalysisBase import AnalysisBase

#=============
# 1-D figures
#=============
class CategoryBinMethodBase(AnalysisBase):
    '''
    Base class used to analyze statistics across binMethods with zero-dimensioned or
      category binValues, e.g., QC flag, named latitude band, cloudiness regime, surface type
    '''

    parallelism = True
    maxBinVarTier = 2

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        # default binVar/binMethod combinations
        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarLat, bu.latbandsMethod): {},
            (vu.obsVarLat, bu.troplatbandsMethod): {},
            (vu.obsVarLat, bu.polarlatbandsMethod): {},
            (vu.obsVarCldFracX, bu.cloudbandsMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
            # vu.modVarLat is redundant with vu.obsVarLat (both have varShort=="lat")
            #(vu.modVarLat, bu.latbandsMethod): {},
            (vu.noBinVar, bu.noBinMethod): {},
            (vu.obsRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.modelRegionBinVar, bu.geoirlatlonboxMethod): {'binVarTier': 2},
            (vu.obsVarPrs, bu.PjetMethod): {'binVarTier': 3},
            (vu.obsVarAlt, bu.altjetMethod): {'binVarTier': 3},
            (vu.obsVarImpact, bu.altjetMethod): {'binVarTier': 3},
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
            while nsubplots%nxplots > 0 and nsubplots%nxplots / nxplots <= 0.5: nxplots += 1
            nyplots = np.int(np.ceil(np.true_divide(nsubplots, nxplots)))

        return nxplots, nyplots, nsubplots

    def analyze_(self, workers = None):
        useWorkers = (not self.blocking and self.parallelism and workers is not None)

        # TODO(JJG): construct member Diagnostic objects (create new class) from
        #            diagnosticConfigs instead of referencing dictionary
        #            entries below.
        # TODO(JJG): use same color, vary line style within diagnosticGroupings
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

            diagLoc = {'diagName': diagnosticNames}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)
            for (fullBinVar, binMethod), options in self.binVarDict.items():
                if options.get('binVarTier', 1) > self.maxBinVarTier: continue
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars or
                    (binMethod is not None and binMethod not in diagBinMethods)): continue

                self.logger.info(diagnosticGroup+', '+binVar+', '+str(binMethod))

                if useWorkers:
                    workers.apply_async(self.innerloopsWrapper,
                        args = (diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options))
                else:
                    self.innerloopsWrapper(
                        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options)

    def innerloopsWrapper(self,
        diagnosticGroup, diagnosticConfigs, fullBinVar, binMethod, selectedStatistics, options):

        binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]

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

        for statName in selectedStatistics:
            if statName not in options.get('onlyStatNames', selectedStatistics): continue

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


def categoryBinValsAttributes(dfw, fullBinVar, binMethod, options):
    '''
    Utility function for providing an ordered list of
    pairs of binVals and associated labels for
    category binMethods in the context of a DFWrapper
    '''

    binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]

    dbSelect1DBinVals = dfw.levels('binVal')
    binUnitss = dfw.uniquevals('binUnits')
    #if (len(binUnitss) == 0 or
    #    len(dbSelect1DBinVals) == 1): return None, None
    assert (len(binUnitss) != 0 and len(dbSelect1DBinVals) > 0), 'ERROR: categoryBinValsAttributes received invalid binVar/binMethod'

    binUnits = binUnitss[0]

    # reorder select1DBinVals to match binMethod definition
    # TODO(JJG): clean up for readability
    select1DBinVals = []
    binMethods = pconf.binVarConfigs.get(fullBinVar, {})
    for method, value in binMethods.items():
      if binMethod is None or method == binMethod:
        tmp = value.get('values', dbSelect1DBinVals)
        if (not isinstance(tmp, Iterable) or
          isinstance(tmp, str)):
          select1DBinVals += [tmp]
        else:
          select1DBinVals += list(tmp)

    for b in dbSelect1DBinVals:
        select1DBinVals.append(b)
    for b in list(select1DBinVals):
        if b not in dbSelect1DBinVals:
            select1DBinVals.remove(b)

    select1DBinVals = pu.uniqueMembers(select1DBinVals)

    binTitles = []
    for binVal in select1DBinVals:
        t = ''
        if pu.isfloat(binVal) or pu.isint(binVal):
            t += ' @ '+binVar+'='+binVal
            if binUnits != vu.miss_s:
                t += ' '+binUnits
        elif binVal not in pconf.goodFlagNames+[bu.allQCMethod]:
            t += ' @ '+binVal
        binTitles.append(t)

    binValsMap = list(zip(select1DBinVals, binTitles))

    return binValsMap


