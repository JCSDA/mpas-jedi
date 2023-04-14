#!/usr/bin/env python3

import binning_utils as bu
import numpy as np
import var_utils as vu

from analysis.AnalysisBase import AnalysisBase
import analysis.StatisticsDatabase as sdb

class GrossValues(AnalysisBase):
    '''
    Calculate gross statistics for specified category binMethods at first forecast length
      NOTE: currently only calculates statistics at self.fcTDeltas[0]
            adjust minimum forecast length in order to calculate
            for non-zero forecast lengths, assuming those lengths
            are present in db
    '''

    requestAggDFW = True

    # Force serial processing so that console output is contiguous
    # TODO(JJG): output to an ascii file and remove this line
    blocking = True

    def __init__(self, db:sdb, analysisType:str, diagnosticGroupings:dict):
        super().__init__(db, analysisType, diagnosticGroupings)

        self.binVarDict = {
            (vu.obsVarQC, bu.goodQCMethod): {},
            (vu.obsVarCldFracY, bu.cloudbandsMethod): {},
        }

    def analyze_(self, workers = None):
        for diagnosticName, diagnosticConfig in self.diagnosticConfigs.items():
            if diagnosticName not in self.db.dfw.levels('diagName'): continue
            selectedStatistics = diagnosticConfig['selectedStatistics']
            availableStatistics = diagnosticConfig['availableStatistics']
            if not set(self.requiredStatistics).issubset(availableStatistics): continue

            diagLoc = {'diagName': diagnosticName}
            diagBinVars = self.db.dfw.levels('binVar', diagLoc)
            diagBinMethods = self.db.dfw.levels('binMethod', diagLoc)

            for (fullBinVar, binMethod), options in self.binVarDict.items():
                binVar = vu.varDictAll.get(fullBinVar, [None, fullBinVar])[1]
                if (binVar not in diagBinVars or
                    binMethod not in diagBinMethods): continue

                # narrow mydfwDict by binVar and binMethod to reduce run-time and memory
                myLoc = {}
                myLoc['binVar'] = binVar
                myLoc['binMethod'] = binMethod

                # reducing to mydfwDict speeds extractions in innerloops
                mydfwDict = {'dfw': self.db.loc(myLoc)}

                if self.requestAggDFW:
                    mydfwDict['agg'] = sdb.DFWrapper.fromAggStats(mydfwDict['dfw'], ['cyDTime'])
                    sdb.createORreplaceDerivedDiagnostics(mydfwDict['agg'], {diagnosticName: diagnosticConfig})

                # further narrow mydfwDict by diagName
                # NOTE: derived diagnostics may require multiple diagName values;
                # can only narrow by diagName after aggregation
                myLoc['diagName'] = diagnosticName
                for key in mydfwDict.keys():
                    mydfwDict[key] = sdb.DFWrapper.fromLoc(mydfwDict[key], myLoc)

                print(' Calculate gross statistics: binVar=>'+binVar+', binMethod=>'+binMethod)

                binValsMap = categoryBinValsAttributes(
                    mydfwDict['dfw'], fullBinVar, binMethod, options)

                print(' at FC length ', self.fcTDeltas[0])
                # Calculate gross statistics for this binVal
                statsLoc = {}
                statsLoc['fcTDelta'] = self.fcTDeltas[0]
                for binVal, binTitle in binValsMap:
                    statsLoc['binVal'] = binVal
                    GrossValues = {}
                    for varName in self.varNames:
                        statsLoc['varName'] = varName
                        for expName in self.expNames:
                            statsLoc['expName'] = expName
                            statsDFW = sdb.DFWrapper.fromLoc(mydfwDict['agg'], statsLoc)
                            for statName in selectedStatistics:
                                GrossValues[(statName, expName, varName)] = statsDFW.var(statName).to_numpy()
                    for expName in self.expNames:
                        print('Gross statistics for')
                        print('experiment=>'+expName)
                        if len(binValsMap) > 1:
                            print('binVal=>'+binVal)
                        print(' variables: ', self.varNames)

                        for statName in selectedStatistics:
                            print(statName)
                            tmp = np.asarray([])
                            for varName in self.varNames:
                                tmp = np.append(tmp, GrossValues[(statName, expName, varName)])
                            print(tmp)
