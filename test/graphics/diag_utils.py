#!/usr/bin/env python3

import binning_utils as bu
from copy import deepcopy
import logging
import numpy as np
import os
import pandas as pd
import stat_utils as su
import var_utils as vu

_logger = logging.getLogger(__name__)

#====================================
# diagnostic functions and dictionary
#====================================

class BiasCorrectedObs:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)

    def evaluate(self, dbVals, caseParams):
        bcdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        mod = dbVals[caseParams['base2db'][vu.selfHofXValue]]
        return np.subtract(mod,bcdep)


class BiasCorrection:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, caseParams):
        bcdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        mod = dbVals[caseParams['base2db'][vu.selfHofXValue]]
        obs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        return np.subtract(np.subtract(mod,bcdep),obs)


class BiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)

    def evaluate(self, dbVals, caseParams):
        bcdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        return np.negative(bcdep)


class RelativeBiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, caseParams):
        bcdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        obs = dbVals[caseParams['base2db'][vu.selfObsValue]]

        OMM = np.negative(bcdep)
        valid = bu.greatBound(np.abs(obs), 0.0, False)
        OMM[valid] = np.divide(OMM[valid],obs[valid])
        OMM[~valid] = np.NaN

        return np.multiply(OMM,100.0)


class ObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, caseParams):
        obs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        mod = dbVals[caseParams['base2db'][vu.selfHofXValue]]
        return np.subtract(obs,mod)


class RelativeObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, caseParams):
        obs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        mod = dbVals[caseParams['base2db'][vu.selfHofXValue]]

        OMM = np.subtract(obs,mod)
        valid = bu.greatBound(np.abs(obs), 0.0, False)
        OMM[valid] = np.divide(OMM[valid],obs[valid])
        OMM[~valid] = np.NaN

        return np.multiply(OMM,100.0)

class STDofHofX:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)

    def evaluate(self, dbVals, caseParams):
        meanVarName = caseParams['base2db'][vu.selfHofXValue]
        memberKeys = []
        for key in dbVals.keys():
            if meanVarName+vu.ensSuffixBase in key:
                memberKeys.append(key)
        nMembers = len(memberKeys)
        if nMembers > 0:
            nLocs = len(dbVals[memberKeys[0]])
            mods = np.full((nMembers, nLocs), np.NaN)
            for member, key in enumerate(memberKeys):
                mods[member,:] = dbVals[key]
            std = np.nanstd(mods, axis=0, ddof=1)
        else:
            std = np.full(nLocs, np.NaN)
        return std

CRStatNames = ['d:(Sh+So)','(d-So):Sh','(d-Sh):So']

def ConsistencyRatioFromDFW(dfw, stateType):
# INPUTS:
# dfw - StatisticsDatabase.DFWrapper object containing
#       a pandas DataFrame with basic statistics for (at least)
#       'om'+stateType, 'sigma'+stateType, and 'sigmao'
# stateType - model state diagnostic type (either 'b' for background or 'a' for analysis)

#OUTPUT: pandas DataFrame containing only 'CRy'+stateType diagName and CRStatNames statistic columns

#For each subpopulation, the consistency ratio can be written differently,
# depending on the quantity of interest in the analysis.  For example:
#(1)
#    CR = SQRT( MS(OBS-MODEL) / [ MS(SIGMAMOD) + MS(SIGMAOBS) ] )
#
# benefit: ensures positive-definite numerator and denominator
#
#(2) Minamide and Zhang, 2019
#    CR = SQRT( [MS(OBS-MODEL) - MS(SIGMAOBS) ] / MS(SIGMAMOD) )
#
# benefit: directly measures the spread of the forecast in observation space
#
#(3)
#    CR = SQRT( [MS(OBS-MODEL) - MS(SIGMAMOD) ] / MS(SIGMAOBS) )
#
# benefit: directly measures the specified ObsError
#
# where MS ≡ Mean Squared Value

    import sys
    np.set_printoptions(threshold=sys.maxsize)

    ommStr = 'om'+stateType
    sigmamStr = 'sigma'+stateType
    sigmaoStr = 'sigmao'

    # get the three variables needed as numpy arrays
    availableDiagNames = dfw.levels('diagName')
    if (ommStr not in availableDiagNames or
        sigmamStr not in availableDiagNames or
        sigmaoStr not in availableDiagNames): return None

    omm = dfw.loc({'diagName': ommStr}, 'MS').to_numpy()
    sigmam = dfw.loc({'diagName': sigmamStr}, 'MS').to_numpy()
    sigmao = dfw.loc({'diagName': sigmaoStr}, 'MS').to_numpy()

    CRStats = {}
    for statName in CRStatNames:
        CRStats[statName] = np.full_like(omm, np.NaN)

    p = np.logical_and(np.isfinite(omm),
          np.logical_and(np.isfinite(sigmam), np.isfinite(sigmao)))
    if p.sum() > 0:
        # perform math operations
        #(1)
        q = np.logical_and(bu.greatBound(np.absolute(sigmam + sigmao), 0.0), p)
        CRStats['d:(Sh+So)'][q] = omm[q] / (sigmam[q] + sigmao[q])
        #(2)
        q = np.logical_and(bu.greatBound(np.absolute(sigmam), 0.0), p)
        CRStats['(d-So):Sh'][q] = (omm[q] - sigmao[q]) / sigmam[q]
        #(3)
        q = np.logical_and(bu.greatBound(np.absolute(sigmao), 0.0), p)
        CRStats['(d-Sh):So'][q] = (omm[q] - sigmam[q]) / sigmao[q]

    # create an abbreviated dictionary with only the new diagName and CR statistic
    CRyDict = dfw.loc({'diagName': [ommStr]}).reset_index().to_dict(orient='list')
    nrows = len(CRStats[CRStatNames[0]])
    CRyDict['diagName'] = ['CRy'+stateType] * nrows
    for statistic, CRStat in CRStats.items():
        CRStat[bu.lessBound(CRStat, 0.0)] = np.NaN
        CRyDict[statistic] = np.sqrt(CRStat)

    for statName in su.allFileStats:
      del CRyDict[statName]

    # convert dict to DataFrame
    CRyDF = pd.DataFrame.from_dict(CRyDict)
    del CRyDict
    CRyDF.set_index(dfw.indexNames, inplace=True)
    CRyDF.sort_index(inplace=True)

    return CRyDF

## Table of available diagnostics
# format of each entry:
#    diagName: {
#        variable: variable name or class used to calculate this diag
#        iter: iteration of "self" type variables ['bg','an',int]
#        vu.mean (optional): whether to apply this diagnostic to the mean state (default: True)
#        vu.ensemble (optional): whether this diagnostic requires variables from ensemble IODA files (default: False)
#        onlyObsSpaces (optional): list of ObsSpaces for which this diag applies
#            see config.DiagSpaceConfig keys for available options
#        analyze (optional): whether to analyze this diagnostic (defualt: True).  Useful for turning off analyses for diagnostics that are only needed for calculating other derived diagnostics.
#        derived (optional): whether the diagnostic is derived from statistics of other diagnostics
#                 if True, then will be calculated in analysis step, not statistics calculation step
#                 (default: False)
#        function (optional): if derived, the python function used to calculate this diagnostic
#                             there will be 2 arguments to this function (1) a StatisticsDatabase.DFWrapper object and (2) 'staticArg' below
#        staticArg (optional): a static argument to function  (default: None)
#        analysisStatistics (optional): statistics selected for this diagnostic (default: application dependent, see AnalyzeStatistics module)
#}
availableDiagnostics = {
    'bc': {
        'variable': BiasCorrection,
        'iter': 'an',
    },
    'omb': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'bg',
    },
    'oma': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'an',
    },
    'omb_nobc': {
        'variable': ObsMinusModel,
        'iter': 'bg',
    },
    'oma_nobc': {
        'variable': ObsMinusModel,
        'iter': 'an',
    },
    'rltv_omb': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'bg',
    },
    'rltv_oma': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'an',
    },
    'rltv_omb_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'bg',
    },
    'rltv_oma_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'an',
    },
    'obs': {
        'variable': vu.selfObsValue,
    },
    'obs_bc': {
        'variable': BiasCorrectedObs,
    },
    'bak': {
        'variable': vu.selfHofXValue,
        'iter': 'bg',
    },
    'ana': {
        'variable': vu.selfHofXValue,
        'iter': 'an',
    },
    'sigmao': {
        'variable': vu.selfErrorValue,
        'iter': 'bg',
        'analyze': False,
    },
    'sigmab': {
        'variable': STDofHofX,
        'iter': 'bg',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
    },
    'sigmaa': {
        'variable': STDofHofX,
        'iter': 'an',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
    },
    'CRyb': {
        'iter': 'bg',
        'DFWFunction': ConsistencyRatioFromDFW,
        'derived': True,
        'staticArg': 'b',
        'analysisStatistics': CRStatNames,
    },
    'CRya': {
        'iter': 'an',
        'DFWFunction': ConsistencyRatioFromDFW,
        'derived': True,
        'staticArg': 'a',
        'analysisStatistics': CRStatNames,
    },
    'SCI': {
        'variable': bu.SCIOkamoto,
        'iter': 'bg',
        'onlyObsSpaces': ['abi_g16'],
    },
}

## groups of diagnostics for which statistics are calculated

## diagnostics for which QC is irrelevant
nonQCedDiags = ['obs']

## difference diagnostics
diffDiagNames = ['omb']

## relative difference diagnostics
relDiagNames = ['rltv_omb']

## absolute diagnostics
absDiagNames = ['obs','bak']

cloudyRadDiagNames = ['SCI']

## STD diagnostics
sigmaDiagNames = ['sigmao','sigmab','CRyb']

# analysis diagnostics
if vu.nOuter > 0:
    diffDiagNames.append('oma')
    relDiagNames.append('rltv_oma')
    absDiagNames.append('ana')
    sigmaDiagNames += ['sigmaa','CRya']

#TODO: have this function return a list of diagnosticConfiguration or Diagnostic (new class) objects
#      instead of a list of dicts
def diagnosticConfigs(diagnosticNames, ObsSpaceName, nMembers=0, analysisStatistics = su.allFileStats):

    diagnosticConfigs = {}

    for diagnosticName in diagnosticNames:
        if diagnosticName in availableDiagnostics:
            config = deepcopy(availableDiagnostics[diagnosticName])
        else:
            _logger.error('diagnosticName is undefined: '+diagnosticName)

        config['analyze'] = config.get('analyze', True)
        config['derived'] = config.get('derived', False)
        config[vu.mean] = config.get(vu.mean, True)
        config[vu.ensemble] = (config.get(vu.ensemble, False) and nMembers > 1)
        config['onlyObsSpaces'] = config.get('onlyObsSpaces',[])
        config['analysisStatistics'] = config.get('analysisStatistics', analysisStatistics)
        config['staticArg'] = config.get('staticArg', None)

        # diagnosticConfig is undefined for the following cases
        if (not config[vu.mean] and not config[vu.ensemble]): continue
        if (len(config['onlyObsSpaces']) > 0 and
            ObsSpaceName not in config['onlyObsSpaces']): continue

        config['osName'] = ObsSpaceName

        # add this diagnosticConfig to the list
        diagnosticConfigs[diagnosticName] = deepcopy(config)

        outerIter = '0'
        outerIterStr = config.get('iter',None)
        if outerIterStr is not None:
            if outerIterStr == 'bg':
                outerIter = vu.bgIter
            elif outerIterStr == 'an':
                outerIter = vu.anIter
            elif pu.isint(outerIterStr):
                outerIter = outerIterStr
            else:
                _logger.error('outerIter is undefined: '+outerIterStr)
        diagnosticConfigs[diagnosticName]['outerIter'] = outerIter

        if config['derived']: continue
        diagnosticConfigs[diagnosticName]['ObsFunction'] = bu.ObsFunctionWrapper(config)

    return diagnosticConfigs
