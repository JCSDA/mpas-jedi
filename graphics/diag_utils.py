#!/usr/bin/env python3

import binning_utils as bu
from predefined_configs import anIter
from copy import deepcopy
import logging
import numpy as np
import pandas as pd
import plot_utils as pu
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

CRStatNames = ['CRd','CRh','CRo']

def ConsistencyRatioFromDFW(dfw, stateType):
# INPUTS:
# dfw - StatisticsDatabase.DFWrapper object containing
#       a pandas DataFrame with basic statistics for (at least)
#       'om'+stateType, 'sigma'+stateType, and 'sigmao'+stateType
# stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
    assert stateType in ['b', 'a', 'f'], 'ConsistencyRatioFromDFW: wrong stateType: '+stateType

#OUTPUT: pandas DataFrame containing only diagName == 'CRy'+stateType and CRStatNames data columns

    ommStr = 'om'+stateType
    sigmamStr = 'sigma'+stateType
    sigmaoStr = 'sigmao'+stateType

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
        #For each subpopulation, the consistency ratio can be written three different
        # ways, depending on the quantity of interest in the analysis. In the
        # following formulae, MS ≡ Mean Squared Value.

        #(1) CRd = SQRT( [MS(SIGMAMOD) + MS(SIGMAOBS)] / MS(OBS-MODEL) )
        # denominator zero check, q
        q = np.logical_and(bu.greatBound(np.absolute(omm), 0.0), p)
        CRStats[CRStatNames[0]][q] = (sigmam[q] + sigmao[q]) / omm[q]
        # benefit: ensures positive-definite numerator and denominator

        #(2) CRh = SQRT( MS(SIGMAMOD) / [MS(OBS-MODEL) - MS(SIGMAOBS)] )
        q = np.logical_and(bu.greatBound(np.absolute(omm - sigmao), 0.0), p)
        CRStats[CRStatNames[1]][q] = sigmam[q] / (omm[q] - sigmao[q])
        # benefit: directly measures the spread of the forecast in observation space

        #(3) CRo = SQRT( MS(SIGMAOBS) / [MS(OBS-MODEL) - MS(SIGMAMOD)] )
        q = np.logical_and(bu.greatBound(np.absolute(omm - sigmam), 0.0), p)
        CRStats[CRStatNames[2]][q] = sigmao[q] / (omm[q] - sigmam[q])
        # benefit: directly measures the specified ObsError


    # create an abbreviated dictionary with only the 'CRy'+stateType diagName and CR statistics
    CRyDict = dfw.loc({'diagName': [ommStr]}).reset_index().to_dict(orient='list')
    nrows = len(CRStats[CRStatNames[0]])
    CRyDict['diagName'] = ['CRy'+stateType] * nrows

    # perform the sqrt operation for all three CR's where positive-semi-definite
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
#        iter: iteration of "self" type variables ['bg', 'an', int, '', default: None=='0']
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
#        analysisStatistics (optional): statistics selected for this diagnostic (default: application dependent, see Analyses module)
#        label (optional): label for plot axes (TODO: hookup in Analyses.py)
#}
availableDiagnostics = {
    'bc': {
        'variable': BiasCorrection,
        'iter': 'an',
        'label': '$y_{bias}$',
    },
    'omb': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'bg',
        'label': '$y - x_b$',

    },
    'oma': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': '$y - x_a$',
    },
    'omf': {
        'variable': ObsMinusModel,
        'iter': '',
        'label': '$y - x_f$',
    },
    'mmgfsan': {
        'offline': True,
        'label': '$x - x_{a,GFS}$',
    },
    'omb_nobc': {
        'variable': ObsMinusModel,
        'iter': 'bg',
        'label': '$y - x_b$',
    },
    'oma_nobc': {
        'variable': ObsMinusModel,
        'iter': 'an',
        'label': '$y - x_a$',
    },
    'rltv_omb': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'bg',
        'label': '$y - x_b$',
    },
    'rltv_oma': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': '$y - x_a$',
    },
    'rltv_omf': {
        'variable': RelativeObsMinusModel,
        'iter': '',
        'label': '$y - x_f$',
    },
    'rltv_omb_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'bg',
        'label': '$y - x_b$',
    },
    'rltv_oma_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'an',
        'label': '$y - x_a$',
    },
    'obs': {
        'variable': vu.selfObsValue,
    },
    'obs_bc': {
        'variable': BiasCorrectedObs,
        'label': '$y$',
    },
    'h(x)': {
        'variable': vu.selfHofXValue,
        'iter': '',
        'label': '$h(x_f)$',
    },
    'bak': {
        'variable': vu.selfHofXValue,
        'iter': 'bg',
        'label': '$h(x_b)$',
    },
    'ana': {
        'variable': vu.selfHofXValue,
        'iter': 'an',
        'label': '$h(x_a)$',
    },
    'sigmaob': {
        'variable': vu.selfErrorValue,
        'iter': 'bg',
        'analyze': False,
        'label': '$\sigma_o$',
    },
    'sigmab': {
        'variable': STDofHofX,
        'iter': 'bg',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_b}$',
    },
    'sigmaoa': {
        'variable': vu.selfErrorValue,
        'iter': 'an',
        'analyze': False,
        'label': '$\sigma_o$',
    },
    'sigmaa': {
        'variable': STDofHofX,
        'iter': 'an',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_a}$',
    },
    'sigmaof': {
        'variable': vu.selfErrorValue,
        'iter': '',
        'analyze': False,
        'label': '$\sigma_o$',
    },
    'sigmaf': {
        'variable': STDofHofX,
        'iter': '',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_f}$',
    },
#TODO: replace 'derived' bool with 'derivedFrom' list.  E.g., for 'CRyb',
#      that would include 'sigmao' and 'sigmab'.  Update writediagstats_obsspace.
#TODO: then only add the derived diagnostic (e.g., CRyb) to conf.DiagSpaceConfig[:]['diagNames']
#      and remove the 'analyze' config entry under availableDiagnostics
    'CRyb': {
        'iter': 'bg',
        'DFWFunction': ConsistencyRatioFromDFW,
        'derived': True,
#        'derivedFrom': ['sigmao','sigmab']
        'staticArg': 'b',
        'analysisStatistics': CRStatNames,
        'label': '$CR_{y,b}$',
    },
    'CRya': {
        'iter': 'an',
        'DFWFunction': ConsistencyRatioFromDFW,
        'derived': True,
#        'derivedFrom': ['sigmao','sigmaa']
        'staticArg': 'a',
        'analysisStatistics': CRStatNames,
        'label': '$CR_{y,a}$',
    },
    'CRyf': {
        'iter': '',
        'DFWFunction': ConsistencyRatioFromDFW,
        'derived': True,
#        'derivedFrom': ['sigmaof','sigmaf']
        'staticArg': 'f',
        'analysisStatistics': CRStatNames,
        'label': '$CR_{y,f}$',
    },
    'SCI': {
        'variable': bu.SCIOkamoto,
        'iter': 'bg',
        'onlyObsSpaces': ['abi_g16'],
    },
}

#TODO: have this function return a list of diagnosticConfiguration or Diagnostic (new class) objects
#      instead of a list of dicts
def diagnosticConfigs(diagnosticNames_, ObsSpaceName, includeEnsembleDiagnostics=True,
                      analysisStatistics = su.allFileStats):

    diagnosticConfigs = {}
    diagnosticNames = list(set(diagnosticNames_))

    for diagnosticName in diagnosticNames:
        if diagnosticName in availableDiagnostics:
            config = deepcopy(availableDiagnostics[diagnosticName])
        elif diagnosticName in diagnosticConfigs:
            _logger.warning('diagnosticName is duplicated: '+diagnosticName)
            continue
        else:
            _logger.error('diagnosticName is undefined: '+diagnosticName)

        config['analyze'] = config.get('analyze', True)
        config['derived'] = config.get('derived', False)
        config['offline'] = config.get('offline', False)
        config[vu.mean] = config.get(vu.mean, True)
        config[vu.ensemble] = (config.get(vu.ensemble, False) and includeEnsembleDiagnostics)
        config['onlyObsSpaces'] = config.get('onlyObsSpaces',[])
        config['analysisStatistics'] = config.get('analysisStatistics', analysisStatistics)
        config['staticArg'] = config.get('staticArg', None)
        config['label'] = config.get('label',diagnosticName)

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
                outerIter = anIter
            elif pu.isint(outerIterStr) or outerIterStr == '':
                outerIter = outerIterStr
            else:
                _logger.error('outerIter is undefined: '+outerIterStr)
        diagnosticConfigs[diagnosticName]['outerIter'] = outerIter

        if config['derived']: continue
        if config['offline']: continue
        diagnosticConfigs[diagnosticName]['ObsFunction'] = bu.ObsFunctionWrapper(config)

    return diagnosticConfigs
