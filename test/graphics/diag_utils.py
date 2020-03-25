import binning_utils as bu
import numpy as np
import os
import var_utils as vu

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


## Table of available diagnostics
# format of each entry:
#    diagName: {
#        variable: variable name or function to calculate this diag
#        iter: iteration of "self" type variables ['bg','an',int]
#        osNames (optional): list of ObsSpaces for which this diag applies
#            see config.DiagSpaceConfig keys for available options
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
    'SCI': {
        'variable': bu.SCIOkamoto,
        'iter': 'bg',
        'osNames': ['abi_g16'],
    },
}

#groups of diagnostics for which statistics are calculated
#analysis outer iteration
#TODO: replace NOUTER with a command-line option or
#      move to central work-flow configuration script
NOUTER=os.getenv('NOUTER',0) #set equal to number of outer iterations

## diagnostics for which QC is irrelevant
nonQCedDiags = ['obs']

## difference diagnostics
diffDiagNames = ['omb']

## relative difference diagnostics
relDiagNames = ['rltv_omb']

## absolute diagnostics
absDiagNames = ['obs','bak']

cloudyRadDiagNames = ['SCI']

# analysis diagnostics
if NOUTER > 0:
    diffDiagNames.append('oma')
    relDiagNames.append('rltv_oma')
    absDiagNames.append('ana')
