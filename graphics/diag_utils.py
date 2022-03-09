#!/usr/bin/env python3

import binning_utils as bu
from predefined_configs import anIter, appName
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
def negativeSafeSqrt(q):
    q_ = deepcopy(q)
    q_[bu.lessBound(q, 0.0)] = np.NaN
    q_ = np.sqrt(q_)

    return q_


class BiasCorrectedObs:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)

    def evaluate(self, dbVals, insituParameters):
        bcdep = dbVals[insituParameters[vu.selfDepValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]
        return np.subtract(mod, bcdep)


class BiasCorrection:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        bcdep = dbVals[insituParameters[vu.selfDepValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]
        return np.subtract(np.subtract(mod, bcdep), obs)


class BiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)

    def evaluate(self, dbVals, insituParameters):
        bcdep = dbVals[insituParameters[vu.selfDepValue]]
        return np.negative(bcdep)


class RelativeBiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        bcdep = dbVals[insituParameters[vu.selfDepValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]

        OMM = np.negative(bcdep)
        valid = bu.greatBound(np.abs(obs), 0.0, False)
        OMM[valid] = np.divide(OMM[valid], obs[valid])
        OMM[~valid] = np.NaN

        return np.multiply(OMM, 100.0)


class ObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        obs = dbVals[insituParameters[vu.selfObsValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]
        return np.subtract(obs, mod)


class RelativeObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        obs = dbVals[insituParameters[vu.selfObsValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]

        OMM = np.subtract(obs, mod)
        valid = bu.greatBound(np.abs(obs), 0.0, False)
        OMM[valid] = np.divide(OMM[valid],obs[valid])
        OMM[~valid] = np.NaN

        return np.multiply(OMM, 100.0)


class AnalysisMinusBackground:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.selfHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        return np.subtract(ana, bak)


class RelativeAnalysisMinusBackground:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.selfHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]

        AMB = np.subtract(ana, bak)
        valid = bu.greatBound(np.abs(obs), 0.0, False)
        AMB[valid] = np.divide(AMB[valid], obs[valid])
        AMB[~valid] = np.NaN

        return np.multiply(AMB, 100.0)


class AnalysisMinusBackgroundOverObsMinusBackground:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.selfHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]

        AMBoOMB = np.subtract(ana, bak)
        OMB = np.subtract(obs, bak)
        valid = bu.greatBound(np.abs(OMB), 0.0, False)

        AMBoOMB[valid] = np.divide(AMBoOMB[valid], OMB[valid])
        AMBoOMB[~valid] = np.NaN

        return AMBoOMB


class DerivedDiagnostic:
    '''
    base class for derived combinations of diagnostics and statistics
    The main functionality is located in the evaluate method, which is specified
    in the derived class.
    '''
    def __init__(self):
        self.availableStatistics = su.allFileStats
        self.diagname = 'DerivedDiagnostic'
        self.requiredDiagnostics = []

    @staticmethod
    def retrieveDiagnosticStat(dfw, diagnostic, statistic):
        return dfw.loc({'diagName': diagnostic}, statistic).to_numpy()

    def evaluate(self, dfw):
        '''
        virtual method
        ARGS
            dfw - StatisticsDatabase.DFWrapper object containing
                  a pandas DataFrame with basic statistics for (at least) self.requiredDiagnostics
        RETURN
            pandas DataFrame containing only diagName == self.diagname and self.availableStatistics data columns
            that has equivalent categorical information to the dfw slice
        '''
        raise NotImplementedError()

    def templateDFfromStats(self, dfw, templateDiagName, Stats):
        '''
        helper function to convert a dict of statistics to a pandas DataFrame
        ARGS
            dfw - source StatisticsDatabase.DFWrapper object
            templateDiagName - diagName on which to template the pandas DataFrame
            Stats - abbreviated dict containing key-value pairs of statistic names and np.array objects, e.g.:
                    Stats = {
                        availableStatistics[0]: np.array,
                        availableStatistics[1]: np.array,
                        etc...
                    }

        RETURN
            pandas.DataFrame templated by slicing dfw with diagName=templateDiagName
            contains self.diagname as diagName and Stats key-value pairs as data columns
        '''

        outDict = dfw.loc({'diagName': [templateDiagName]}).reset_index().to_dict(orient='list')
        nrows = len(Stats[self.availableStatistics[0]])
        outDict['diagName'] = [self.diagname] * nrows

        # transfer the statistics to outDict
        for key, value in Stats.items():
            outDict[key] = deepcopy(value)

        for statName in su.allFileStats:
          if statName not in Stats:
              del outDict[statName]

        # convert dict to DataFrame
        outDF = pd.DataFrame.from_dict(outDict)
        del outDict
        outDF.set_index(dfw.indexNames, inplace=True)
        outDF.sort_index(inplace=True)

        return outDF


class idealSigmao(DerivedDiagnostic):
    '''
    Calculates ideal sigmao for the
    background, analysis, or forecast states
    '''
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'idealSigmao: wrong stateType => '+stateType
        self.diagname = 'ideal-sigmao'+stateType
        self.availableStatistics = ['MS', 'RMS']

        self.label = '$ideal-sigmao_{'+stateType+'}$',
        self.omm = 'om'+stateType
        self.sigmam = 'sigma'+stateType

        self.requiredDiagnostics = [self.omm, self.sigmam]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        omm = self.retrieveDiagnosticStat(dfw, self.omm, 'MS')
        sigmam = self.retrieveDiagnosticStat(dfw, self.sigmam, 'MS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(omm, np.NaN)

        p = np.logical_and(np.isfinite(omm), np.isfinite(sigmam))
        if p.sum() > 0:
            # idealSigmao(MS) = MS(OBS-MODEL) - MS(SIGMAMOD)
            Stats['MS'][p] = omm[p] - sigmam[p]

            # idealSigmao(RMS) = SQRT( [MS(OBS-MODEL) - MS(SIGMAMOD)] )
            Stats['RMS'][p] = omm[p] - sigmam[p]
            Stats['RMS'] = negativeSafeSqrt(Stats['RMS'])

        return self.templateDFfromStats(dfw, self.omm, Stats)


class ObsSpaceConsistencyRatio(DerivedDiagnostic):
    '''
    Calculates multiple forms of the observation-space ensemble consistency ratio for the
    background, analysis, or forecast states
    '''
    CRd = 'CRd'
    CRo = 'CRo'
    CRh = 'CRh'
    invCRo = 'invCRo'
    invCRh = 'invCRh'
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'ObsSpaceConsistencyRatio: wrong stateType => '+stateType
        self.diagname = 'CRy'+stateType
        self.availableStatistics = [self.CRd, self.CRo, self.CRh, self.invCRo, self.invCRh]

        self.label = '$CR_{y,'+stateType+'}$',
        self.omm = 'om'+stateType
        self.sigmao = 'sigmao'+stateType
        self.sigmam = 'sigma'+stateType

        self.requiredDiagnostics = [self.omm, self.sigmao, self.sigmam]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        omm = self.retrieveDiagnosticStat(dfw, self.omm, 'MS')
        sigmam = self.retrieveDiagnosticStat(dfw, self.sigmam, 'MS')
        sigmao = self.retrieveDiagnosticStat(dfw, self.sigmao, 'MS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(omm, np.NaN)

        p = np.logical_and(np.isfinite(omm),
              np.logical_and(np.isfinite(sigmam), np.isfinite(sigmao)))
        if p.sum() > 0:
            #For each subpopulation, the consistency ratio can be written three different
            # ways, depending on the quantity of interest in the analysis. In the
            # following formulae, MS ≡ Mean Squared Value.

            #(1) CRd = SQRT( [MS(SIGMAMOD) + MS(SIGMAOBS)] / MS(OBS-MODEL) )
            # denominator zero check, q
            q = np.logical_and(bu.greatBound(np.absolute(omm), 0.0), p)
            Stats[self.CRd][q] = (sigmam[q] + sigmao[q]) / omm[q]
            # utility: ensures positive-definite numerator and denominator

            #(2) CRo = SQRT( MS(SIGMAOBS) / [MS(OBS-MODEL) - MS(SIGMAMOD)] )
            q = np.logical_and(bu.greatBound(np.absolute(omm - sigmam), 0.0), p)
            Stats[self.CRo][q] = sigmao[q] / (omm[q] - sigmam[q])
            # utility: directly measures the specified ObsError

            #(3) CRh = SQRT( MS(SIGMAMOD) / [MS(OBS-MODEL) - MS(SIGMAOBS)] )
            q = np.logical_and(bu.greatBound(np.absolute(omm - sigmao), 0.0), p)
            Stats[self.CRh][q] = sigmam[q] / (omm[q] - sigmao[q])
            # utility: directly measures the spread of the forecast in observation space

            # Additional diagnostic quantities useful for tuning spread in real-data experiments
            #(4) invCRo = SQRT( [MS(OBS-MODEL) - MS(SIGMAMOD)] / MS(SIGMAOBS) )
            q = np.logical_and(bu.greatBound(np.absolute(sigmao), 0.0), p)
            Stats[self.invCRo][q] = (omm[q] - sigmam[q]) / sigmao[q]
            # utility: gives scaling factors for adjusting the specified ObsError

            #(5) invCRh = SQRT( [MS(OBS-MODEL) - MS(SIGMAOBS)] / MS(SIGMAMOD) )
            q = np.logical_and(bu.greatBound(np.absolute(sigmam), 0.0), p)
            Stats[self.invCRh][q] = (omm[q] - sigmao[q]) / sigmam[q]
            # utility: gives scaling factors for adjusting the spread of the forecast in observation
            #          space.  May not translate directly to model-space inflation due to
            #          obs operator nonlinearities.

            # perform the sqrt operation for all CR's where positive-semi-definite
            for statistic, Stat in Stats.items():
                Stats[statistic] = negativeSafeSqrt(Stat)

        return self.templateDFfromStats(dfw, self.omm, Stats)

class SumDiagnosticsStatisticsPairs(DerivedDiagnostic):
    '''
    As a base class, SumDiagnosticsStatisticsPairs enables the renaming and summing of N non-derived
    diagnostics and statistic pairs.  The primary purpose is to plot multiple unrelated statistics
    on the same figure by giving them identical statistic names.
    '''
    statName = 'NA'
    diagName = 'NA'
    sourceStatistics = ['NA']
    sourceDiagnosticPrefixes = ['NA']
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'SumDiagnosticsStatisticsPairs: wrong stateType => '+stateType
        self.availableStatistics = [self.statName]
        self.diagname = self.diagName+'_'+stateType
        self.label = '$'+self.diagName+'$',
        self.requiredDiagnostics = [d+stateType for d in self.sourceDiagnosticPrefixes]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        diagStats = []
        for d, s in list(zip(self.requiredDiagnostics, self.sourceStatistics)):
            diagStats.append(self.retrieveDiagnosticStat(dfw, d, s))

        # re-store them in a dictionary with the new statName
        Stats = {}
        Stats[self.statName] = np.full_like(diagStats[0], np.NaN)

        p = np.full_like(diagStats[0], True, bool)
        for ds in diagStats:
            p = np.logical_and(p, np.isfinite(ds))

        if p.sum() > 0:
            Stats[self.statName][p] = diagStats[0][p]
            for ds in diagStats[1:]:
                Stats[self.statName][p] += ds[p]

            # perform the postOperator when provided
            if self.postOperator is not None:
                Stats[self.statName] = self.postOperator(Stats[self.statName])

        # convrt to a DataFrame
        return self.templateDFfromStats(dfw, self.requiredDiagnostics[0], Stats)

    @staticmethod
    def postOperator(q):
        return q


class DepartureSTD(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the LHS term (with mean removed) of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'RMSq'
    diagName = vu.DiagnosticVars[vu.EddT]
    sourceStatistics = ['STD']
    sourceDiagnosticPrefixes = ['om']


class DepartureRMS(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the LHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'RMSq'
    diagName = vu.DiagnosticVars[vu.EddT]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['om']

class EnsembleSpread(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the 1st RHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'RMSq'
    diagName = vu.DiagnosticVars[vu.HBHT]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['sigma']


class ObsError(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the 2nd RHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'RMSq'
    diagName = vu.DiagnosticVars[vu.R]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['sigmao']


class TotalSpread(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the full RHS of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'RMSq'
    diagName = vu.DiagnosticVars[vu.HBHTplusR]
    sourceStatistics = ['MS', 'MS']
    sourceDiagnosticPrefixes = ['sigma', 'sigmao']

    @staticmethod
    def postOperator(q):
        return negativeSafeSqrt(q)


class ObsErrorNormalizedInnovation(DerivedDiagnostic):
    '''
    Calculates the ratio of [y - h(x)]/ObsError for the
    background, analysis, or forecast states
    '''
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'ObsErrorNormalizedInnovation: wrong stateType => '+stateType
        self.diagname = 'OENI'+stateType
        self.availableStatistics = ['Mean', 'MS', 'RMS']


        self.label = '$OENI_{'+stateType+'}$',
        self.omm = 'om'+stateType
        self.sigmao = 'sigmao'+stateType
        self.requiredDiagnostics = [self.omm, self.sigmao]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        ommMean = self.retrieveDiagnosticStat(dfw, self.omm, 'Mean')
        ommMS = self.retrieveDiagnosticStat(dfw, self.omm, 'MS')
        ommRMS = self.retrieveDiagnosticStat(dfw, self.omm, 'RMS')

        sigmao = self.retrieveDiagnosticStat(dfw, self.sigmao, 'RMS')
        varianceo = self.retrieveDiagnosticStat(dfw, self.sigmao, 'MS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(ommMean, np.NaN)

        p = np.logical_and(np.isfinite(ommMean),
              np.logical_and(np.isfinite(ommMS), np.isfinite(sigmao)))
        if p.sum() > 0:
            # denominator zero check, q
            q = np.logical_and(bu.greatBound(np.absolute(sigmao), 0.0), p)

            #(1) ObsErrorNormalizedInnovation['Mean'] = Mean(OBS-MODEL) / RMS(SIGMAOBS)
            Stats['Mean'][q] = ommMean[q] / sigmao[q]

            #(2) ObsErrorNormalizedInnovation['MS'] = MS(OBS-MODEL) / MS(SIGMAOBS)
            Stats['MS'][q] = ommMS[q] / varianceo[q]

            #(3) ObsErrorNormalizedInnovation['RMS'] = RMS(OBS-MODEL) / RMS(SIGMAOBS)
            Stats['RMS'][q] = ommRMS[q] / sigmao[q]

        return self.templateDFfromStats(dfw, self.omm, Stats)


class InnovationRatio(DerivedDiagnostic):
    '''
    Calculates the ratio of Statistic(OMA)/Statistic(OMB)
    note: for Mean, the absolute ratio is calculated in order to make plots easier to decipher
    '''
    def __init__(self):
        self.diagname = 'InnovationRatio'
        self.availableStatistics = ['AbsMean', 'MS', 'RMS']

        self.label = '$\frac{OMA}{OMB}$',
        self.omb = 'omb'
        self.oma = 'oma'
        self.requiredDiagnostics = [self.omb, self.oma]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        ombMean = self.retrieveDiagnosticStat(dfw, self.omb, 'Mean')
        ombMS = self.retrieveDiagnosticStat(dfw, self.omb, 'MS')
        ombRMS = self.retrieveDiagnosticStat(dfw, self.omb, 'RMS')

        omaMean = self.retrieveDiagnosticStat(dfw, self.oma, 'Mean')
        omaMS = self.retrieveDiagnosticStat(dfw, self.oma, 'MS')
        omaRMS = self.retrieveDiagnosticStat(dfw, self.oma, 'RMS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(ombMean, np.NaN)

        p = np.logical_and(np.isfinite(omaMean),
              np.logical_and(np.isfinite(omaMS), np.isfinite(omaRMS)))
        if p.sum() > 0:
            #(1) InnovationRatio['Mean'] = |Mean(OBS-ANALYSIS) / Mean(OBS-BACKGROUND)|
            # denominator zero check, q
            q = np.logical_and(bu.greatBound(np.absolute(ombMean), 0.0), p)
            #Stats['AbsMean'][q] = omaMean[q] / ombMean[q]
            Stats['AbsMean'][q] = np.abs(omaMean[q] / ombMean[q])

            #(2) InnovationRatio['MS'] = MS(OBS-ANALYSIS) / MS(OBS-BACKGROUND)
            q = np.logical_and(bu.greatBound(np.absolute(ombMS), 0.0), p)
            Stats['MS'][q] = omaMS[q] / ombMS[q]

            #(3) InnovationRatio['RMS'] = RMS(OBS-ANALYSIS) / RMS(OBS-BACKGROUND)
            q = np.logical_and(bu.greatBound(np.absolute(ombRMS), 0.0), p)
            Stats['RMS'][q] = omaRMS[q] / ombRMS[q]

        return self.templateDFfromStats(dfw, self.omb, Stats)


class SpreadRatio(DerivedDiagnostic):
    '''
    Calculates the observation-space ensemble spread reduction ratio
    both for variance (MS) and standard deviation (RMS)
    '''
    def __init__(self):
        self.diagname = 'SpreadRatio'
        self.availableStatistics = ['VAR', 'STD']

        #self.label = '$\frac{\sigma_{h_a}}{\sigma_{h_b}}$',
        self.label = '$SpreadRatio$',
        self.sigmab = 'sigmab'
        self.sigmaa = 'sigmaa'
        self.requiredDiagnostics = [self.sigmab, self.sigmaa]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        sigmabMS = self.retrieveDiagnosticStat(dfw, self.sigmab, 'MS')
        sigmabRMS = self.retrieveDiagnosticStat(dfw, self.sigmab, 'RMS')

        sigmaaMS = self.retrieveDiagnosticStat(dfw, self.sigmaa, 'MS')
        sigmaaRMS = self.retrieveDiagnosticStat(dfw, self.sigmaa, 'RMS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(sigmabMS, np.NaN)

        p = np.logical_and(np.isfinite(sigmaaMS), np.isfinite(sigmabMS))
        if p.sum() > 0:
            #(1) SpreadRatio['VAR'] = MS(SIGMAA) / MS(SIGMAB)
            q = np.logical_and(bu.greatBound(np.absolute(sigmabMS), 0.0), p)
            Stats['VAR'][q] = sigmaaMS[q] / sigmabMS[q]

            #(2) SpreadRatio['STD'] = RMS(SIGMAA) / RMS(SIGMAB)
            q = np.logical_and(bu.greatBound(np.absolute(sigmabRMS), 0.0), p)
            Stats['STD'][q] = sigmaaRMS[q] / sigmabRMS[q]

        return self.templateDFfromStats(dfw, self.sigmab, Stats)


## Table of available diagnostics
# format of each entry:
#    diagName: {
#        variable: variable name or class used to calculate this diag
#        iter: iteration of "self" type variables
#          ('bg', 'an', int, default: None OR vu.bgIter depending on appName)
#        vu.mean (optional): whether to apply this diagnostic to the mean state (default: True)
#        vu.ensemble (optional): whether this diagnostic requires variables from ensemble IODA files (default: False)
#        onlyObsSpaces (optional): list of ObsSpaces for which this diag applies
#            see config.DiagSpaceConfig keys for available options
#        analyze (optional): whether to analyze this diagnostic (defualt: True).  Useful for turning
#            off analyses for diagnostics that are not of particular interest or that are not available
#            for particular experiments, e.g., sigmaf for a single-state cycling experiment
#        derived (optional): whether the diagnostic is derived from statistics of other diagnostics
#            if True, then will be calculated in analysis step, not statistics calculation step
#            (default: False)
#        DerivedDiagnostic object (optional): if derived, an object defining the function used to
#            calculate this diagnostic. Must have an 'evaluate' method that takes a
#            StatisticsDatabase.DFWrapper object as the only argument.
#        requiredDiagnostics (optional): list of other availableDiagnostics that are required in order
#            to generate the DerivedDiagnostic
#        availableStatistics (optional): limited set of statistics provided by this diagnostic
#            (default: application dependent, see Analyses module)
#        selectedStatistics (optional): statistics selected for this diagnostic
#            (default: same as availableStatistics)
#        label (optional): label for plot axes (TODO: hookup in Analyses.py)
#}
#TODO: refactor the 'analyze' config entry under availableDiagnostics so that it is automated and
#      configurable at a top level

# ideal σ_o
biso_ = idealSigmao('b')
aiso_ = idealSigmao('a')
fiso_ = idealSigmao('f')

# consistency ratio
bcr_ = ObsSpaceConsistencyRatio('b')
acr_ = ObsSpaceConsistencyRatio('a')
fcr_ = ObsSpaceConsistencyRatio('f')

# innovation (omb, oma, omf) normalized by σ_o
boeni_ = ObsErrorNormalizedInnovation('b')
aoeni_ = ObsErrorNormalizedInnovation('a')
foeni_ = ObsErrorNormalizedInnovation('f')

# ratio of oma / omb
innovratio_ = InnovationRatio()

# ratio of σ_h(x_a) / σ_h(x_b)
spreadratio_ = SpreadRatio()

# components of departure spread equation:
#   E[dd^T] = HBH^T + R,
#departurespread_ = DepartureSTD('f')
departurespread_ = DepartureRMS('f')
ensspread_ = EnsembleSpread('f')
obserror_ = ObsError('f')
totalspread_ = TotalSpread('f')

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
        'analyze': True,
        'variable': ObsMinusModel,
        'label': '$y - x_f$',
    },
    'amb': {
        'variable': AnalysisMinusBackground,
        'iter': 'an',
        'label': '$x_a - x_b$',
    },
    'mmgfsan': {
        'offline': True,
        'label': '$x - x_{a,GFS}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'rltv_mmgfsan': {
        'offline': True,
        'label': '$\frac{x - x_{a,GFS}}{x_{a,GFS}}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'log_mogfsan': {
        'offline': True,
        'label': '$\log{\frac{x}{x_{a,GFS}}}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
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
        'label': '$(y - x_b) / y$',
    },
    'rltv_oma': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': '$(y - x_a) / y$',
    },
    'rltv_omf': {
        'analyze': True,
        'variable': RelativeObsMinusModel,
        'label': '$(y - x_f) / y$',
    },
    'rltv_amb': {
        'variable': RelativeAnalysisMinusBackground,
        'iter': 'an',
        'label': '$(x_a - x_b) / y$',
    },
    'rltv_omb_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'bg',
        'label': '$(y - x_b) / y$',
    },
    'rltv_oma_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'an',
        'label': '$(y - x_a) / y$',
    },
    'amb_o_omb': {
        'variable': AnalysisMinusBackgroundOverObsMinusBackground,
        'iter': 'an',
        'label': '$(x_a - x_b) / (y - x_b)$',
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
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmab': {
        'variable': bu.STDofHofX,
        'iter': 'bg',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_b}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaoa': {
        'variable': vu.selfErrorValue,
        'iter': 'an',
        'analyze': False,
        'label': '$\sigma_o$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaa': {
        'variable': bu.STDofHofX,
        'iter': 'an',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_a}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaof': {
        'variable': vu.selfErrorValue,
        'analyze': False,
        'label': '$\sigma_o$',
        'selectedStatistics': ['RMS', 'STD'],
    },
    #NOTE: a failure results when 'analyze' is True under sigmaf, any one of the experiments
    # does not have sigmaf (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmaf
    'sigmaf': {
        'variable': bu.STDofHofX,
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': '$\sigma_{h_f}$',
        'selectedStatistics': su.sigmaStatistics,
    },
# cloud-related diagnostics
    'SCI-'+bu.OkamotoMethod: {
        'variable': bu.SCIOkamoto,
        'analyze': False,
        'onlyObsSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'ACI-'+bu.MZ19Method: {
        'variable': bu.ACIMZ19,
        'analyze': False,
        'onlyObsSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'MCI': {
        'variable': bu.MCI,
        'analyze': False,
        'onlyObsSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'CFy': {
        'variable': vu.cldfracMeta,
        'analyze': False,
        'onlyObsSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'ABEILambda': {
        'variable': bu.ABEILambda,
        'onlyObsSpaces': ['abi_g16', 'ahi_himawari8'],
        'label': '$\lambda_{ABEI}$',
    },
# DerivedDiagnostics
    'ideal-sigmaob': {
        'iter': 'bg',
        'analyze': False,
        'DerivedDiagnostic': biso_,
        'requiredDiagnostics': biso_.requiredDiagnostics,
        'availableStatistics': biso_.availableStatistics,
        'label': biso_.label,
    },
    'ideal-sigmaoa': {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': aiso_,
        'requiredDiagnostics': aiso_.requiredDiagnostics,
        'availableStatistics': aiso_.availableStatistics,
        'label': aiso_.label,
    },
    #NOTE: a failure results when 'analyze' is True under ideal-sigmaof, any one of the experiments
    # does not have sigmaf (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmaf
    'ideal-sigmaof': {
        'DerivedDiagnostic': fiso_,
        'analyze': False,
        'requiredDiagnostics': fiso_.requiredDiagnostics,
        'availableStatistics': fiso_.availableStatistics,
        'label': fiso_.label,
    },
    'CRyb': {
        'iter': 'bg',
        'analyze': False,
        'DerivedDiagnostic': bcr_,
        'requiredDiagnostics': bcr_.requiredDiagnostics,
        'availableStatistics': bcr_.availableStatistics,
        'label': bcr_.label,
    },
    'CRya': {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': acr_,
        'requiredDiagnostics': acr_.requiredDiagnostics,
        'availableStatistics': acr_.availableStatistics,
        'label': acr_.label,
    },
    #NOTE: a failure results when 'analyze' is True under CRyf, any one of the experiments
    # does not have sigmaf (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmaf
    'CRyf': {
        'DerivedDiagnostic': fcr_,
        'analyze': False,
        'requiredDiagnostics': fcr_.requiredDiagnostics,
        'availableStatistics': fcr_.availableStatistics,
        'label': fcr_.label,
    },
    'OENIb': {
        'iter': 'bg',
        'DerivedDiagnostic': boeni_,
        'analyze': False,
        'requiredDiagnostics': boeni_.requiredDiagnostics,
        'availableStatistics': boeni_.availableStatistics,
        'label': boeni_.label,
    },
    'OENIa': {
        'iter': 'an',
        'DerivedDiagnostic': aoeni_,
        'analyze': False,
        'requiredDiagnostics': aoeni_.requiredDiagnostics,
        'availableStatistics': aoeni_.availableStatistics,
        'label': aoeni_.label,
    },
    'OENIf': {
        'DerivedDiagnostic': foeni_,
        'analyze': False,
        'requiredDiagnostics': foeni_.requiredDiagnostics,
        'availableStatistics': foeni_.availableStatistics,
        'label': foeni_.label,
    },
    'InnovationRatio': {
        'DerivedDiagnostic': innovratio_,
        'analyze': False,
        'requiredDiagnostics': innovratio_.requiredDiagnostics,
        'availableStatistics': innovratio_.availableStatistics,
        'label': innovratio_.label,
    },
    'SpreadRatio': {
        'DerivedDiagnostic': spreadratio_,
        'analyze': False,
        'requiredDiagnostics': spreadratio_.requiredDiagnostics,
        'availableStatistics': spreadratio_.availableStatistics,
        'label': spreadratio_.label,
    },
    departurespread_.diagname: {
        'DerivedDiagnostic': departurespread_,
        'analyze': False,
        'requiredDiagnostics': departurespread_.requiredDiagnostics,
        'availableStatistics': departurespread_.availableStatistics,
        'label': departurespread_.label,
    },
    ensspread_.diagname: {
        'DerivedDiagnostic': ensspread_,
        'analyze': False,
        'requiredDiagnostics': ensspread_.requiredDiagnostics,
        'availableStatistics': ensspread_.availableStatistics,
        'label': ensspread_.label,
    },
    obserror_.diagname: {
        'DerivedDiagnostic': obserror_,
        'analyze': False,
        'requiredDiagnostics': obserror_.requiredDiagnostics,
        'availableStatistics': obserror_.availableStatistics,
        'label': obserror_.label,
    },
    totalspread_.diagname: {
        'DerivedDiagnostic': totalspread_,
        'analyze': False,
        'requiredDiagnostics': totalspread_.requiredDiagnostics,
        'availableStatistics': totalspread_.availableStatistics,
        'label': totalspread_.label,
    },
}

#TODO: have this function return a list of diagnosticConfiguration or Diagnostic (new class) objects
#      instead of a list of dicts
def diagnosticConfigs(diagnosticNames_, ObsSpaceName, includeEnsembleDiagnostics=True,
                      selectedStatistics = su.allFileStats, fileFormat=vu.hdfFileFormat):

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
        config['offline'] = config.get('offline', False)
        config[vu.mean] = config.get(vu.mean, True)
        config[vu.ensemble] = (config.get(vu.ensemble, False) and includeEnsembleDiagnostics)
        config['onlyObsSpaces'] = config.get('onlyObsSpaces',[])
        config['availableStatistics'] = config.get('availableStatistics', selectedStatistics)
        config['selectedStatistics'] = config.get('selectedStatistics', config['availableStatistics'])
        config['label'] = config.get('label',diagnosticName)
        config['requiredDiagnostics'] = set(config.get('requiredDiagnostics',[]))

        # diagnosticConfig is undefined for the following cases
        if (not config[vu.mean] and not config[vu.ensemble]): continue
        if (len(config['onlyObsSpaces']) > 0 and
            ObsSpaceName not in config['onlyObsSpaces']): continue

        config['osName'] = ObsSpaceName
        config['fileFormat'] = fileFormat
        # add this diagnosticConfig to the list
        diagnosticConfigs[diagnosticName] = deepcopy(config)

        outerIterStr = config.get('iter', None)
        if outerIterStr is  None or outerIterStr == '':
            if appName == 'hofx':
                outerIter = None
            else:
                outerIter = vu.bgIter
        else:
            if outerIterStr == 'bg':
                outerIter = vu.bgIter
            elif outerIterStr == 'an':
                outerIter = anIter
            elif pu.isint(outerIterStr):
                outerIter = outerIterStr
            else:
                _logger.error('outerIter is undefined: '+outerIterStr)

        diagnosticConfigs[diagnosticName]['outerIter'] = outerIter

        if ('DerivedDiagnostic' in config): continue
        if config['offline']: continue
        diagnosticConfigs[diagnosticName]['BinFunction'] = bu.BinFunctionWrapper(config)

    return diagnosticConfigs
