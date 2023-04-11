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
        self.baseVars.append(vu.selfOMMValue)
        self.baseVars.append(vu.selfHofXValue)

    def evaluate(self, dbVals, insituParameters):
        omm_bc = dbVals[insituParameters[vu.selfOMMValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]
        return np.add(mod, omm_bc)


class BiasCorrection:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMMValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        omm_bc = dbVals[insituParameters[vu.selfOMMValue]]
        mod = dbVals[insituParameters[vu.selfHofXValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]
        return np.subtract(np.add(mod, omm_bc), obs)


class BiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMMValue)

    def evaluate(self, dbVals, insituParameters):
        omm_bc = deepcopy(dbVals[insituParameters[vu.selfOMMValue]])
        return omm_bc


class ExpectedDoaDob:
    '''Sqrt of diagonal of Doesroziers et al. (2005)'s E[d^o_a d^o_b^T] ~ R'''
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMBValue)
        self.baseVars.append(vu.selfOMAValue)

    def evaluate(self, dbVals, insituParameters):
        omb_bc = dbVals[insituParameters[vu.selfOMBValue]]
        oma_bc = dbVals[insituParameters[vu.selfOMAValue]]
        return negativeSafeSqrt(np.multiply(oma_bc, omb_bc))


class ExpectedDabDob:
    '''Sqrt of diagonal of Doesroziers et al. (2005)'s E[d^a_b d^o_b^T] ~ HBH^T'''
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMBValue)
        self.baseVars.append(vu.anHofXValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.anHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        dab = np.subtract(ana, bak)

        omb_bc = dbVals[insituParameters[vu.selfOMBValue]]
        return negativeSafeSqrt(np.multiply(dab, omb_bc))


class ExpectedRelativeDoaDob:
    '''Sqrt of diagonal of Doesroziers et al. (2005)'s E[d^o_a d^o_b^T] ~ R, normalized by observed value'''
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMBValue)
        self.baseVars.append(vu.selfOMAValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        omb_bc = dbVals[insituParameters[vu.selfOMBValue]]
        oma_bc = dbVals[insituParameters[vu.selfOMAValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]
        return np.multiply(np.divide(negativeSafeSqrt(np.multiply(oma_bc, omb_bc)), obs), 100.)


class ExpectedRelativeDabDob:
    '''Sqrt of diagonal of Doesroziers et al. (2005)'s E[d^a_b d^o_b^T] ~ HBH^T, normalized by observed value'''
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMBValue)
        self.baseVars.append(vu.anHofXValue)
        self.baseVars.append(vu.bgHofXValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.anHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        dab = np.subtract(ana, bak)
        obs = dbVals[insituParameters[vu.selfObsValue]]

        omb_bc = dbVals[insituParameters[vu.selfOMBValue]]
        return np.multiply(np.divide(negativeSafeSqrt(np.multiply(dab, omb_bc)), obs), 100.)


class NonBiasCorrectedObsMinusModel(bu.InsituLocFunction):
    def _initBaseVars(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMMValue)
        self.baseVars.append(vu.selfBCValue)

    def _get(self):
        omm_bc = self._variables[vu.selfOMMValue]
        bc = self._variables[vu.selfBCValue]

        return np.subtract(omm_bc, bc)


class RelativeBiasCorrectedObsMinusModel:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfOMMValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        omm_bc = deepcopy(dbVals[insituParameters[vu.selfOMMValue]])
        obs = dbVals[insituParameters[vu.selfObsValue]]

        valid = bu.greatBound(np.abs(obs), 0.0, False)
        omm_bc[valid] = np.divide(omm_bc[valid], obs[valid])
        omm_bc[~valid] = np.NaN

        return np.multiply(omm_bc, 100.0)


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
        self.baseVars.append(vu.anHofXValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.anHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        return np.subtract(ana, bak)


class RelativeAnalysisMinusBackground:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.anHofXValue)
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.anHofXValue]]
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
        self.baseVars.append(vu.anHofXValue)
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.bgHofXValue)

    def evaluate(self, dbVals, insituParameters):
        ana = dbVals[insituParameters[vu.anHofXValue]]
        bak = dbVals[insituParameters[vu.bgHofXValue]]
        obs = dbVals[insituParameters[vu.selfObsValue]]

        AMBoOMB = np.subtract(ana, bak)
        OMB = np.subtract(obs, bak)
        valid = bu.greatBound(np.abs(OMB), 0.0, False)

        AMBoOMB[valid] = np.divide(AMBoOMB[valid], OMB[valid])
        AMBoOMB[~valid] = np.NaN

        return AMBoOMB


class RelativeError:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfErrorValue)
        self.baseVars.append(vu.selfObsValue)

    def evaluate(self, dbVals, insituParameters):
        err = deepcopy(dbVals[insituParameters[vu.selfErrorValue]])
        obs = dbVals[insituParameters[vu.selfObsValue]]

        valid = bu.greatBound(np.abs(obs), 0.0, False)
        err[valid] = np.divide(err[valid], obs[valid])
        err[~valid] = np.NaN

        return np.multiply(err, 100.0)


class DerivedDiagnostic:
    '''
    base class for derived combinations of diagnostics and statistics
    The main functionality is located in the evaluate method, which is specified
    in the derived class.
    '''
    def __init__(self):
        '''
        virtual method
        '''
        raise NotImplementedError()
        # Child.__init__() must define the following:
        #self.availableStatistics = su.allFileStats
        #self.diagname = 'DerivedDiagnostic'

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
    availableStatistics = ['RMS']
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'idealSigmao: wrong stateType => '+stateType
        self.diagname = 'ideal-sigmao'+stateType

        self.label = r'$\sigma_{o'
        if stateType != 'f':
          self.label += r','+stateType
        self.label += r'}^*$'

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
            # idealSigmao(RMS) = SQRT( [MS(OBS-MODEL) - MS(SIGMAMOD)] )
            Stats['RMS'][p] = omm[p] - sigmam[p]
            Stats['RMS'] = negativeSafeSqrt(Stats['RMS'])

        return self.templateDFfromStats(dfw, self.omm, Stats)


class idealRelativeSigmao(DerivedDiagnostic):
    '''
    Calculates ideal sigmao for the
    background, analysis, or forecast states
    '''
    availableStatistics = ['E']
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'idealRelativeSigmao: wrong stateType => '+stateType
        self.diagname = 'ideal-rltv_sigmao'+stateType

        self.label = r'$\frac{\sigma_{o'
        if stateType != 'f':
          self.label += r','+stateType
        self.label += r'}^*}{\bar{y}} \times 100$'

        self.omm = 'om'+stateType
        self.rltvomm = 'rltv_om'+stateType
        self.sigmam = 'sigma'+stateType

        self.requiredDiagnostics = [self.omm, self.rltvomm, self.sigmam]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        omm = self.retrieveDiagnosticStat(dfw, self.omm, 'MS')
        rltvomm = self.retrieveDiagnosticStat(dfw, self.rltvomm, 'RMS')
        sigmam = self.retrieveDiagnosticStat(dfw, self.sigmam, 'MS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(omm, np.NaN)

        p = np.logical_and(np.isfinite(omm), np.isfinite(sigmam))
        if p.sum() > 0:
            # yy is approximately equal to (obs)^2
            yy = np.divide(omm, np.square(np.divide(rltvomm, 100.)))

            # idealRelativeSigmao(E) = SQRT( [MS(OBS-MODEL) - MS(SIGMAMOD)] / OBS^2 )
            Stats['E'][p] = np.divide(omm[p] - sigmam[p], yy[p])
            Stats['E'] = np.multiply(negativeSafeSqrt(Stats['E']), 100.)

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
    availableStatistics = [CRd, CRo, CRh, invCRo, invCRh]

    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'ObsSpaceConsistencyRatio: wrong stateType => '+stateType
        self.diagname = 'CRy'+stateType

        self.label = r'$CR_{y'
        if stateType != 'f':
          self.label += r','+stateType
        self.label += r'}$'

        self.label = r'$CR_{y,'+stateType+'}$'
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


class ModelSpaceConsistencyRatio(DerivedDiagnostic):
    '''
    Calculates multiple forms of the model-space ensemble consistency ratio for the
    background, analysis, inflated, or forecast states
    '''
    CR = 'CR'
    invCR = 'invCR'
    availableStatistics = [CR, invCR]

    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'inf', 'f'], 'ModelSpaceConsistencyRatio: wrong stateType => '+stateType
        self.diagname = 'CRx'+stateType

        self.label = r'$CR_{x'
        if stateType != 'f':
          self.label += r','+stateType
        self.label += r'}$'

        self.mmref = 'mmgfsan'
        self.sigmam = 'sigmax'+stateType

        self.requiredDiagnostics = [self.mmref, self.sigmam]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        mmref = self.retrieveDiagnosticStat(dfw, self.mmref, 'MS')
        sigmam = self.retrieveDiagnosticStat(dfw, self.sigmam, 'MS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(mmref, np.NaN)

        p = np.logical_and(np.isfinite(mmref), np.isfinite(sigmam))
        if p.sum() > 0:
            #For each subpopulation, the consistency ratio can be written two different
            # ways, depending on the quantity of interest in the analysis. In the
            # following formulae, MS ≡ Mean Squared Value.

            #(1) CR = SQRT( MS(SIGMAMOD) / MS(MODEL-REF) )
            # denominator zero check, q
            q = np.logical_and(bu.greatBound(np.absolute(mmref), 0.0), p)
            Stats[self.CR][q] = sigmam[q] / mmref[q]
            # utility: reveals actual spread skill (too large or too small by specific factor)

            #(2) invCRh = SQRT( [MS(OBS-MODEL) - MS(SIGMAOBS)] / MS(SIGMAMOD) )
            q = np.logical_and(bu.greatBound(np.absolute(sigmam), 0.0), p)
            Stats[self.invCR][q] = mmref[q] / sigmam[q]
            # utility: gives scaling factors for adjusting the spread of self.stateType

            # perform the sqrt operation for all CR's where positive-semi-definite
            for statistic, Stat in Stats.items():
                Stats[statistic] = negativeSafeSqrt(Stat)

        return self.templateDFfromStats(dfw, self.mmref, Stats)


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
        self.diagname = self.diagName
        if stateType != 'f':
          self.diagname += '_'+stateType
        self.label = r'$'+self.diagName+'$'
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
    statName = 'rms'
    diagName = vu.DiagnosticVars[vu.EddT]
    sourceStatistics = ['STD']
    sourceDiagnosticPrefixes = ['om']


class DepartureRMS(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the LHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'rms'
    diagName = vu.DiagnosticVars[vu.EddT]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['om']

class EnsembleSpread(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the 1st RHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'rms'
    diagName = vu.DiagnosticVars[vu.HBHT]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['sigma']


class ObsError(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the 2nd RHS term of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'rms'
    diagName = vu.DiagnosticVars[vu.R]
    sourceStatistics = ['RMS']
    sourceDiagnosticPrefixes = ['sigmao']


class TotalSpread(SumDiagnosticsStatisticsPairs):
    '''
    Calculates the sqrt of the mean of the full RHS of the spread expectation equation:
      E[dd^T] = HBH^T + R,
    '''
    statName = 'rms'
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
    availableStatistics = ['Mean', 'MS', 'RMS']
    def __init__(self, stateType):
    # stateType - model state diagnostic type (either 'b' for background, 'a' for analysis, or 'f' for forecast)
        assert stateType in ['b', 'a', 'f'], 'ObsErrorNormalizedInnovation: wrong stateType => '+stateType
        self.diagname = 'OENI'+stateType

        self.label = r'$OENI'
        if stateType != 'f':
          self.label += r'_{'+stateType+'}'
        self.label += r'$'

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
    diagname = 'InnovationRatio'
    availableStatistics = ['AbsMean', 'MS', 'RMS']
    def __init__(self):
        self.label = r'$\frac{O-A}{O-B}$'
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
    Calculates the ensemble spread ratio
    both for variance (MS) and standard deviation (RMS)
    '''
    availableStatistics = ['STD']
    def __init__(self, diagname, sigmab, sigmaa):
        self.diagname = diagname
        self.label = r'$'+diagname+'$'
        self.sigmab = sigmab
        self.sigmaa = sigmaa
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
            #q = np.logical_and(bu.greatBound(np.absolute(sigmabMS), 0.0), p)
            #Stats['VAR'][q] = sigmaaMS[q] / sigmabMS[q]

            #(2) SpreadRatio['STD'] = RMS(SIGMAA) / RMS(SIGMAB)
            q = np.logical_and(bu.greatBound(np.absolute(sigmabRMS), 0.0), p)
            Stats['STD'][q] = sigmaaRMS[q] / sigmabRMS[q]

        return self.templateDFfromStats(dfw, self.sigmab, Stats)


class ApproxDoaDob(DerivedDiagnostic):
    '''
    Calculates sqrt(RMS(OMA) * RMS(OMB)), rough approximation of Doesroziers et al. (2005)'s E[d^o_a d^o_b^T] ~ R
    '''
    diagname = 'ApproxDoaDob'
    availableStatistics = ['RMS']
    def __init__(self):
        self.label = r'$\sqrt{ rms\left(O-A\right)\times rms\left(O-B\right)}$'
        self.omb = 'omb'
        self.oma = 'oma'
        self.requiredDiagnostics = [self.omb, self.oma]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        ombRMS = self.retrieveDiagnosticStat(dfw, self.omb, 'RMS')
        omaRMS = self.retrieveDiagnosticStat(dfw, self.oma, 'RMS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(ombRMS, np.NaN)

        p = np.logical_and(np.isfinite(ombRMS), np.isfinite(omaRMS))

        if p.sum() > 0:
            # ApproxDoaDob['RMS'] = sqrt(RMS(OMA) * RMS(OMB))
            Stats['RMS'][p] = np.sqrt(omaRMS[p] * ombRMS[p])

        return self.templateDFfromStats(dfw, self.omb, Stats)


class ApproxRelativeDoaDob(DerivedDiagnostic):
    '''
    Calculates sqrt(RMS(OMA) * RMS(OMB)), rough approximation of Doesroziers et al. (2005)'s E[d^o_a d^o_b^T] ~ R, normalized by approx observed value
    '''
    diagname = 'ApproxRelativeDoaDob'
    availableStatistics = ['RMS']
    def __init__(self):
        self.label = r'$\frac{\sqrt{ rms\left(O-A\right)\times rms\left(O-B\right)}}{\bar{O}} \times 100$'
        self.omb = 'omb'
        self.rltvomb = 'rltv_omb'
        self.oma = 'oma'
        self.requiredDiagnostics = [self.omb, self.oma, self.rltvomb]

    def evaluate(self, dfw):
        diagNamesAvailableInSlice = dfw.levels('diagName')
        for diag in self.requiredDiagnostics:
            if diag not in diagNamesAvailableInSlice: return None

        # get the statistics needed as numpy arrays
        ombRMS = self.retrieveDiagnosticStat(dfw, self.omb, 'RMS')
        rltvombRMS = self.retrieveDiagnosticStat(dfw, self.rltvomb, 'RMS')
        omaRMS = self.retrieveDiagnosticStat(dfw, self.oma, 'RMS')

        Stats = {}
        for statName in self.availableStatistics:
            Stats[statName] = np.full_like(ombRMS, np.NaN)

        p = np.logical_and(np.isfinite(ombRMS), np.isfinite(omaRMS))

        if p.sum() > 0:
            # yy is approximately equal to (obs)^2
            yy = np.divide(np.square(ombRMS), np.square(np.divide(rltvombRMS, 100.)))

            # ApproxRelativeDoaDob['RMS'] = sqrt(RMS(OMA) * RMS(OMB))
            Stats['RMS'][p] = np.multiply(np.sqrt(np.divide(omaRMS[p] * ombRMS[p], yy[p])), 100.)

        return self.templateDFfromStats(dfw, self.omb, Stats)

## Table of available diagnostics
# format of each entry:
#    diagName: {
#        variable: variable name or class used to calculate this diag
#        iter: iteration of "self" type variables
#          ('bg', 'an', int, default: None OR vu.bgIter depending on appName)
#        vu.mean (optional): whether to apply this diagnostic to the mean state (default: True)
#        vu.ensemble (optional): whether this diagnostic requires variables from ensemble IODA files (default: False)
#        onlyDiagSpaces (optional): list of DiagSpaces for which this diag applies
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

birso_ = idealRelativeSigmao('b')
airso_ = idealRelativeSigmao('a')
firso_ = idealRelativeSigmao('f')

# consistency ratio
bcry_ = ObsSpaceConsistencyRatio('b')
acry_ = ObsSpaceConsistencyRatio('a')
fcry_ = ObsSpaceConsistencyRatio('f')

bcrx_ = ModelSpaceConsistencyRatio('b')
acrx_ = ModelSpaceConsistencyRatio('a')
icrx_ = ModelSpaceConsistencyRatio('inf')
fcrx_ = ModelSpaceConsistencyRatio('f')

# innovation (omb, oma, omf) normalized by σ_o
boeni_ = ObsErrorNormalizedInnovation('b')
aoeni_ = ObsErrorNormalizedInnovation('a')
foeni_ = ObsErrorNormalizedInnovation('f')

# ratio of oma / omb
innovratio_ = InnovationRatio()

# ratio of σ_h(x_a) / σ_h(x_b)
sry_ = SpreadRatio('SRy', 'sigmab', 'sigmaa')

# ratio of σ_{x_a} / σ_{x_b}
srxeda_ = SpreadRatio('SRx-eda', 'sigmaxb', 'sigmaxa')

# ratio of σ_{x_inf} / σ_{x_a}
srxrtpp_ = SpreadRatio('SRx-rtpp', 'sigmaxa', 'sigmaxinf')

# approximate desroziers diagnostic: sqrt(diag((O-A)(O-B)^T))
approxdoadob_ = ApproxDoaDob()
approxrltvdoadob_ = ApproxRelativeDoaDob()

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
        'label': r'$y_{bc}$',
    },
    'omb': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'bg',
        'label': r'$O-B$',
    },
    'oma': {
        'variable': BiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': r'$O-A$',
    },
    'doadob': {
        'variable': ExpectedDoaDob,
        'label': r'$\sqrt{(O-A)^T(O-B)}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'dabdob': {
        'variable': ExpectedDabDob,
        'label': r'$\sqrt{(A-B)^T(O-B)}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'omf': {
        'analyze': True,
        'variable': ObsMinusModel,
        'label': r'$O-F$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'amb': {
        'variable': AnalysisMinusBackground,
        'label': r'$A-B$',
    },
    'mmgfsan': {
        'offline': True,
        'label': r'$\delta{x_{GFSa}}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'rltv_mmgfsan': {
        'offline': True,
        'analyze': False,
        'label': r'$\frac{\delta{x_{GFSa}}}{x_{GFSa}}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'log_mogfsan': {
        'offline': True,
        'analyze': False,
        'label': r'$\log{\frac{x}{x_{GFSa}}}$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    #NOTE: a failure results when 'analyze' is True under sigmax*, any one of the experiments
    # does not have sigmax* (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmax*
    'sigmaxf': {
        'offline': True,
        'analyze': True,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{x}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaxb': {
        'offline': True,
        'analyze': True,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{x_b}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaxa': {
        'offline': True,
        'analyze': True,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{x_a}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaxinf': {
        # inf == inflated analysis, e.g., after RTPP
        'offline': True,
        'analyze': True,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{x_{inf}}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'omb_nobc': {
        'variable': NonBiasCorrectedObsMinusModel,
        'iter': 'bg',
        'label': r'$O_{orig}-B$',
    },
    'oma_nobc': {
        'variable': NonBiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': r'$O_{orig}-A$',
    },
    'rltv_omb': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'bg',
        'label': r'$\frac{O-B}{O} \times 100$',
    },
    'rltv_oma': {
        'variable': RelativeBiasCorrectedObsMinusModel,
        'iter': 'an',
        'label': r'$\frac{O-A}{O} \times 100$',
    },
    'rltv_doadob': {
        'variable': ExpectedRelativeDoaDob,
        'label': r'$\frac{\sqrt{(O-A)^T(O-B)}}{O} \times 100$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'rltv_dabdob': {
        'variable': ExpectedRelativeDabDob,
        'label': r'$\frac{\sqrt{(A-B)^T(O-B)}}{O} \times 100$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'rltv_omf': {
        'analyze': True,
        'variable': RelativeObsMinusModel,
        'label': r'$\frac{O-F}{O} \times 100$',
        'selectedStatistics': ['Mean', 'RMS', 'STD'],
    },
    'rltv_amb': {
        'variable': RelativeAnalysisMinusBackground,
        'label': r'$\frac{A-B}{y} \times 100$',
    },
    'rltv_omb_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'bg',
        'label': r'$\frac{O-B}{O} \times 100$',
    },
    'rltv_oma_nobc': {
        'variable': RelativeObsMinusModel,
        'iter': 'an',
        'label': r'$\frac{O-A}{O} \times 100$',
    },
    'amb_o_omb': {
        'variable': AnalysisMinusBackgroundOverObsMinusBackground,
        'label': r'$\frac{A-B}{O-B}$',
    },
    'obs': {
        'variable': vu.selfObsValue,
        'label': r'$y$',
    },
    'obs_bc': {
        'variable': BiasCorrectedObs,
        'label': r'$y$',
    },
    'h(x)': {
        'variable': vu.selfHofXValue,
        'label': r'$h(x)$',
    },
    'bak': {
        'variable': vu.selfHofXValue,
        'iter': 'bg',
        'label': r'$h(x_b)$',
    },
    'ana': {
        'variable': vu.selfHofXValue,
        'iter': 'an',
        'label': r'$h(x_a)$',
    },
    'sigmaob': {
        'variable': vu.selfErrorValue,
        'iter': 'bg',
        'analyze': True,
        'label': r'$\sigma_{o,b}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'rltv_sigmaob': {
        'variable': RelativeError,
        'iter': 'bg',
        'analyze': True,
        'label': r'$\frac{\sigma_{o,b}}{y} \times 100$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmab': {
        'variable': bu.STDofHofX,
        'iter': 'bg',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{h(x_b)}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaoa': {
        'variable': vu.selfErrorValue,
        'iter': 'an',
        'analyze': False,
        'label': r'$\sigma_{o,a}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'rltv_sigmaoa': {
        'variable': RelativeError,
        'iter': 'an',
        'analyze': False,
        'label': r'$\frac{\sigma_{o,a}}{y} \times 100$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaa': {
        'variable': bu.STDofHofX,
        'iter': 'an',
        'analyze': False,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{h(x_a)}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'sigmaof': {
        'variable': vu.selfErrorValue,
        'analyze': True,
        'label': r'$\sigma_{o}$',
        'selectedStatistics': su.sigmaStatistics,
    },
    'rltv_sigmaof': {
        'variable': RelativeError,
        'analyze': True,
        'label': r'$\frac{\sigma_{o}}{y} \times 100$',
        'selectedStatistics': su.sigmaStatistics,
    },
    #NOTE: a failure results when 'analyze' is True under sigmaf, any one of the experiments
    # does not have sigmaf (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmaf
    'sigmaf': {
        'variable': bu.STDofHofX,
        'analyze': True,
        vu.mean: False,
        vu.ensemble: True,
        'label': r'$\sigma_{h(x_f)}$',
        'selectedStatistics': su.sigmaStatistics,
    },
# cloud-related diagnostics
    'SCI-'+bu.OkamotoMethod: {
        'variable': bu.SCIOkamoto,
        'analyze': False,
        'onlyDiagSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'ACI-'+bu.MZ19Method: {
        'variable': bu.ACIMZ19,
        'analyze': False,
        'onlyDiagSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'MCI': {
        'variable': bu.MCI,
        'analyze': False,
        'onlyDiagSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'CFy': {
        'variable': vu.cldfracMeta,
        'analyze': False,
        'onlyDiagSpaces': ['abi_g16', 'ahi_himawari8'],
        'selectedStatistics': ['Mean', 'STD'],
    },
    'ABEILambda': {
        'variable': bu.ABEILambda,
        'onlyDiagSpaces': ['abi_g16', 'ahi_himawari8'],
        'label': r'$\lambda_{ABEI}$',
    },
# DerivedDiagnostics
    biso_.diagname: {
        'iter': 'bg',
        'analyze': False,
        'DerivedDiagnostic': biso_,
        'requiredDiagnostics': biso_.requiredDiagnostics,
        'availableStatistics': biso_.availableStatistics,
        'label': biso_.label,
    },
    birso_.diagname: {
        'iter': 'bg',
        'analyze': False,
        'DerivedDiagnostic': birso_,
        'requiredDiagnostics': birso_.requiredDiagnostics,
        'availableStatistics': birso_.availableStatistics,
        'label': birso_.label,
    },
    aiso_.diagname: {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': aiso_,
        'requiredDiagnostics': aiso_.requiredDiagnostics,
        'availableStatistics': aiso_.availableStatistics,
        'label': aiso_.label,
    },
    airso_.diagname: {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': airso_,
        'requiredDiagnostics': airso_.requiredDiagnostics,
        'availableStatistics': airso_.availableStatistics,
        'label': airso_.label,
    },
    #NOTE: a failure results when 'analyze' is True under ideal-sigmaof, any one of the experiments
    # does not have sigmaf (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have sigmaf
    fiso_.diagname: {
        'DerivedDiagnostic': fiso_,
        'analyze': False,
        'requiredDiagnostics': fiso_.requiredDiagnostics,
        'availableStatistics': fiso_.availableStatistics,
        'label': fiso_.label,
    },
    firso_.diagname: {
        'DerivedDiagnostic': firso_,
        'analyze': False,
        'requiredDiagnostics': firso_.requiredDiagnostics,
        'availableStatistics': firso_.availableStatistics,
        'label': firso_.label,
    },
    #NOTE: a failure results when 'analyze' is True under CR*, any one of the experiments
    # does not have CR* (i.e., a non-ensemble-DA experiment), and any one of the experiments
    # does have CR*
    #CRy*
    'CRyb': {
        'iter': 'bg',
        'analyze': False,
        'DerivedDiagnostic': bcry_,
        'requiredDiagnostics': bcry_.requiredDiagnostics,
        'availableStatistics': bcry_.availableStatistics,
        'label': bcry_.label,
    },
    'CRya': {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': acry_,
        'requiredDiagnostics': acry_.requiredDiagnostics,
        'availableStatistics': acry_.availableStatistics,
        'label': acry_.label,
    },
    'CRyf': {
        'DerivedDiagnostic': fcry_,
        'analyze': True,
        'requiredDiagnostics': fcry_.requiredDiagnostics,
        'availableStatistics': fcry_.availableStatistics,
        'label': fcry_.label,
    },
    # CRx*
    'CRxb': {
        'iter': 'bg',
        'analyze': True,
        'DerivedDiagnostic': bcrx_,
        'requiredDiagnostics': bcrx_.requiredDiagnostics,
        'availableStatistics': bcrx_.availableStatistics,
        'label': bcrx_.label,
    },
    'CRxa': {
        'iter': 'an',
        'analyze': False,
        'DerivedDiagnostic': acrx_,
        'requiredDiagnostics': acrx_.requiredDiagnostics,
        'availableStatistics': acrx_.availableStatistics,
        'label': acrx_.label,
    },
    'CRxf': {
        'DerivedDiagnostic': fcrx_,
        'analyze': False,
        'requiredDiagnostics': fcrx_.requiredDiagnostics,
        'availableStatistics': fcrx_.availableStatistics,
        'label': fcrx_.label,
    },
    'CRxinf': {
        'DerivedDiagnostic': icrx_,
        'analyze': False,
        'requiredDiagnostics': icrx_.requiredDiagnostics,
        'availableStatistics': icrx_.availableStatistics,
        'label': icrx_.label,
    },
    boeni_.diagname: {
        'iter': 'bg',
        'DerivedDiagnostic': boeni_,
        'analyze': False,
        'requiredDiagnostics': boeni_.requiredDiagnostics,
        'availableStatistics': boeni_.availableStatistics,
        'label': boeni_.label,
    },
    aoeni_.diagname: {
        'iter': 'an',
        'DerivedDiagnostic': aoeni_,
        'analyze': False,
        'requiredDiagnostics': aoeni_.requiredDiagnostics,
        'availableStatistics': aoeni_.availableStatistics,
        'label': aoeni_.label,
    },
    foeni_.diagname: {
        'DerivedDiagnostic': foeni_,
        'analyze': False,
        'requiredDiagnostics': foeni_.requiredDiagnostics,
        'availableStatistics': foeni_.availableStatistics,
        'label': foeni_.label,
    },
    innovratio_.diagname: {
        'DerivedDiagnostic': innovratio_,
        'analyze': False,
        'requiredDiagnostics': innovratio_.requiredDiagnostics,
        'availableStatistics': innovratio_.availableStatistics,
        'label': innovratio_.label,
    },
    sry_.diagname: {
        'DerivedDiagnostic': sry_,
        'analyze': False,
        'requiredDiagnostics': sry_.requiredDiagnostics,
        'availableStatistics': sry_.availableStatistics,
        'label': sry_.label,
    },
    srxeda_.diagname: {
        'DerivedDiagnostic': srxeda_,
        'analyze': True,
        'requiredDiagnostics': srxeda_.requiredDiagnostics,
        'availableStatistics': srxeda_.availableStatistics,
        'label': srxeda_.label,
    },
    srxrtpp_.diagname: {
        'DerivedDiagnostic': srxrtpp_,
        'analyze': True,
        'requiredDiagnostics': srxrtpp_.requiredDiagnostics,
        'availableStatistics': srxrtpp_.availableStatistics,
        'label': srxrtpp_.label,
    },
    approxdoadob_.diagname: {
        'DerivedDiagnostic': approxdoadob_,
        'analyze': False,
        'requiredDiagnostics': approxdoadob_.requiredDiagnostics,
        'availableStatistics': approxdoadob_.availableStatistics,
        'label': approxdoadob_.label,
    },
    approxrltvdoadob_.diagname: {
        'DerivedDiagnostic': approxrltvdoadob_,
        'analyze': False,
        'requiredDiagnostics': approxrltvdoadob_.requiredDiagnostics,
        'availableStatistics': approxrltvdoadob_.availableStatistics,
        'label': approxrltvdoadob_.label,
    },
    departurespread_.diagname: {
        'DerivedDiagnostic': departurespread_,
        'analyze': False,
        'requiredDiagnostics': departurespread_.requiredDiagnostics,
        'availableStatistics': departurespread_.availableStatistics,
        'label': r'$d$',
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
        'label': r'total spread',
    },
}

# classifications of diagnostics useful for plotting conventions
# TODO: incorporate these into availableDiagnostics dict

diagnosticIndependentStatistics = ['Count'] \
  + ObsSpaceConsistencyRatio.availableStatistics \
  + ModelSpaceConsistencyRatio.availableStatistics

statisticDependentDiagnostics = [
  biso_.diagname,
  aiso_.diagname,
  fiso_.diagname,
  boeni_.diagname,
  aoeni_.diagname,
  foeni_.diagname,
  innovratio_.diagname,
  sry_.diagname,
  srxeda_.diagname,
  srxrtpp_.diagname,
]

absoluteOnlyDiagnostics = set(statisticDependentDiagnostics + [
  'bc', 'obs', 'bak', 'ana', 'h(x)', 'amb',
  'CRyb', 'CRya', 'CRyf',
  'CRxb', 'CRxa', 'CRxinf',
  'doadob', 'dabdob', 'rltv_doadob',
  approxdoadob_.diagname,
  approxrltvdoadob_.diagname,
  birso_.diagname, airso_.diagname, firso_.diagname,
  #'sigmaoa', 'sigmaob', 'sigmaof',
])

#TODO: have this function return a list of diagnosticConfiguration or Diagnostic (new class) objects
#      instead of a list of dicts
def diagnosticConfigs(diagnosticNames_, DiagSpaceName, includeEnsembleDiagnostics=True,
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
        config['onlyDiagSpaces'] = config.get('onlyDiagSpaces',[])
        config['availableStatistics'] = config.get('availableStatistics', selectedStatistics)
        config['selectedStatistics'] = config.get('selectedStatistics', config['availableStatistics'])
        config['label'] = config.get('label', diagnosticName)
        config['requiredDiagnostics'] = set(config.get('requiredDiagnostics',[]))

        # diagnosticConfig is undefined for the following cases
        if (not config[vu.mean] and not config[vu.ensemble]): continue
        if (len(config['onlyDiagSpaces']) > 0 and
            DiagSpaceName not in config['onlyDiagSpaces']): continue

        config['dsName'] = DiagSpaceName
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
