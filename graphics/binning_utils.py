#!/usr/bin/env python3

from collections import defaultdict
from collections.abc import Iterable
from copy import deepcopy
import inspect
import itertools
import logging
import numpy as np
import os
import plot_utils as pu
from binning_params import allCIErrParams, ABEIParams
import var_utils as vu

_logger = logging.getLogger(__name__)

#===========================
# binning method descriptors
#===========================
# TODO(JJG): should be moved to binning_params

# identity
identityBinMethod = 'identity'

# none
noBinMethod = identityBinMethod

#QC
passQCMethod = 'pass'
goodQCMethod = 'good'
badQCMethod = 'bad'
allQCMethod = 'all'

#jet-stream pressure
P_jet_min = 250.0
P_jet_max = 350.0
P_jet_val = '{:.0f}'.format(0.5 * (P_jet_min + P_jet_max))
PjetMethod = 'P='+P_jet_val+'hPa'

#jet-stream altitude
alt_jet_min = 9500.0
alt_jet_max = 10500.0
alt_jet_val = '{:.0f}'.format(0.5 * (alt_jet_min + alt_jet_max))
altjetMethod = 'alt='+alt_jet_val+'m'
impactjetMethod = 'alt='+alt_jet_val+'m'

#LocalHour
LH0  = 0.0
LH1  = 23.0
LHDT = 1.0

#named latitude-bands
latbandsMethod = 'LatBands'
troplatbandsMethod = 'TropicalLatBands'
polarlatbandsMethod = 'PolarLatBands'

#named surface-type-bands
seasurfMethod = 'sea'
landsurfMethod = 'land'
mixsurfMethod = 'mixed-land-sea'
allsurfMethod = 'all-surface'
surfbandsMethod = 'surface-type'

#named cloudiness-bands
clrskyMethod = 'clear'
cldskyMethod = 'cloudy'
mixskyMethod = 'mixed-clrcld'
allskyMethod = 'allsky'
cloudbandsMethod = 'cloudiness'
allskyBTMethod = vu.obsVarBT+'_'+allskyMethod
#clrskyBTMethod = vu.obsVarBT+'_'+clrskyMethod

clrskyCFThreshold = 0.05
cldskyCFThreshold = 1.0 - clrskyCFThreshold

ClearCloudModeMethod = 'ClearCloudMode'

# cloud impact (CI)
OkamotoMethod       = 'Okamoto'
ModHarnischMethod      = 'ModHarnisch'
MZ19Method = 'MZ19'
QuadratureMethod = 'Quadrature'
MCIMethod = 'MCI'
OCIMethod = 'OCI'
CFQuadratureMethod = 'CFQuadrature'

# geostationary IR instrument longitude parameters
geoirlatlonboxMethod = 'geoirBox'
geoirMaxZenith = 65.0

geoirlatlonBoxParams = defaultdict(list)

# activate desired instruments by uncommenting them below

#seviri_m11 = 'seviri_m11'
#geoirlatlonBoxParams['values'] += [seviri_m11]
#geoirlatlonBoxParams['centerLon'] += [0.]

#seviri_m08 = 'seviri_m08'
#geoirlatlonBoxParams['values'] += [seviri_m08]
#geoirlatlonBoxParams['centerLon'] += [41.5]

ahi_himawari8 = 'ahi_himawari8'
geoirlatlonBoxParams['values'] += [ahi_himawari8]
geoirlatlonBoxParams['centerLon'] += [140.7]

#abi_g17 = 'abi_g17'
#geoirlatlonBoxParams['values'] += [abi_g17]
#geoirlatlonBoxParams['centerLon'] += [360. - 137.2]

abi_g16 = 'abi_g16'
geoirlatlonBoxParams['values'] += [abi_g16]
geoirlatlonBoxParams['centerLon'] += [360. - 75.2]

# glint angle
maxGlint = 90.0

#========================
# binning where functions
#========================

# Note: NaN and Inf values have mask set to true by default

def equalBound(x, bound, missingValue=True):
    if isinstance(bound, str):
        mask = np.char.equal(x, bound)
    else:
        mask = np.full_like(x, missingValue, bool)
        finite = np.isfinite(x)
        mask[finite] = np.equal(x[finite], bound)
        #mask[~finite] = missingValue
    return mask

def notEqualBound(x, bound, missingValue=True):
    if isinstance(bound, str):
        mask = np.char.not_equal(x, bound)
    else:
        mask = np.full_like(x, missingValue, bool)
        finite = np.isfinite(x)
        mask[finite] = np.not_equal(x[finite], bound)
    return mask

def notEqualAnyBound(x, bounds, missingValue=True):
    assert isinstance(bounds, Iterable), \
        ('ERROR, bounds must be Iterable for notEqualAnyBound')

    mask = np.full_like(x, True, bool)
    if isinstance(bounds[0], str):
        for bound in bounds:
            mask = np.logical_and(mask,
                     np.char.not_equal(x, bound))
    else:
        finite = np.isfinite(x)
        for bound in bounds:
            mask[finite] = np.logical_and(mask[finite],
                           np.not_equal(x[finite], bound))
        mask[~finite] = missingValue
    return mask

def lessEqualBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.full_like(finite, missingValue, bool)
    mask[finite] = np.less_equal(x[finite], bound)
    return mask

def lessBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.full_like(finite, missingValue, bool)
    mask[finite] = np.less(x[finite], bound)
    return mask

def greatEqualBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.full_like(finite, missingValue, bool)
    mask[finite] = np.greater_equal(x[finite], bound)
    return mask

def greatBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.full_like(finite, missingValue, bool)
    mask[finite] = np.greater(x[finite], bound)
    return mask

class multiBounds():
  def __init__(self, x, bounds, missingValue=True):
    bound1 = bounds[0]
    bound2 = bounds[1]
    self.__finite = np.isfinite(x)
    self.__mask = np.full_like(self.__finite, missingValue, bool)

    below = lessBound(x[self.__finite], bound1)
    above = greatBound(x[self.__finite], bound2)
    self.__outside = np.logical_or(below, above)

    belowEqual = lessEqualBound(x[self.__finite], bound1)
    aboveEqual = greatEqualBound(x[self.__finite], bound2)
    self.__outsideEqual = np.logical_or(belowEqual, aboveEqual)

  def outside(self):
    self.__mask[self.__finite] = self.__outside
    return self.__mask

  def outsideEqual(self):
    self.__mask[self.__finite] = self.__outsideEqual
    return self.__mask

  def inside(self):
    self.__mask[self.__finite] = np.logical_not(self.__outsideEqual)
    return self.__mask

  def insideEqual(self):
    self.__mask[self.__finite] = np.logical_not(self.__outside)
    return self.__mask


def insideBounds(x, bounds=[], missingValue=True):
    test =  multiBounds(x, bounds, missingValue)
    return test.inside()
#    bound1 = bounds[0]
#    bound2 = bounds[1]
#    finite = np.isfinite(x)
#    mask = np.empty_like(finite, bool)
#    belowbounds = lessEqualBound(x[finite], bound1)
#    abovebounds = greatEqualBound(x[finite], bound2)
#    mask[finite] = np.logical_not(
#        np.logical_or(belowbounds, abovebounds ))
#    mask[~finite] = missingValue
#    return mask

def insideEqualBounds(x, bounds=[], missingValue=True):
    test =  multiBounds(x, bounds, missingValue)
    return test.insideEqual()
#    bound1 = bounds[0]
#    bound2 = bounds[1]
#    finite = np.isfinite(x)
#    mask = np.empty_like(finite, bool)
#    belowbounds = lessBound(x[finite], bound1)
#    abovebounds = greatBound(x[finite], bound2)
#    mask[finite] = np.logical_not(
#        np.logical_or(belowbounds, abovebounds))
#    mask[~finite] = missingValue
#    return mask

def outsideBounds(x, bounds=[], missingValue=True):
    test =  multiBounds(x, bounds, missingValue)
    return test.outside()
#    bound1 = bounds[0]
#    bound2 = bounds[1]
#    finite = np.isfinite(x)
#    mask = np.empty_like(finite, bool)
#    belowbounds = lessBound(x[finite], bound1)
#    abovebounds = greatBound(x[finite], bound2)
#    mask[finite] = np.logical_or(
#        belowbounds, abovebounds )
#    mask[~finite] = missingValue
#    return mask

def outsideEqualBounds(x, bounds=[], missingValue=True):
    test =  multiBounds(x, bounds, missingValue)
    return test.outsideEqual()
#    bound1 = bounds[0]
#    bound2 = bounds[1]
#    finite = np.isfinite(x)
#    mask = np.empty_like(finite, bool)
#    belowbounds = lessEqualBound(x[finite], bound1)
#    abovebounds = greatEqualBound(x[finite], bound2)
#    mask[finite] = np.logical_or(
#        belowbounds, abovebounds )
#    mask[~finite] = missingValue
#    return mask

#=========================================================
# BinFunction classes to be accessed w/ BinFunctionWrapper
# i.e., functions of variables contained in the database
#=========================================================

################
## Base Classes
################

class BaseLocFunction:
  '''base class for evaluating location-specific quantities'''
  def __init__(self):
    self.baseVars = []
    self._initBaseVars()

    for baseVar in self.baseVars:
      if isinstance(baseVar, BaseLocFunction):
        for vv in baseVar.baseVars:
          if vv not in self.baseVars:
            self.baseVars.append(vv)
    self.baseVars = pu.uniqueMembers(self.baseVars)
    self._variables = {}

  def _initBaseVars(self):
    '''virtual method, used to populate the required self.baseVars'''
    raise NotImplementedError()

  def __calculateBaseVars(self, dbVals, insituParameters):
    '''populate self._variables with all baseVar values'''
    # store shallow copies of raw/independent variables
    for baseVar in self.baseVars:
      if isinstance(baseVar, str):
        self._variables[baseVar] = dbVals[insituParameters[baseVar]]

    # evaluate dependent variables
    for baseVar in self.baseVars:
      if isinstance(baseVar, BaseLocFunction):
        baseVar.evaluate(dbVals, insituParameters)
        # TODO: avoid many copies of the same derived variables
        self._variables[baseVar.__class__.__name__] = baseVar._get()

  def evaluate(self, dbVals, insituParameters):
    coords = self._coords(insituParameters)
    if coords not in self._variables:
      self.__calculateBaseVars(dbVals, insituParameters)
      # TODO: avoid many copies of the same derived variables, possibly store in dbVals
      self._variables[coords] = self._get()

    return self._variables[coords]

  def _coords(self, insituParameters):
    '''
    virtual method, should return tuple of insituParameters over which this function is uniform
    '''
    raise NotImplementedError()

  def _get(self):
    '''
    virtual method, should return 1D location-resolved array of function values
    '''
    raise NotImplementedError()


class UniformLocFunction(BaseLocFunction):
  '''LocFunction that is uniform for each dsName across all other insituParameters'''
  def _coords(self, insituParameters):
    return (self.__class__.__name__, insituParameters['dsName'])


class InsituLocFunction(BaseLocFunction):
  '''LocFunction that depends on insituParameters'''
  def _coords(self, insituParameters):
    c = list(insituParameters['coords'])
    c.append(self.__class__.__name__) #only necessary if storing in dbVals
    return tuple(c)


###################
## Derived Classes
###################

class GlintAngle(UniformLocFunction):
    '''satellite sun glint angle'''
    def _initBaseVars(self):
        self.baseVars.append(vu.senzenMeta)
        self.baseVars.append(vu.senaziMeta)
        self.baseVars.append(vu.solzenMeta)
        self.baseVars.append(vu.solaziMeta)

    def _get(self):
        senazi = self._variables[vu.senaziMeta]
        solazi = self._variables[vu.solaziMeta]

        relazi = np.abs(np.subtract(solazi, senazi))
        p = greatBound(relazi, 180.0, False)
        relazi[p] = np.subtract(360.0, relazi[p])
        relazi = np.multiply(np.subtract(180.0, relazi), vu.deg2rad)

        senzen = np.multiply(self._variables[vu.senzenMeta], vu.deg2rad)
        solzen = np.multiply(self._variables[vu.solzenMeta], vu.deg2rad)

        glint = np.add(np.multiply(np.cos(solzen), np.cos(senzen)),
                    np.multiply(np.sin(solzen),
                        np.multiply(np.sin(senzen), np.cos(relazi))))

        glint[greatBound(glint, 1.0)] = np.NaN
        glint[lessBound(glint, -1.0)] = np.NaN

        glint = np.multiply(np.arccos(glint), vu.rad2deg)
        glint[greatBound(glint, maxGlint, False)] = maxGlint

        return glint


class LocalHour(UniformLocFunction):
    '''local hour of the day based on longitude offset from 0'''
    def _initBaseVars(self):
        self.baseVars.append(vu.datetimeMeta)
        self.baseVars.append(vu.lonMeta)

    def _get(self):
        TimeStr = self._variables[vu.datetimeMeta]
        tzOffset = np.divide(self._variables[vu.lonMeta], 15.0)

        hh = np.empty_like(tzOffset, dtype=np.float32)
        mmi = np.empty_like(tzOffset, dtype=np.float32)
        ss = np.empty_like(tzOffset, dtype=np.float32)

        t0 = LH0
        t1 = LH1
        dt = LHDT
        for ii, Time in enumerate(TimeStr):
            ## Expecting Time to fit YYYY-MM-DDThh:mm:ssZ
            # YYYY[ii] = float(Time[0:4])
            # MMo[ii]  = float(Time[5:7])
            # DD[ii]   = float(Time[8:10])
            hh[ii]   = float(Time[11:13])
            mmi[ii]  = float(Time[14:16]) / 60.0
            ss[ii]   = float(Time[17:19]) / 3600.0

        LH = hh + mmi + ss + tzOffset

        yesterday = (LH < t0 - 0.5*dt)
        LH[yesterday] = LH[yesterday] + 24.0

        LH = np.mod(LH, 24.0)

        tomorrow = (LH >= t1 + 0.5*dt)
        LH[tomorrow] = LH[tomorrow] - 24.0

        return LH


## cloud impact parameteric functions

class MaximumCFx(UniformLocFunction):
    '''maximum cloud overlap based cloud fraction in model profile'''
    def _initBaseVars(self):
        self.baseVars.append(vu.cldfracGeo)

    def _get(self):
        cfAllLevels = self._variables[vu.cldfracGeo]

        # take maximum along second axis (levels)
        cfMax = np.max(cfAllLevels, 1)

        return cfMax


class CFQuadrature(UniformLocFunction):
    '''quadrature sum of simulated and retrieved cloud fractions'''
    def _initBaseVars(self):
        self.baseVars.append(MaximumCFx())
        self.baseVars.append(vu.cldfracMeta)

    def _get(self):
        CFx = self._variables['MaximumCFx']
        CFy = self._variables[vu.cldfracMeta]
        wx = 1.
        wy = 1.
        weights = np.asarray([wx, wy])
        return np.divide(
                 np.sqrt(
                   np.add(
                     np.square(weights[0] * CFx),
                     np.square(weights[1] * CFy)
                   )
                 ),
                 np.sqrt(weights.sum())
               )


class BTclr(InsituLocFunction):
    '''clear sky brightness temperature'''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def _get(self):
        BTbak = self._variables[vu.selfHofXValue]
        out = deepcopy(self._variables[vu.clrskyBTDiag])
        p = lessBound(out, 1.0, False)
        out[p] = BTbak[p]

        return out


class UnscaledOCI(InsituLocFunction):
    '''observed cloud impact'''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(BTclr())

    def _get(self):
        BTobs = self._variables[vu.selfObsValue]
        BTclr = self._variables['BTclr']

        #Scale Co by retrieved cloud fraction
        out = np.subtract(BTclr, BTobs)

        return out


class OCI(InsituLocFunction):
    '''observed cloud impact scaled by cloud fraction'''
    def _initBaseVars(self):
        self.baseVars.append(vu.cldfracMeta)
        self.baseVars.append(UnscaledOCI())

    def _get(self):
        OCI = self._variables['OCI']
        CldFracY = self._variables[vu.cldfracMeta]

        #Scale Co by retrieved cloud fraction
        out = np.multiply(CldFracY, OCI)

        return out


class MCI(InsituLocFunction):
    '''modeled cloud impact'''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(BTclr())

    def _get(self):
        BTbak = self._variables[vu.selfHofXValue]
        BTclr = self._variables['BTclr']

        #Scale only Co by retrieved cloud fraction
        out = np.subtract(BTclr, BTbak)

        return out


class SqrtMCI(InsituLocFunction):
    '''sqrt of modeled cloud impact'''
    def _initBaseVars(self):
        self.baseVars = [MCI()]

    def _get(self):
        MCI = self._variables['MCI']
        return np.sqrt(np.abs(MCI))


## functions related to the Asymmetric Cloud Impact (ACI)
class ACIMZ19(InsituLocFunction):
    '''Minamide and Zhang (2019) ACI'''
    def _initBaseVars(self):
        self.baseVars.append(UnscaledOCI())
        self.baseVars.append(MCI())

    def _get(self):
        MCI = self._variables['MCI']
        OCI = self._variables['UnscaledOCI']
        # Minamide and Zhang (2018)
        ACI = np.subtract(np.abs(OCI), np.abs(MCI))

        return ACI


class ACIQuadrature(InsituLocFunction):
    '''quadrature-based ACI'''
    def _initBaseVars(self):
        self.baseVars.append(OCI())
        self.baseVars.append(MCI())

    def _get(self):
        MCI = np.square(self._variables['MCI'])
        OCI = np.square(self._variables['OCI'])

        # MCI = np.square(np.maximum(MCI, 0.))

        d = OCI - MCI

        # subtract OCI and MCI in quadrature, keeping sign outside sqrt
        ACI = np.full_like(d, np.NaN)
        p = np.isfinite(d)
        ACI[p] = np.multiply(np.sqrt(np.abs(d[p])), np.sign(d[p]))

        return ACI


class ABEILambda:
    '''Minamide and Zhang (2019) ACI-parameterized inflation factor'''
    minLambda = 1.0
    maxLambda = 1.4
    def __init__(self):
        self.ACI = ACIMZ19()
        self.baseVars = []
        self.baseVars = pu.uniqueMembers(self.baseVars + self.ACI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        dsName = insituParameters['dsName']
        if dsName is None or dsName not in ABEIParams:
            _logger.error('dsName not available in ABEIParams => '+dsName)
            os._exit(1)

        # varName, ch = vu.splitIntSuffix(insituParameters[vu.selfDepValue])
        varName, ch = vu.splitIntSuffix(insituParameters[vu.selfHofXValue])
        LambdaOverACI = ABEIParams[dsName][(int(ch))]['LambdaOverACI']

        ACI = self.ACI.evaluate(dbVals, insituParameters)
        out = np.ones(ACI.shape)
        crit = (ACI > 0.0)
        out[crit] = np.multiply(ACI[crit], LambdaOverACI) + self.minLambda
        crit = (ACI >= (self.maxLambda - self.minLambda) / LambdaOverACI)
        out[crit] = self.maxLambda

        return out


## functions related to the Symmetric Cloud Impact (SCI)
class SCIOkamoto(InsituLocFunction):
    '''Okamoto et al. (2014) SCI'''
    def _initBaseVars(self):
        self.baseVars.append(UnscaledOCI())
        self.baseVars.append(MCI())

    def _get(self):
        MCI = self._variables['MCI']
        OCI = self._variables['UnscaledOCI']

        # Okamoto, et al. (2014)
        # Co = abs(Bias-Corrected BTobs - BTclr)
        # Cm = abs(BTbak - BTclr)
        SCI = np.multiply(0.5, np.add(np.abs(OCI), np.abs(MCI)))

        return SCI


class SCIQuadrature(InsituLocFunction):
    '''quadrature sum of MCI and cloud-fracion-scaled OCI'''
    def _initBaseVars(self):
        self.baseVars.append(OCI())
        self.baseVars.append(MCI())

    def _get(self):
        MCI = self._variables['MCI']
        OCI = self._variables['OCI']

        # add OCI and MCI in quadrature
        SCI = np.sqrt(np.square(OCI) + np.square(MCI))

        return SCI


class SCIModHarnisch(InsituLocFunction):
    '''Harnisch et al. (2016) SCI, using BTclr in place of BT threshold'''
    def _initBaseVars(self):
        self.baseVars.append(OCI())
        self.baseVars.append(MCI())

    def _get(self):
        MCI = self._variables['MCI']
        OCI = self._variables['OCI']

        # Modified Harnisch, et al.
        zeros = np.full_like(MCI,0.0)
        SCI = np.multiply( 0.5,
                 np.add(np.maximum(zeros, OCI),
                        np.maximum(zeros, MCI)) )
        return SCI


## functions related to errors and spread

class STDofHofX:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfHofXValue)
        #TODO: enable ensembleHofXValue to request ensemble HofX values
        #self.baseVars.append(vu.ensembleHofXValue)

    def evaluate(self, dbVals, insituParameters):
        meanVarName = insituParameters[vu.selfHofXValue]
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
            nLocs = len(dbVals[meanVarName])
            std = np.full(nLocs, np.NaN)
        return std


class TotalSpread:
    '''
    Calculates the total or some component of the spread, σ, where
    σ could be σ_o, σ_h, or σ_total = √(σ_o^2 + σ_h^2)
    The default behavior is to evaluate σ_total; however STDofHofX.evaluate returns all NaN
    values when ensemble HofX values are unavailable, which defaults back to only returning
    the ObsError.  Dependent functions may override the errortype in order to choose a particular
    error component for both ensemble and deterministic applications.
    '''
    def __init__(self, errortype='total'):
        self.baseVars = []
        self.baseVars.append(vu.selfErrorValue)
        self.ensembleSpread = STDofHofX()

        ## Does not technically work to add self.ensembleSpread.baseVars until DiagnoseObsStatistics
        ## explicitly handles binned ensemble variables.  Only works under the assumption
        ## that sigmab, simaa, or sigmaf diagnostics are already being calculated for
        ## ensemble-based applications.
        self.baseVars = pu.uniqueMembers(self.baseVars + self.ensembleSpread.baseVars)

        # select error type, can override in derived classes
        self.errortype = errortype
        self.errortypes = {
            'obs': self.getObsError,
            'bak': self.getEnsSpread,
            'total': self.getTotalSpread,
        }

    def getObsError(self, dbVals, insituParameters):
        err = dbVals[insituParameters[vu.selfErrorValue]]
        err[lessEqualBound(err, 0.0)] = np.NaN
        return err

    def getEnsSpread(self, dbVals, insituParameters):
        return self.ensembleSpread.evaluate(dbVals, insituParameters)

    def getTotalSpread(self, dbVals, insituParameters):
        sigmao = self.getObsError(dbVals, insituParameters)
        sigmah = self.getEnsSpread(dbVals, insituParameters)

        validO = greatBound(sigmao, 0., False)
        validH = greatBound(sigmah, 0., False)
        bothValid = np.logical_and(validO, validH)
        onlyO = np.logical_and(validO, np.logical_not(validH))
        onlyH = np.logical_and(validH, np.logical_not(validO))
        validvalues = np.full_like(sigmao, np.NaN)
        validvalues[onlyO] = sigmao[onlyO]
        validvalues[onlyH] = sigmah[onlyH]

        spread = np.sqrt(np.square(sigmao) + np.square(sigmah),
                         out=validvalues, where=bothValid)

        return spread

    def evaluate(self, dbVals, insituParameters):
        return self.errortypes[self.errortype](dbVals, insituParameters)


class CITotalSpread(TotalSpread):
    '''
    Derived class that overrides getObsError in TotalSpread base class
    with a CI-parameterized ObsError
    '''
    mpasFCRes = 30 # km [30, 120]
    biasCorrectType = {} #[None, 'constant', 'varbc']
    #biasCorrectType['abi_g16'] = 'constant'
    biasCorrectType['abi_g16'] = None
    biasCorrectType['ahi_himawari8'] = None

    def __init__(self, CIName, CIClass, CIVariable=vu.obsVarCI, errortype='total'):
        super().__init__(errortype)
        self.CI = CIClass()
        self.CIName = CIName
        self.CIVariable = CIVariable

        self.baseVars.append(vu.selfHofXValue)
        self.baseVars = pu.uniqueMembers(self.baseVars + self.CI.baseVars)

    def getObsError(self, dbVals, insituParameters):
        # Parameterize BTerr as a ramped step function
        # --------------------------------------------
        # STD1 |. . . ____
        #      |     /.
        #      |    / .
        #      |   /  .
        # STD0 |__/   .
        #      |  .   .
        #      |__.___.___
        #         .   .
        #        CI0 CI1
        #---------------------------------------------
        dsName = insituParameters['dsName']
        CIErrParams = deepcopy(allCIErrParams[
                           (self.CIVariable, self.mpasFCRes, self.biasCorrectType.get(dsName,None))
                       ])

        if dsName is None or dsName not in CIErrParams:
            _logger.error('dsName not available in CIErrParams => '+dsName)
            os._exit(1)

        varName, ch = vu.splitIntSuffix(insituParameters[vu.selfHofXValue])
        STD0 = CIErrParams[dsName][(int(ch), self.CIName)]['ERR'][0]
        STD1 = CIErrParams[dsName][(int(ch), self.CIName)]['ERR'][1]
        CI0  = CIErrParams[dsName][(int(ch), self.CIName)]['X'][0]
        CI1  = CIErrParams[dsName][(int(ch), self.CIName)]['X'][1]
        slope = (STD1 - STD0) / (CI1 - CI0)

        CI = self.CI.evaluate(dbVals, insituParameters)
        belowramp = lessEqualBound(CI, CI0, False)
        aboveramp = greatEqualBound(CI, CI1, False)
        onramp    = insideBounds(CI, [CI0, CI1], False)

        err = np.full_like(CI, np.NaN)
        err[belowramp] = STD0
        err[onramp]    = STD0 + slope * (CI[onramp] - CI0)
        err[aboveramp] = STD1

        return err


## innovation/departure functions

class Departure(InsituLocFunction):
    '''
    Departures -- y - h(x)
    '''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfHofXValue)

    def _get(self):
        obs = self._variables[vu.selfObsValue]
        bak = self._variables[vu.selfHofXValue]
        dep = np.subtract(obs, bak)

        return dep


class ClearSkyDeparture(InsituLocFunction):
    '''
    Departures -- y - h_clear(x)
    '''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(BTclr())

    def _get(self):
        obs = self._variables[vu.selfObsValue]
        BTclr = self._variables['BTclr']
        dep = np.subtract(obs, BTclr)

        return dep


class LogDepartureRatio(InsituLocFunction):
    '''
    Log of departure ratio -- log(y / h(x))
    '''
    def _initBaseVars(self):
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfHofXValue)

    def _get(self):
        obs = self._variables[vu.selfObsValue]
        bak = self._variables[vu.selfHofXValue]
        logdep = np.log(np.divide(obs, bak))

        return logdep


## functions for departures normalized by spread
class NormalizedDeparture:
    '''
    Error-normalized departures -- (y - h(x)) / σ
    Derived classes can override the TotalSpread behavior by re-constructing the
    self.TotalSpread object.
    '''
    TotalSpread = TotalSpread('total')
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfHofXValue)

        ## Does not fully work to add self.TotalSpread.baseVars until DiagnoseObsStatistics
        ## explicitly handles binned ensemble variables.  Only works under the assumption
        ## that sigmab, simaa, or sigmaf diagnostics are already being calculated for
        ## ensemble-based applications.
        self.baseVars = pu.uniqueMembers(self.baseVars + self.TotalSpread.baseVars)

    def evaluate(self, dbVals, insituParameters):
        err = self.TotalSpread.evaluate(dbVals, insituParameters)
        obs = dbVals[insituParameters[vu.selfObsValue]]
        bak = dbVals[insituParameters[vu.selfHofXValue]]
        dep = np.subtract(obs, bak)

        return np.divide(dep, err)


class OkamotoNormalizedDeparture(NormalizedDeparture):
    TotalSpread = CITotalSpread(OkamotoMethod, SCIOkamoto)


class QuadratureNormalizedDeparture(NormalizedDeparture):
    TotalSpread = CITotalSpread(QuadratureMethod, SCIQuadrature, vu.obsVarCI)


#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(dbVals, REGION_NAME):
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = dbVals['longitdue']
#    lats = dbVals['latitude']
#    nlocs = len(lons)
#    REGIONS = isinsideregions(lons, lats, shp)
#    return REGIONS


#=========================================
# generic wrappers for BinFunction classes
#=========================================
class BaseBinFunction:
    def __init__(self, baseVars):
        # only expose string baseVars to higher-level code
        self.baseVarsStr = [v for v in baseVars if isinstance(v, str)]

        #for baseVar in baseVars:
        #    if isinstance(baseVar, str):
        #        self.baseVarsStr.append(baseVar)
        pass

    def dbVars(self, varName, fileFormat, outerIters_):
        dbVars = []

        if (not isinstance(outerIters_, Iterable)
           or isinstance(outerIters_,str)):
            outerIters = [outerIters_]
        else:
            outerIters = outerIters_

        for baseVar in self.baseVarsStr:
            for outerIter in outerIters:
                dbVar = vu.base2dbVar(
                    baseVar, varName, fileFormat, outerIter)
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)


class IdBinFunction(BaseBinFunction):
    def __init__(self, variable):
        super().__init__([variable])

    def evaluate(self, dbVals, insituParameters):
        return dbVals[insituParameters[self.baseVarsStr[0]]]


class BinFunction(BaseBinFunction):
    def __init__(self, function):
        self.function = function()
        assert hasattr(self.function, 'baseVars'), \
            ("ERROR, function class must have the baseVars attribute:", function)
        super().__init__(self.function.baseVars)

    def evaluate(self, dbVals, insituParameters):
        return self.function.evaluate(dbVals, insituParameters)


class BinFunctionWrapper:
    def __init__(self, config):
        self.dsName = config['dsName']
        self.fileFormat = config['fileFormat']

        variable = config['variable']
        varIsString = isinstance(variable,str)
        varIsClass = inspect.isclass(variable)
        assert varIsString ^ varIsClass, \
            ("ERROR: 'variable' must either be str or class", config)

        if varIsString:
            self.function = IdBinFunction(variable)

        if varIsClass:
            self.function = BinFunction(variable)

    def dbVars(self, varName, outerIters):
        return self.function.dbVars(varName, self.fileFormat, outerIters)

    def evaluate(self, dbVals, varName, outerIter):
        # setup context-specific insituParameters for the evaluation
        insituParameters = {}
        for baseVar in self.function.baseVarsStr:
            insituParameters[baseVar] = vu.base2dbVar(
                    baseVar, varName, self.fileFormat, outerIter)
        insituParameters['dsName'] = self.dsName
        insituParameters['coords'] = (varName, self.fileFormat, outerIter, self.dsName)

        # evaluate the function
        self.result = self.function.evaluate(dbVals, insituParameters)


#========================
# generic binning classes
#========================
class BinningAxis:
    '''
    Backend support for an individual binning axis
    '''
    # The constructor takes a configuration dict with the following elements:
    # start - first bin central value
    # stop - final bin central value
    # values - array of bin central values; starts and stops will be
    #   interpolated/extrapolated
    # Either (starts and stops) or values MUST be provided
    # step - delta between bin central values (overrides nsteps)
    # nsteps - number of bins, can be used instead of step
    # f - format string for storing and plotting values
    # period - period of repetition for bins, starting at 0.
    # transforms - two transformation functions for nonlinear stepping intervals
    # transforms[0] transforms (start, stop, and step) or values before
    #   generating arrays
    # transforms[1] transforms intermediate starts, stops, and values to final arrays
    #     E.g., one could specify (start, stop, and step) or values in a
    #     natural log space, then specify transforms=[None, np.exp] in order to
    #     create bins in a linear space

    nsteps = 16

    def __init__(self, config):

        # parse config for first guesses of all parameters
        validConfig = (
          (('start' in config and 'stop' in config) and 'values' not in config)
          or (('values' in config or 'value' in config)
            and 'start' not in config and 'stop' not in config)
        )
        assert validConfig, 'BinningAxis.__init__: missing config component'

        start = config.get('start', None)
        stop = config.get('stop', None)
        values = config.get('values', None)
        if values is None:
          value = config.get('value', None)
          if value is not None:
            values = [value]

        step = config.get('step', None)
        nsteps = config.get('nsteps', self.nsteps)
        f = config.get('format', '{:.0f}')
        period = config.get('period', None)
        transforms = config.get('transforms', [self.identityFunc, self.identityFunc])

        starts = []
        stops = []
        if values is None:
            # determine starts, stops, values
            values = []
            if period is not None and (stop<start or stop>=period):
                assert start >= 0.0, 'BinningAxis.__init__: can only handle positive periodic bin values'
                assert start < period, 'BinningAxis.__init__: can only handle start < period'
                if np.abs(stop - start) >= period:
                    range_ = period
                    if step is None:
                        step = range_ / float(nsteps)

                    start = 0.5*period + 0.5*step
                    stop = 0.5*period
                else:
                    while stop >= period: stop-=period
                    range_ = (stop + (period - start))
                    if step is None:
                        step = range_ / float(nsteps)

                assert step > 0.0, 'BinningAxis.__init__: can only handle positive steps for periodic bins'

                # handle regions above and below period independently
                allBinBounds = []
                allBinBounds.append(np.flip(np.arange(period, start-0.5*step, -step)))
                allBinBounds.append(np.arange(0.0, stop+1.5*step, step))
                for binBounds in allBinBounds:
                    for ibin in range(len(binBounds)-1):
                        starts.append(binBounds[ibin])
                        stops.append(binBounds[ibin+1])

                values = [0.5*(i+e) for i, e in zip(starts, stops)]

            else:
                range_ = stop - start
                if step is None:
                    step = range_ / float(nsteps)

                nsteps = np.trunc((stop - start) / step)
                start = transforms[0](start)
                stop = transforms[0](stop)
                step = (stop - start) / nsteps

                if nsteps == 1:
                    values = np.array([(stop+start)*0.5])
                    starts = np.array([start])
                    stops = np.array([stop])
                else:
                    values = np.arange(start, stop+step, step)
                    starts = np.subtract(values, 0.5*step)
                    stops = np.add(values, 0.5*step)

                values = transforms[1](values)
                starts = transforms[1](starts)
                stops = transforms[1](stops)
        elif len(values) == 1:
            starts = deepcopy(values)
            stops = deepcopy(values)
        else:
            assert self.isMonotonic(values), 'BinningAxis.__init__: can only handle monotonically increasing or decreasing values'
            ## initial transform
            values = transforms[0](values)

            ## determine starts and stops
            # extraploate first bound
            bound = values[0] - 0.5*(values[1] - values[0])

            # interpolate inner bounds
            for ii in range(len(values)-1):
              starts.append(bound)
              bound = values[ii] + 0.5*(values[ii+1] - values[ii])
              stops.append(bound)
            starts.append(bound)

            # extrapolate final bound
            bound = values[-1] + 0.5*(values[-1] - values[-2])
            stops.append(bound)

            ## final transform
            values = transforms[1](values)
            starts = transforms[1](starts)
            stops = transforms[1](stops)

        values = [f.format(v) for v in values]

        self.starts = np.asarray(starts)
        self.stops = np.asarray(stops)
        self.values = np.asarray(values)

        self.n = len(starts)
        self.index = np.arange(self.n)

    @staticmethod
    def identityFunc(x):
        return x

    # Check if given array is Monotonic
    @staticmethod
    def isMonotonic(A):
        return (all(A[i] <= A[i + 1] for i in range(len(A) - 1)) or
                all(A[i] >= A[i + 1] for i in range(len(A) - 1)))

    def indices(self):
        return self.index

    def resetIndex(self, index):
        self.index = np.asarray(index)

    def getstarts(self):
        return self.starts[self.index]

    def getstops(self):
        return self.stops[self.index]

    def getcentrals(self, astype=float):
        return self.values[self.index].astype(astype)

    def getvalues(self):
        return self.values[self.index]

class BinningAxes:
    '''
    Provides an interface to 1 or more BinningAxis objects for
    the purposes of univariate and multi-dimensional binning
    '''
    # constructor takes a dictionary object with variable:config key/value pairs:
    # axes = {
    #   (var1, 0): {'start': start1, 'stop': stop1, 'step': step1, 'format': fmt1}
    #   (var2, 1): {'start': start2, 'stop': stop2, 'step': step2, 'format': fmt2}
    #   etc...
    # }
    # The first component of each dictionary key is the variable name
    # The second component of each dictionary key is the order of that variable among
    # all axes. When plotting, the axes will be accessed in that order.
    # E.g., for 2D BinningAxes, the order==0 will be on the x-axis, and order==1 will
    # be on the y-axis
    # Each axes dictionary value is a configuration used to construct a BinningAxis object
    def __init__(self, axes):
        assert isinstance(axes, dict), 'BinningAxes.__init__: axes must be a dict'

        self.axes = {}
        self.naxes = len(axes)
        self.order = {}
        self.variables = [None]*self.naxes
        coordinates = [None]*self.naxes
        for (var, order), ax in axes.items():
            # initialize axis
            self.axes[var] = BinningAxis(ax)

            # initialize order
            assert order not in self.order.values(), 'BinningAxes.__init__: order must be unique'
            assert order in np.arange(self.naxes), 'BinningAxes.__init__: order must be sequential'
            self.order[var] = order
            self.variables[order] = var

            # initialize ordered coordinates
            coordinates[order] = self.axes[var].indices()

        if self.naxes == 1:
            coordGenerator = itertools.product(
                coordinates[0])
        elif self.naxes == 2:
            coordGenerator = itertools.product(
                coordinates[0],
                coordinates[1])
        elif self.naxes == 3:
            coordGenerator = itertools.product(
                coordinates[0],
                coordinates[1],
                coordinates[2])
        elif self.naxes == 4:
            coordGenerator = itertools.product(
                coordinates[0],
                coordinates[1],
                coordinates[2],
                coordinates[3])

        axesIndices = defaultdict(list)
        for coord in coordGenerator:
            for var, c in zip(self.variables, coord):
                axesIndices[var].append(c)

        for var, axisIndices in axesIndices.items():
            self.axes[var].resetIndex(axisIndices)

    def variable(self, var=None, order=None):
        # returns a variable associated with this BinningAxes
        # if self.naxes == 1, returns singular variable
        # elif var is not None, checks that var is one of the axes keys and returns it back
        # elif var is None and order is not None, returns self.variables[order]
        # else causes an error
        var_ = var
        if var_ is None:
            if self.naxes==1:
                var_ = list(self.axes.keys())[0]
            elif order is None or order not in np.arange(self.naxes):
                _logger.error('BinningAxes.variable: need to provide valid var or order for multi-variate axes')
                os._exit(1)
            else:
                var_ = self.variables[order]
        elif var_ not in self.axes:
            _logger.error('BinningAxes.variable: var not in axes: '+var_)
            os._exit(1)

        return var_

    def starts(self, var=None):
        return self.axes[self.variable(var)].getstarts()

    def stops(self, var=None):
        return self.axes[self.variable(var)].getstops()

    def centrals(self, var=None, astype=float):
        return self.axes[self.variable(var)].getcentrals(astype)

#    def valuesBYvar(self, var=None):
#        return self.axes[self.variable(var)].getvalues()

#    def index(self, index)
#        assert (index >= 0 and index < self.n), ('CompoundFilter.index: invalid index: ', index)
#
#        index = np.empty(self.naxes)
#        for var, axis in self.axes.items():
#            index[self.order[var]] = axis.index(index)
#        return tuple(index)

#    def indices(self):
#        indices = [None]*self.naxes
#        for var, axis in self.axes.items():
#            indices[self.order[var]] = axis.index
#
#        if self.naxes == 1:
#            return list(zip(indices[0]))
#        elif self.naxes == 2:
#            return list(zip(indices[0], indices[1]))
#        elif self.naxes == 3:
#            return list(zip(indices[0], indices[1], indices[2]))
#        elif self.naxes == 4:
#            return list(zip(indices[0], indices[1], indices[2], indices[3]))

    def values(self):
        values = [None]*self.naxes
        for var, axis in self.axes.items():
            values[self.order[var]] = axis.getvalues()

        if self.naxes == 1:
            return values[0]
        elif self.naxes == 2:
            coordTuples = list(zip(values[0], values[1]))
        elif self.naxes == 3:
            coordTuples = list(zip(values[0], values[1], values[2]))
        elif self.naxes == 4:
            coordTuples = list(zip(values[0], values[1], values[2], values[3]))

        #coords = np.asarray([','.join(list(c)) for c in coordTuples])
        coords = [','.join(list(c)) for c in coordTuples]

        return coords

#    def value(self, index):
#        value = [None]*self.naxes
#        for var, axis in self.axes.items():
#            ii = self.order[var]
#            value[ii] = axis.value(index)
#
#        return tuple(value)


class BinFilter:
    def __init__(self, config):
        self._where  = config['where']
        tmp         = config['bounds']
        nBins       = config['nBins']

        # note: for where functions that take multiple bounds (e.g., notEqualAnyBound)
        #   config['bounds'] can wrap an inner Iterable, e.g.,
        #   * config['bounds'] = [[0,1]] would apply the bounds 0 and 1 to all nBins bins
        #   * config['bounds'] = [[0,1], [5,6]] would apply the bounds 0 and 1 to the first bin,
        #                                                   the bounds 5 and 6 to the second bin,
        #                                                   and nBins must be equal to 2

        ibins = list(range(nBins))

        # allow for scalar and Iterable bounds
        if (not isinstance(tmp, Iterable) or
            isinstance(tmp, str)):
            self._bounds = np.empty(nBins, dtype=type(tmp))
            for ii in ibins:
                self._bounds[ii] = tmp
        else:
            self._bounds = np.empty(nBins, dtype=type(tmp[0]))

            # single element Iterable is applied uniformly to all bins
            if len(tmp) == 1:
                for ii in ibins:
                    self._bounds[ii] = tmp[0]
            # multiple element Iterable must be same length as nBins
            elif len(tmp) == nBins:
                for ii in ibins:
                    self._bounds[ii] = tmp[ii]
            else:
                _logger.error("'bounds' must be a scalar, single-member Iterable, or an Iterable with the same length as 'values'!")
                os._exit(1)
        self._function = BinFunctionWrapper(config)
        self.except_diags = config.get('except_diags', [])
        self.include_diags = config.get('include_diags', [])

        self.mask_value = config.get('mask_value', np.NaN)
        #TODO: add other actions besides mask_value/exclude

#    def baseVars(self):
#        return pu.uniqueMembers(self._function.baseVars)

    def dbVars(self, varName, outerIters):
        dbVars = self._function.dbVars(
            varName, outerIters)
        return pu.uniqueMembers(dbVars)

    def evaluate(self, dbVals, varName, outerIter):
        self._function.evaluate(dbVals, varName, outerIter)

    def updateMask(self, maskIn, diagName, ibin, maskValue=True):
        maskOut = deepcopy(maskIn)
        if diagName not in self.except_diags:
            # modify maskOut where mask is True
            mask = self._where(self._function.result,(self._bounds)[ibin])
            ashape = maskOut.shape
            mshape = mask.shape
            if mshape == ashape:
                maskOut[mask] = maskValue
            elif (len(ashape)==2 and len(mshape)==1):
                if ashape[0]==mshape[0]:
                    maskOut[mask,:] = maskValue
                elif ashape[1]==mshape[0]:
                    maskOut[:,mask] = maskValue
                else:
                    message = 'BinFilter.apply: shape mismatch in mask: '+str(mshape)+','+str(ashape)
                    _logger.error(message)
                    os._exit(1)
            elif len(mshape)==1 and mshape[0]==1:
                if mask[0]:
                    maskOut[:] = maskValue
            else:
                message = 'BinFilter.apply: mask is incorrectly defined: '+str(mshape)+','+str(ashape)
                _logger.error(message)
                os._exit(1)

        return maskOut

#    def maskArray(self, arrayLike, diagName, ibin):
#        array = np.asarray(deepcopy(arrayLike))
#
#        if diagName not in self.except_diags:
#            # remove locations where the mask is True
#            mask = self._where(self._function.result,(self._bounds)[ibin])
#            ashape = array.shape
#            mshape = mask.shape
#            if mshape == ashape:
#                array[mask] = self.mask_value
#            elif (len(ashape)==2 and len(mshape)==1):
#                if ashape[0]==mshape[0]:
#                    array[mask,:] = self.mask_value
#                elif ashape[1]==mshape[0]:
#                    array[:,mask] = self.mask_value
#                else:
#                    message = 'BinFilter.apply: shape mismatch in mask: '+str(mshape)+','+str(ashape)
#                    _logger.error(message)
#                    os._exit(1)
#            elif len(mshape)==1 and mshape[0]==1:
#                if mask[0]:
#                    array[:] = self.mask_value
#            else:
#                message = 'BinFilter.apply: mask is incorrectly defined: '+str(mshape)+','+str(ashape)
#                _logger.error(message)
#                os._exit(1)
#
#        return array


class BinMethod:
    ## exclusiveDiags
    # list of diagnostics to skip for a particular binMethod, unless they appear
    # in override_exclusiveDiags
    exclusiveDiags = [
        'obs','bak','ana',
        'SCI-'+OkamotoMethod,
        'ACI-'+MZ19Method,
        'MCI',
        'CFy',
    ]
    ## commonDiags
    # subset of exclusiveDiags that are commonly desired in plots
    commonDiags = [
        'obs','bak','ana',
    ]
    def __init__(self, config):
        #allows for scalar, str, and Iterable 'values'
        tmp = config['values']
        self.__values = []
        if (not isinstance(tmp, Iterable) or
            isinstance(tmp, str)):
            self.__values += [tmp]
        else:
            self.__values = list(tmp)

        # list of diagnostics to exclude
        self._excludeDiags = deepcopy(self.exclusiveDiags)
        override = config.get('override_exclusiveDiags',[])
        for diag in override:
            if diag in self._excludeDiags:
                self._excludeDiags.remove(diag)

        # list of variables to include
        self._includeVars = config.get('include variables', None)
        assert self._includeVars is None or isinstance(self._includeVars, Iterable), (
          self.__class__.__name__+'.__init__: invalid self._includeVars: '+str(self._includeVars))

        # list of variables to exclude
        self._excludeVars = config.get('exclude variables', None)
        assert self._excludeVars is None or isinstance(self._excludeVars, Iterable), (
          self.__class__.__name__+'.__init__: invalid self._excludeVars: '+str(self._excludeVars))

        fconf = {}
        fconf['dsName'] = config['dsName']
        fconf['fileFormat'] = config['fileFormat']
        fconf['nBins'] = len(self.__values)

        self.filters = []
        for filterConf in config['filters']:
            filterConf.update(fconf)
            self.filters.append(BinFilter(filterConf))

        enoughBounds = False
        for Filter in self.filters:
            if len(Filter._bounds) == len(self.__values):
                enoughBounds = True
        assert enoughBounds, '\n\nERROR: BinMethod : at least one filter must have len(bounds) == len(values)!'

#    def baseVars(self):
#        baseVars = []
#        for Filter in self.filters:
#            for variable in Filter.baseVars():
#                baseVars.append(variable)
#        return pu.uniqueMembers(baseVars)

    def dbVars(self, varName, outerIters=None):
        dbVars = []
        for Filter in self.filters:
            dbVars += Filter.dbVars(
                varName, outerIters)
        return pu.uniqueMembers(dbVars)

    def evaluate(self, dbVals, varName, outerIter=None):
        for ii in list(range(len(self.filters))):
            self.filters[ii].evaluate(
                dbVals, varName, outerIter)

    def apply(self, array, diagName, binVal):
        ibin = self.__values.index(binVal)
        mask = np.full_like(array, False, bool)
        for Filter in self.filters:
            mask = Filter.updateMask(
                mask, diagName, ibin, maskValue=True)
        masked_array = np.asarray(deepcopy(array))
        masked_array[mask] = np.NaN

#        masked_array = np.asarray(deepcopy(array))
#        for Filter in self.filters:
#            masked_array = Filter.maskArray(
#                masked_array, diagName, ibin)

        return masked_array

    def getvalues(self):
        return deepcopy(self.__values)

    def excludeDiag(self, diagName):
        return diagName in self._excludeDiags

    def excludeVariable(self, varName):
        exclude = False

        if self._includeVars is not None:
          exclude = varName not in self._includeVars

        if self._excludeVars is not None:
          exclude = exclude or (varName in self._excludeVars)

        return exclude
