#!/usr/bin/env python3

from collections import defaultdict
from collections.abc import Iterable
from copy import deepcopy
import inspect
import logging
import numpy as np
import os
import plot_utils as pu
from binning_params import allSCIErrParams, ABEIParams
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

#LocalHour
LH0  = 0.0
LH1  = 23.0
LHDT = 1.0

#named latitude-bands
latbandsMethod = 'LatBands'

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

# symmetric cloud impact (SCI)
OkamotoMethod       = 'Okamoto'
ScaleOkamotoMethod  = 'ScaledOkamoto'
ModHarnischMethod      = 'ModHarnisch'
ScaleModHarnischMethod = 'ScaledModHarnisch'

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
    mask = np.empty_like(x, bool)
    if isinstance(bound, str):
        mask = np.char.equal(x, bound)
    else:
        finite = np.isfinite(x)
        mask[finite] = np.equal(x[finite], bound)
        mask[~finite] = missingValue
    return mask

def notEqualBound(x, bound, missingValue=True):
    mask = np.empty_like(x, bool)
    if isinstance(bound, str):
        mask = np.char.not_equal(x, bound)
    else:
        finite = np.isfinite(x)
        mask[finite] = np.not_equal(x[finite], bound)
        mask[~finite] = missingValue
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
    mask = np.empty_like(finite, bool)
    mask[finite] = np.less_equal(x[finite], bound)
    mask[~finite] = missingValue
    return mask

def lessBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.empty_like(finite, bool)
    mask[finite] = np.less(x[finite], bound)
    mask[~finite] = missingValue
    return mask

def greatEqualBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.empty_like(finite, bool)
    mask[finite] = np.greater_equal(x[finite], bound)
    mask[~finite] = missingValue
    return mask

def greatBound(x, bound, missingValue=True):
    finite = np.isfinite(x)
    mask = np.empty_like(finite, bool)
    mask[finite] = np.greater(x[finite], bound)
    mask[~finite] = missingValue
    return mask

def betweenBounds(x, bound1, bound2, missingValue=True):
    finite = np.isfinite(x)
    belowbounds = lessEqualBound(x, bound1)
    abovebounds = greatEqualBound(x, bound2)
    mask        = np.logical_not(np.logical_or(
                      belowbounds, abovebounds ))
    mask[~finite] = missingValue
    return mask


#=========================================================
# ObsFunction classes to be accessed w/ ObsFunctionWrapper
# i.e., functions of variables contained in the database
#=========================================================

class GlintAngle:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.senzenMeta)
        self.baseVars.append(vu.senaziMeta)
        self.baseVars.append(vu.solzenMeta)
        self.baseVars.append(vu.solaziMeta)

    def evaluate(self, dbVals, insituParameters):
        senazi = dbVals[insituParameters[vu.senaziMeta]]
        solazi = dbVals[insituParameters[vu.solaziMeta]]

        relazi = np.abs(np.subtract(solazi,senazi))
        p = greatBound(relazi, 180.0, False)
        relazi[p] = np.subtract(360.0, relazi[p])
        relazi = np.multiply(np.subtract(180.0, relazi), vu.deg2rad)

        senzen = np.multiply(dbVals[insituParameters[vu.senzenMeta]], vu.deg2rad)
        solzen = np.multiply(dbVals[insituParameters[vu.solzenMeta]], vu.deg2rad)

        glint = np.add(np.multiply(np.cos(solzen), np.cos(senzen)),
                    np.multiply(np.sin(solzen),
                        np.multiply(np.sin(senzen), np.cos(relazi))))

        glint[greatBound(glint, 1.0)] = np.NaN
        glint[lessBound(glint, -1.0)] = np.NaN

        glint = np.multiply(np.arccos(glint), vu.rad2deg)
        glint[greatBound(glint, maxGlint, False)] = maxGlint

        return glint


class LocalHour:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.datetimeMeta)
        self.baseVars.append(vu.lonMeta)

    def evaluate(self, dbVals, insituParameters):
        TimeStr = dbVals[insituParameters[vu.datetimeMeta]]
        tzOffset = np.divide(dbVals[insituParameters[vu.lonMeta]], 15.0)

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


#classes related to the Asymmetric Cloud Impact (ACI)
class AsymmetricCloudImpact:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,insituParameters):
        # Minamide and Zhang (2018)
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        # BTbak = np.add(BTdep,BTobs)
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTclr = deepcopy(dbVals[insituParameters[vu.clrskyBTDiag]])
        p = lessBound(BTclr, 1.0, False)
        BTclr[p] = BTbak[p]
        # Minamide and Zhange (2018), maybe run into problems with non-clear-sky-bias-corrected
        ACI = np.subtract(np.abs(np.subtract(BTobs, BTclr)),
                          np.abs(np.subtract(BTbak, BTclr)))
# some ideas, but still causes bias in ACI for non-clear-sky-bias-corrected
#        zeros = np.zeros(BTbak.shape)
#        ACI =  np.subtract(np.maximum(zeros, np.subtract(BTclr, BTobs)),
#                           np.abs(np.subtract(BTbak, BTclr)))
##                           np.maximum(zeros, np.subtract(BTclr, BTbak)))

        return ACI


class ABEILambda:
    minLambda = 1.0
    maxLambda = 1.4
    def __init__(self):
        self.ACI = AsymmetricCloudImpact()
        self.baseVars = []
        self.baseVars = pu.uniqueMembers(self.baseVars + self.ACI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        osName = insituParameters['osName']
        if osName is None or osName not in ABEIParams:
            _logger.error('osName not available in ABEIParams => '+osName)
            os._exit(1)

        # varName, ch = vu.splitIntSuffix(insituParameters[vu.selfDepValue])
        varName, ch = vu.splitIntSuffix(insituParameters[vu.selfHofXValue])
        LambdaOverACI = ABEIParams[osName][(int(ch))]['LambdaOverACI']

        ACI = self.ACI.evaluate(dbVals, insituParameters)
        lambdaOut = np.ones(ACI.shape)
        crit = (ACI > 0.0)
        lambdaOut[crit] = np.multiply(ACI[crit], LambdaOverACI) + self.minLambda
        crit = (ACI >= (self.maxLambda - self.minLambda) / LambdaOverACI)
        lambdaOut[crit] = self.maxLambda

        return lambdaOut


#classes related to the Symmetric Cloud Impact (SCI)
class SCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self, dbVals, insituParameters):
        # Okamoto, et al.
        # Co = abs(Bias-Corrected BTobs - BTclr)
        # Cm = abs(BTbak - BTclr)

        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        # BTbak = np.add(BTdep, BTobs)
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTclr = deepcopy(dbVals[insituParameters[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        SCI = np.multiply( 0.5,
                 np.add(np.abs(np.subtract(BTobs, BTclr)),
                        np.abs(np.subtract(BTbak, BTclr))) )
        return SCI


class ScaledSCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self, dbVals, insituParameters):
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        # BTbak = np.add(BTdep, BTobs)
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTclr = deepcopy(dbVals[insituParameters[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[insituParameters[vu.cldfracMeta]]

        #Scale both Co and Cm by retrieved cloud fraction
        # SCI = np.multiply( 0.5,
        #          np.multiply(CldFrac,
        #              np.add(np.abs(np.subtract(BTobs, BTclr)),
        #                     np.abs(np.subtract(BTbak, BTclr))) ) )

        #Scale only Co by retrieved cloud fraction
        SCI = np.multiply( 0.5,
                  np.add(
                     np.multiply(CldFrac,
                            np.abs(np.subtract(BTobs, BTclr))),
                            np.abs(np.subtract(BTbak, BTclr)) ) )

        return SCI


class SCIModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self, dbVals, insituParameters):
        # Modified Harnisch, et al.
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        # BTbak = np.add(BTdep, BTobs)
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTclr = deepcopy(dbVals[insituParameters[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        zeros = np.full_like(BTbak,0.0)
        SCI = np.multiply( 0.5,
                 np.add(np.maximum(zeros, np.subtract(BTclr, BTobs)),
                        np.maximum(zeros, np.subtract(BTclr, BTbak))) )
        return SCI


class ScaledSCIModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self, dbVals, insituParameters):
        # Modified Harnisch, et al.
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        # BTbak = np.add(BTdep, BTobs)
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTclr = deepcopy(dbVals[insituParameters[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[insituParameters[vu.cldfracMeta]]

        zeros = np.full_like(BTbak,0.0)
        SCI = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.maximum(zeros, np.subtract(BTclr, BTobs)),
                            np.maximum(zeros, np.subtract(BTclr, BTbak))) ) )
        return SCI


class NormalizedError:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfHofXValue)
        self.baseVars.append(vu.selfErrorValue)

    def evaluate(self, dbVals, insituParameters):
        BTerr = dbVals[insituParameters[vu.selfErrorValue]]
        BTerr[BTerr==0.0] = np.NaN
        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTdep = np.subtract(BTbak, BTobs)

        return np.divide(BTdep, BTerr)


mpasFCRes = 120 # km [30, 120]
biasCorrectType = {} #[None, 'constant', 'varbc']
biasCorrectType['abi_g16'] = 'constant'
biasCorrectType['ahi_himawari8'] = None

class SCINormalizedError:
    def __init__(self):
        self.baseVars = []
        # self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfHofXValue)

    def evaluate(self, dbVals, insituParameters, SCISTDName, SCI):
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
        #        SCI0 SCI1
        #---------------------------------------------
        osName = insituParameters['osName']
        SCIErrParams = deepcopy(allSCIErrParams[(mpasFCRes,biasCorrectType.get(osName,None))])

        if osName is None or osName not in SCIErrParams:
            _logger.error('osName not available in SCIErrParams => '+osName)
            os._exit(1)

        # varName, ch = vu.splitIntSuffix(insituParameters[vu.selfDepValue])
        varName, ch = vu.splitIntSuffix(insituParameters[vu.selfHofXValue])
        STD0 = SCIErrParams[osName][(int(ch), SCISTDName)]['ERR'][0]
        STD1 = SCIErrParams[osName][(int(ch), SCISTDName)]['ERR'][1]
        SCI0  = SCIErrParams[osName][(int(ch), SCISTDName)]['X'][0]
        SCI1  = SCIErrParams[osName][(int(ch), SCISTDName)]['X'][1]
        slope = (STD1 - STD0) / (SCI1 - SCI0)

        belowramp = lessEqualBound(SCI, SCI0, False)
        aboveramp = greatEqualBound(SCI, SCI1, False)
        onramp    = betweenBounds(SCI, SCI0, SCI1, False)

        BTerr = np.full_like(SCI, np.NaN)
        BTerr[belowramp] = STD0
        BTerr[onramp]    = STD0 + slope * (SCI[onramp] - SCI0)
        BTerr[aboveramp] = STD1

        # BTdep = dbVals[insituParameters[vu.selfDepValue]]
        BTobs = dbVals[insituParameters[vu.selfObsValue]]
        BTbak = dbVals[insituParameters[vu.selfHofXValue]]
        BTdep = np.subtract(BTbak, BTobs)

        return np.divide(BTdep, BTerr)


class OkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.SCI = SCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        SCI = self.SCI.evaluate(dbVals, insituParameters)
        return super().evaluate(dbVals, insituParameters, OkamotoMethod, SCI)


class ScaledOkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.SCI = ScaledSCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        SCI = self.SCI.evaluate(dbVals, insituParameters)
        return super().evaluate(dbVals, insituParameters, ScaleOkamotoMethod, SCI)


class ModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.SCI = SCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        SCI = self.SCI.evaluate(dbVals, insituParameters)
        return super().evaluate(dbVals, insituParameters, ModHarnischMethod, SCI)


class ScaledModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.SCI = ScaledSCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, insituParameters):
        SCI = self.SCI.evaluate(dbVals, insituParameters)
        return super().evaluate(dbVals, insituParameters, ScaleModHarnischMethod, SCI)

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
# generic wrappers for ObsFunction classes
#=========================================
class BaseObsFunction:
    def __init__(self, baseVars):
        self.baseVars = deepcopy(baseVars)
        pass

    def dbVars(self, varName, fileFormat, outerIters_):
        dbVars = []

        if (not isinstance(outerIters_, Iterable)
           or isinstance(outerIters_,str)):
            outerIters = [outerIters_]
        else:
            outerIters = outerIters_

        for baseVar in self.baseVars:
            for outerIter in outerIters:
                dbVar = vu.base2dbVar(
                    baseVar, varName, fileFormat, outerIter)
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)


class IdObsFunction(BaseObsFunction):
    def __init__(self, variable):
        super().__init__([variable])

    def evaluate(self, dbVals, insituParameters):
        return dbVals[insituParameters[self.baseVars[0]]]


class ObsFunction(BaseObsFunction):
    def __init__(self, function):
        self.function = function()
        assert hasattr(self.function, 'baseVars'), \
            ("ERROR, function class must have the baseVars attribute:", function)
        super().__init__(self.function.baseVars)

    def evaluate(self, dbVals, insituParameters):
        return self.function.evaluate(dbVals, insituParameters)


class ObsFunctionWrapper:
    def __init__(self, config):
        self.osName = config['osName']
        self.fileFormat = config['fileFormat']

        variable = config['variable']
        varIsString = isinstance(variable,str)
        varIsClass = inspect.isclass(variable)
        assert varIsString ^ varIsClass, \
            ("ERROR: 'variable' must either be a String or a Class", config)

        if varIsString:
            self.function = IdObsFunction(variable)

        if varIsClass:
            self.function = ObsFunction(variable)

    def dbVars(self, varName, outerIters):
        return self.function.dbVars(varName, self.fileFormat, outerIters)

    def evaluate(self, dbVals, varName, outerIter):
        # setup context-specific insituParameters for the evaluation
        insituParameters = {}
        for baseVar in self.function.baseVars:
            insituParameters[baseVar] = vu.base2dbVar(
                    baseVar, varName, self.fileFormat, outerIter)
        insituParameters['osName'] = self.osName

        # evaluate the function
        self.result = self.function.evaluate(dbVals, insituParameters)


#========================
# generic binning classes
#========================

class BinFilter:
    def __init__(self, config):
        self.where  = config['where']
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
            self.bounds = np.empty(nBins, dtype=type(tmp))
            for ii in ibins:
                self.bounds[ii] = tmp
        else:
            self.bounds = np.empty(nBins, dtype=type(tmp[0]))

            # single element Iterable is applied uniformly to all bins
            if len(tmp) == 1:
                for ii in ibins:
                    self.bounds[ii] = tmp[0]
            # multiple element Iterable must be same length as nBins
            elif len(tmp) == nBins:
                for ii in ibins:
                    self.bounds[ii] = tmp[ii]
            else:
                _logger.error("'bounds' must be a scalar, single-member Iterable, or an Iterable with the same length as 'values'!")
                os._exit(1)
        self.function = ObsFunctionWrapper(config)
        self.except_diags = config.get('except_diags', [])
        self.mask_value = config.get('mask_value', np.NaN)
        #TODO: add other actions besides mask_value/exclude

#    def baseVars(self):
#        return pu.uniqueMembers(self.function.baseVars)

    def dbVars(self, varName, outerIters):
        dbVars = self.function.dbVars(
            varName, outerIters)
        return pu.uniqueMembers(dbVars)

    def evaluate(self, dbVals, varName, outerIter):
        self.function.evaluate(dbVals, varName, outerIter)

    def apply(self, array, diagName, ibin):
        newArray = deepcopy(array)

        if diagName not in self.except_diags:
            # remove locations where the mask is True
            mask = self.where(self.function.result,(self.bounds)[ibin])

            if len(mask) == len(newArray):
                newArray[mask] = self.mask_value
            else:
                _logger.error('BinFilter mask is incorrectly defined!')
                os._exit(1)

        return newArray


exclusiveDiags = ['obs','bak','ana','SCI']

class BinMethod:
    def __init__(self, config):
        #allows for scalar, str, and Iterable 'values'
        tmp = config['values']
        self.values = []
        if (not isinstance(tmp, Iterable) or
            isinstance(tmp, str)):
            self.values += [tmp]
        else:
            self.values += tmp

        self.excludeDiags = deepcopy(exclusiveDiags)
        override = config.get('override_exclusiveDiags',[])
        for diag in override:
            if diag in self.excludeDiags:
                self.excludeDiags.remove(diag)

        fconf = {}
        fconf['osName'] = config['osName']
        fconf['fileFormat'] = config['fileFormat']
        fconf['nBins'] = len(self.values)

        self.filters = []
        for filterConf in config['filters']:
            filterConf.update(fconf)
            self.filters.append(BinFilter(filterConf))

        enoughBounds = False
        for Filter in self.filters:
            if len(Filter.bounds) == len(self.values):
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

    def evaluate(self, dbVals, varName, outerIter):
        for ii in list(range(len(self.filters))):
            self.filters[ii].evaluate(
                dbVals, varName, outerIter)

    def apply(self, array, diagName, binVal):
        ibin = self.values.index(binVal)
        masked_array = deepcopy(array)
        for Filter in self.filters:
            masked_array = Filter.apply(
                masked_array, diagName, ibin)
        return masked_array
