#!/usr/bin/env python3

from collections.abc import Iterable
from copy import deepcopy
import inspect
import numpy as np
import os
import plot_utils as pu
from binning_params import allSCIErrParams
import var_utils as vu


#===========================
# binning method descriptors
#===========================
# TODO(JJG): should be moved to binning_params

# identity
identityBinMethod = 'identity'

#QC
goodQCMethod = 'good'
badQCMethod = 'bad'

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

# glint angle
maxGlint = 90.0

#========================
# binning where functions
#========================

# Note: NaN values have mask set to true by default

def equalBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.equal(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def notEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.not_equal(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def lessEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.less_equal(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def lessBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.less(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def greatEqualBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.greater_equal(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def greatBound(x, bound, nanmask=True):
    nanlocs = np.isnan(x)
    mask = np.empty_like(x, dtype=bool)
    mask[~nanlocs] = np.greater(x[~nanlocs], bound)
    mask[nanlocs] = nanmask
    return mask

def betweenBounds(x, bound1, bound2, nanmask=True):
    nanlocs = np.isnan(x)
    belowbounds = lessEqualBound(x, bound1)
    abovebounds = greatEqualBound(x, bound2)
    mask        = np.logical_not(np.logical_or(
                      belowbounds, abovebounds ))
    mask[nanlocs] = nanmask
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

    def evaluate(self, dbVals, caseParams):
        senazi = dbVals[vu.senaziMeta]
        solazi = dbVals[vu.solaziMeta]

        relazi = np.abs(np.subtract(solazi,senazi))
        relazi[relazi > 180.0] = np.subtract(360.0,relazi[relazi > 180.0])
        relazi = np.multiply(np.subtract(180.0,relazi), vu.deg2rad)

        senzen = np.multiply(dbVals[vu.senzenMeta], vu.deg2rad)
        solzen = np.multiply(dbVals[vu.solzenMeta], vu.deg2rad)

        glint = np.add(np.multiply(np.cos(solzen), np.cos(senzen)),
                    np.multiply(np.sin(solzen),
                        np.multiply(np.sin(senzen), np.cos(relazi))))

        glint[glint >  1.0] = np.NaN
        glint[glint < -1.0] = np.NaN

        glint = np.multiply(np.arccos(glint), vu.rad2deg)
        glint[glint > maxGlint] = maxGlint

        return glint


class LocalHour:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.dtMeta)
        self.baseVars.append(vu.lonMeta)

    def evaluate(self, dbVals, caseParams):
        TimeStr = dbVals[vu.dtMeta]
        tzOffset = np.divide(dbVals[vu.lonMeta],15.0)

        LH = np.full_like(tzOffset,0.0)
        t0 = LH0
        t1 = LH1
        dt = LHDT
        for ii, time in enumerate(TimeStr):
            ## Expecting time to fit YYYY-MM-DDThh:mm:ssZ
            # YYYY = float(time[0:4])
            # MMo  = float(time[5:7])
            # DD   = float(time[8:10])
            hh   = float(time[11:13])
            mmi  = float(time[14:16]) / 60.0
            ss   = float(time[17:19]) / 3600.0

            LH[ii] = hh + mmi + ss + tzOffset[ii]
            if LH[ii] < t0 - 0.5*dt: LH[ii] += 24.0
            LH[ii] = LH[ii] % 24.0
            if LH[ii] >= t1 + 0.5*dt: LH[ii] -= 24.0

        return LH


#classes related to the Asymmetric Cloud Impact (ACI)
class AsymmetricCloudImpact:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,caseParams):
        # Minamide and Zhang, 2018
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        ACE = np.subtract(np.abs(np.subtract(BTobs,BTclr)),
                          np.abs(np.subtract(BTbak,BTclr)))
        return ACE


#classes related to the Symmetric Cloud Impact (SCI)
class SCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self, dbVals, caseParams):
        # Okamoto, et al.
        # Co = abs(Bias-Corrected BTobs - BTclr)
        # Cm = abs(BTbak - BTclr)

        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep, BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        SCI = np.multiply( 0.5,
                 np.add(np.abs(np.subtract(BTobs, BTclr)),
                        np.abs(np.subtract(BTbak, BTclr))) )
        return SCI


class ScaledSCIOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self, dbVals, caseParams):
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep, BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[caseParams['base2db'][vu.cldfracMeta]]

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
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self, dbVals, caseParams):
        # Modified Harnisch, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep, BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
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
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self, dbVals, caseParams):
        # Modified Harnisch, et al.
        BTobs = dbVals[caseParams['base2db'][vu.selfObsValue]]
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]
        BTbak = np.add(BTdep, BTobs)
        BTclr = deepcopy(dbVals[caseParams['base2db'][vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[caseParams['base2db'][vu.cldfracMeta]]

        zeros = np.full_like(BTbak,0.0)
        SCI = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.maximum(zeros, np.subtract(BTclr, BTobs)),
                            np.maximum(zeros, np.subtract(BTclr, BTbak))) ) )
        return SCI


class NormalizedError:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.selfErrorValue)

    def evaluate(self, dbVals, caseParams):
        BTerr = dbVals[caseParams['base2db'][vu.selfErrorValue]]
        BTerr[BTerr==0.0] = np.NaN
        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]

        return np.divide(BTdep, BTerr)


mpasFCRes = 120 # km [30, 120]
biasCorrectType = {} #[None, 'constant', 'varbc']
biasCorrectType['abi_g16'] = 'constant'
biasCorrectType['ahi_himawari8'] = None

class SCINormalizedError:
    def __init__(self):
        pass

    def evaluate(self, dbVals, caseParams, SCISTDName, SCI):
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
        osName = caseParams['osName']
        SCIErrParams = deepcopy(allSCIErrParams[(mpasFCRes,biasCorrectType.get(osName,None))])

        if osName is None or osName not in SCIErrParams:
            print("ERROR: osName not available in SCIErrParams => "+osName)
            os._exit(1)

        varName, ch = vu.splitIntSuffix(caseParams['base2db'][vu.selfDepValue])
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

        BTdep = dbVals[caseParams['base2db'][vu.selfDepValue]]

        return np.divide(BTdep, BTerr)


class OkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = SCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, caseParams):
        SCI = self.SCI.evaluate(dbVals, caseParams)
        return super().evaluate(dbVals, caseParams, OkamotoMethod, SCI)


class ScaledOkamotoNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = ScaledSCIOkamoto()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, caseParams):
        SCI = self.SCI.evaluate(dbVals, caseParams)
        return super().evaluate(dbVals, caseParams, ScaleOkamotoMethod, SCI)


class ModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = SCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, caseParams):
        SCI = self.SCI.evaluate(dbVals, caseParams)
        return super().evaluate(dbVals, caseParams, ModHarnischMethod, SCI)


class ScaledModHarnischNormalizedError(SCINormalizedError):
    def __init__(self):
        super().__init__()
        self.baseVars = []
        self.baseVars.append(vu.selfDepValue)
        self.SCI = ScaledSCIModHarnisch()
        self.baseVars = pu.uniqueMembers(self.baseVars + self.SCI.baseVars)

    def evaluate(self, dbVals, caseParams):
        SCI = self.SCI.evaluate(dbVals, caseParams)
        return super().evaluate(dbVals, caseParams, ScaleModHarnischMethod, SCI)

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

    def dbVars(self, varName, outerIters_):
        dbVars = []

        if not isinstance(outerIters_, Iterable):
            outerIters = [outerIters_]
        else:
            outerIters = outerIters_

        for base in self.baseVars:
            for outerIter in outerIters:
                dbVar = vu.base2dbVar(
                    base, varName, outerIter)
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)


class IdObsFunction(BaseObsFunction):
    def __init__(self, variable):
        super().__init__([variable])

    def evaluate(self, dbVals, caseParams):
        return dbVals[caseParams['base2db'][self.baseVars[0]]]


class ObsFunction(BaseObsFunction):
    def __init__(self, function):
        self.function = function()
        assert hasattr(self.function, 'baseVars'), \
            ("ERROR, function class must have the baseVars attribute:", function)
        super().__init__(self.function.baseVars)

    def evaluate(self, dbVals, caseParams):
        return self.function.evaluate(dbVals, caseParams)


class ObsFunctionWrapper:
    def __init__(self, config):
        self.osName = config.get('osName', None)
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
        return self.function.dbVars(varName, outerIters)

    def evaluate(self, dbVals, varName, outerIter):
        caseParams = {}
        caseParams['base2db'] = {}
        for base in self.function.baseVars:
            caseParams['base2db'][base] = vu.base2dbVar(
                    base, varName, outerIter)
        caseParams['osName'] = self.osName
        self.result = self.function.evaluate(dbVals, caseParams)


#========================
# generic binning classes
#========================

class BinFilter:
    def __init__(self, config):
        self.where  = config['where']
        tmp         = config['bounds']

        # allow for scalar and iterable bounds
        self.bounds = np.empty(config['nBins'])
        if (not isinstance(tmp, Iterable) or
            len(tmp) == 1):
            self.bounds[:] = tmp
        elif len(tmp) == config['nBins']:
            for ii in list(range(config['nBins'])):
                self.bounds[ii] = tmp[ii]
        else:
            print("ERROR: 'bounds' need to be a scalar or an Iterable with the same length as 'values'!")
            os._exit(1)

        self.function = ObsFunctionWrapper(config)
        self.except_diags = config.get('except_diags', [])
        self.mask_value = config.get('mask_value', np.NaN)
        #TODO: add other actions besides mask_value/blacklist

#    def baseVars(self):
#        return pu.uniqueMembers(self.function.baseVars)

    def dbVars(self, varName, outerIters):
        dbVars = self.function.dbVars(
            varName, outerIters)
        return pu.uniqueMembers(dbVars)

    def evaluate(self, dbVals, varName, outerIter):
        self.function.evaluate(dbVals, varName, outerIter)

    def apply(self, array, diagName, ibin):
        # blacklist locations where the mask is True
        mask = self.where(self.function.result,(self.bounds)[ibin])

        if diagName not in self.except_diags:
            if len(mask) == len(array):
                array[mask] = self.mask_value
            else:
                print('\n\nERROR: BinFilter mask is incorrectly defined!')
                os._exit(1)

        return array


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

    def dbVars(self, varName, outerIters=['0']):
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
