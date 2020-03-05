from collections.abc import Iterable
from copy import deepcopy
import numpy as np
import os
import plot_utils as pu
import var_utils as vu

#==========================
# names and values of bins
#==========================

# heterogeneous/named bins
# latitude bands, north to south
allNamedLatBands = {}
allNamedLatBands['values']    = ['NPol','NXTro','Tro','SXTro','SPol']
allNamedLatBands['minBounds'] = [60.0, 30.0, -30.0, -90.0, -90.0]
allNamedLatBands['maxBounds'] = [90.0, 90.0,  30.0, -30.0, -60.0]

namedLatBands = {}
namedLatBands['values'] = ['NXTro','Tro','SXTro']
namedLatBands['minBounds'] = []
namedLatBands['maxBounds'] = []

for latBand in namedLatBands['values']:
    iband = allNamedLatBands['values'].index(latBand)
    namedLatBands['minBounds'].append(allNamedLatBands['minBounds'][iband])
    namedLatBands['maxBounds'].append(allNamedLatBands['maxBounds'][iband])

# homogeneous bins
binLims = {}

binLims[vu.obsVarPrs] = {}
binLims[vu.obsVarPrs]['start']  = 1100.0
binLims[vu.obsVarPrs]['finish'] = 0.0
binLims[vu.obsVarPrs]['step']   = -100.0
binLims[vu.obsVarPrs]['format'] = '{:.0f}'

binLims[vu.obsVarAlt] = {}
binLims[vu.obsVarAlt]['start']  = 1000.0
binLims[vu.obsVarAlt]['finish'] = 50000.0
binLims[vu.obsVarAlt]['step']   = 2000.0
binLims[vu.obsVarAlt]['format'] = '{:.0f}'

binLims[vu.obsVarLat] = {}
binLims[vu.obsVarLat]['start']  = 90.0
binLims[vu.obsVarLat]['finish'] = -90.0
binLims[vu.obsVarLat]['step']   = -10.0
binLims[vu.obsVarLat]['format'] = '{:.0f}'

binLims[vu.obsVarSatZen] = {}
binLims[vu.obsVarSatZen]['start']  = 0.0
binLims[vu.obsVarSatZen]['finish'] = 70.0
binLims[vu.obsVarSatZen]['step']   = 5.0
binLims[vu.obsVarSatZen]['format'] = '{:.0f}'

binLims[vu.obsVarCldFrac] = {}
binLims[vu.obsVarCldFrac]['start']  = 0.0
binLims[vu.obsVarCldFrac]['finish'] = 1.0
binLims[vu.obsVarCldFrac]['step']   = 0.2
binLims[vu.obsVarCldFrac]['format'] = '{:.2f}'

binLims[vu.obsVarCa] = {}
binLims[vu.obsVarCa]['start']  = 0.0
binLims[vu.obsVarCa]['finish'] = 60.0
binLims[vu.obsVarCa]['step']   = 1.0
binLims[vu.obsVarCa]['format'] = '{:.0f}'

#binLims[vu.modVarLat] = {}
#binLims[vu.modVarLat]['start']  = 90.0
#binLims[vu.modVarLat]['finish'] = -90.0
#binLims[vu.modVarLat]['step']   = -30.0
#binLims[vu.modVarLat]['format'] = '{:.0f}'

#binLims[vu.modVarAlt] = {}
#binLims[vu.modVarAlt]['start']  = 0.0
#binLims[vu.modVarAlt]['finish'] = 50000.0
#binLims[vu.modVarAlt]['step']   = 5000.0
#binLims[vu.modVarAlt]['format'] = '{:.0f}'

for binType, param in binLims.items():
    binBounds = list(np.arange(
        param['start']-0.5*np.abs(param['step']),
        param['finish']+0.5*param['step'],
        param['step']))
    binLims[binType]['minBounds'] = []
    binLims[binType]['maxBounds'] = []
    binLims[binType]['values'] = []
    for ibin in list(range(len(binBounds)-1)):
        if param['step'] > 0.0:
            binLims[binType]['minBounds'].append(binBounds[ibin])
            binLims[binType]['maxBounds'].append(binBounds[ibin+1])
        else:
            binLims[binType]['minBounds'].append(binBounds[ibin+1])
            binLims[binType]['maxBounds'].append(binBounds[ibin])

        binVal = 0.5 * (binBounds[ibin+1] + binBounds[ibin])
        binLims[binType]['values'].append(param['format'].format(binVal))


#jet-stream bin method specifications
P_jet_min = 250.0
P_jet_max = 350.0
P_jet_val = '{:.0f}'.format(0.5 * (P_jet_min + P_jet_max))
PjetMethod = 'P='+P_jet_val+'hPa'

alt_jet_min = 9500.0
alt_jet_max = 10500.0
alt_jet_val = '{:.0f}'.format(0.5 * (alt_jet_min + alt_jet_max))
altjetMethod = 'alt='+alt_jet_val+'m'


# clear/cloudy bin method specifications
clrskyMethod = 'clear-sky'
clrskyThresh = 0.2

cldskyMethod = 'cloud-sky'
cldskyThresh = 0.8


#@EffectiveQC* values:
# pass    = 0;   // we like that one!
# missing = 1;   // missing values prevent use of observation
# preQC   = 2;   // observation rejected by pre-processing
# bounds  = 3;   // observation value out of bounds
# domain  = 4;   // observation not within domain of use
# black   = 5;   // observation black listed
# Hfailed = 6;   // H(x) computation failed
# thinned = 7;   // observation removed due to thinning
# diffref = 8;   // metadata too far from reference
# clw     = 9;   // observation removed due to cloud field
# fguess  = 10;  // observation too far from guess
# seaice  = 11;  // observation based sea ice detection, also flags land points
#
# Static list above copied on 3 Oct 2019
# see ufo/src/ufo/filters/QCflags.h for up-to-date list

goodFlag = 0
goodFlagName = 'pass'

badFlags     = [1,         2,         3,         4,
                5,         6,         7,         8,
                9,         10,        11]
badFlagNames = ['missing', 'preQC',   'bounds',  'domain',
                'black',   'Hfailed', 'thinned', 'diffref',
                'clw',     'fguess',  'seaice']


#========================
# binning where functions
#========================

# Note: NaN values have mask set to true for blacklisting

def equalBound(x, bound):
    nonnan = ~np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[nonnan] = np.equal(x[nonnan],bound)
    mask[~nonnan] = True
    return mask

def notEqualBound(x, bound):
    nonnan = ~np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[nonnan] = np.not_equal(x[nonnan],bound)
    mask[~nonnan] = True
    return mask

#def lessEqualBound(x, bound):
#    nonnan = ~np.isnan(x)
#    mask = np.empty_like(x,dtype=bool)
#    mask[nonnan] = np.less_equal(x[nonnan],bound)
#    mask[~nonnan] = True
#    return mask

def lessBound(x, bound):
    nonnan = ~np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[nonnan] = np.less(x[nonnan],bound)
    mask[~nonnan] = True
    return mask

def greatEqualBound(x, bound):
    nonnan = ~np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[nonnan] = np.greater_equal(x[nonnan],bound)
    mask[~nonnan] = True
    return mask

def greatBound(x, bound):
    nonnan = ~np.isnan(x)
    mask = np.empty_like(x,dtype=bool)
    mask[nonnan] = np.greater(x[nonnan],bound)
    mask[~nonnan] = True
    return mask


#========================================================
# specific function classes
# i.e., functions of variables contained in the database
#========================================================

class CaOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,base2db):
        # Okamoto, et al.
        BTobs = dbVals[base2db[vu.selfObsValue]]
        BTdep = dbVals[base2db[vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[base2db[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        Ca = np.multiply( 0.5,
                 np.add(np.abs(np.subtract(BTobs,BTclr)),
                        np.abs(np.subtract(BTbak,BTclr))) )
        return Ca


class ScaledCaOkamoto:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self,dbVals,base2db):
        # Okamoto, et al.
        BTobs = dbVals[base2db[vu.selfObsValue]]
        BTdep = dbVals[base2db[vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[base2db[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[base2db[vu.cldfracMeta]]

        Ca = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.abs(np.subtract(BTobs,BTclr)),
                            np.abs(np.subtract(BTbak,BTclr))) ) )
        return Ca


class CaModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)

    def evaluate(self,dbVals,base2db):
        # Modified Harnisch, et al.
        BTobs = dbVals[base2db[vu.selfObsValue]]
        BTdep = dbVals[base2db[vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[base2db[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        zeros = np.full_like(BTbak,0.0)
        Ca = np.multiply( 0.5,
                 np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                        np.maximum(zeros,np.subtract(BTclr,BTbak))) )
        return Ca


class ScaledCaModHarnisch:
    def __init__(self):
        self.baseVars = []
        self.baseVars.append(vu.selfObsValue)
        self.baseVars.append(vu.selfDepValue)
        self.baseVars.append(vu.clrskyBTDiag)
        self.baseVars.append(vu.cldfracMeta)

    def evaluate(self,dbVals,base2db):
        # Modified Harnisch, et al.
        BTobs = dbVals[base2db[vu.selfObsValue]]
        BTdep = dbVals[base2db[vu.selfDepValue]]
        BTbak = np.add(BTdep,BTobs)
        BTclr = deepcopy(dbVals[base2db[vu.clrskyBTDiag]])
        BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
        CldFrac = dbVals[base2db[vu.cldfracMeta]]

        zeros = np.full_like(BTbak,0.0)
        Ca = np.multiply( 0.5,
                 np.multiply(CldFrac,
                     np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                            np.maximum(zeros,np.subtract(BTclr,BTbak))) ) )
        return Ca

#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(dbVals,REGION_NAME):
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = dbVals['longitdue']
#    lats = dbVals['latitude']
#    nlocs = len(lons)
#    REGIONS = isinsideregions(lons,lats,shp)
#    return REGIONS


#======================
# generic bin classes
#======================
class BaseFilterFunc:
    def __init__(self, baseVars):
        self.baseVars = deepcopy(baseVars)
        pass

    def dbVars(self,varName,outerIters):
        dbVars = []
        for base in self.baseVars:
            for outerIter in outerIters:
                dbVar = vu.base2dbVar(
                    base,varName,outerIter)
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)


class IdFilterFunc(BaseFilterFunc):
    def __init__(self,variable):
        super().__init__([variable])
        self.result = []

    def evaluate(self,dbVals,base2db):
        self.result = dbVals[base2db[self.baseVars[0]]]


class FilterFuncWrapper(BaseFilterFunc):
    def __init__(self,function):
        self.function = function()
        super().__init__(self.function.baseVars)
        self.result = []

    def evaluate(self,dbVals,base2db):
        self.result = self.function.evaluate(dbVals,base2db)


class BinFilter:
    def __init__(self,config,nBins):
        self.where          = config['where']
        tmp                 = config['bounds']

        # allow for scalar and iterable inputs
        if not isinstance(tmp, Iterable):
            self.bounds = [tmp]*nBins
        elif len(tmp) == 1:
            self.bounds = tmp*nBins
        elif len(tmp) == nBins:
            self.bounds = tmp
        else:
            print("ERROR: 'bounds' needs to be a scalar or same length as 'values'")
            os._exit(1)

        self.apply_to_diags = config.get('apply_to_diags',vu.allDiags)
        self.mask_value     = config.get('mask_value',np.NaN)
        #TODO: add other actions besides mask_value/blacklist

        function = config.get('function',None)
        variable = config.get('variable',None)
        if variable is not None and function is None:
            # print("VARIABLE ",variable)
            self.function = IdFilterFunc(variable)
        elif function is not None and variable is None:
            # print("FUNC ",function)
            self.function = FilterFuncWrapper(function)
        else:
            print("ERROR: either 'variable' or 'function' must be provided to BinFilter constructor, but never both")
            os._exit(1)

#    def baseVars(self):
#        return pu.uniqueMembers(self.function.baseVars)

    def dbVars(self,varName,outerIters):
        dbVars = []
        for dbVar in self.function.dbVars(
            varName,outerIters):
            dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)

    def evaluate(self,dbVals,varName,outerIter):
        base2db = {}
        for base in self.function.baseVars:
            base2db[base] = vu.base2dbVar(
                    base,varName,outerIter)

        self.function.evaluate(dbVals,base2db)

    def apply(self,array,diagName,ibin):
        # blacklist locations where the mask is True
        mask = self.where(self.function.result,(self.bounds)[ibin])

        if diagName in self.apply_to_diags:
            if len(mask) == len(array):
                array[mask] = self.mask_value
            else:
                print('\n\nERROR: BinFilter : mask is incorrectly defined!')
                os._exit(1)

        return array


class BinMethod:
    def __init__(self, config):
        #handle scalar/str and iterable values inputs
        tmp = config['values']
        if (not isinstance(tmp, Iterable) or
            isinstance(tmp,str)):
            self.values = [tmp]
        else:
            self.values = tmp

        self.filters = []
        for filterConf in config['filters']:
            self.filters.append(
                BinFilter(filterConf,len(self.values)) )

        enoughBounds = False
        for Filter in self.filters:
            if len(Filter.bounds) == len(self.values):
                enoughBounds = True
        if not enoughBounds:
            print('\n\nERROR: BinMethod : at least one filter must have len(bounds) == len(values)!')
            os._exit(1)

#    def baseVars(self):
#        baseVars = []
#        for Filter in self.filters:
#            for variable in Filter.baseVars():
#                baseVars.append(variable)
#        return pu.uniqueMembers(baseVars)

    def dbVars(self,varName,outerIters=['0']):
        dbVars = []
        for Filter in self.filters:
            for dbVar in Filter.dbVars(
                varName,outerIters):
                dbVars.append(dbVar)
        return pu.uniqueMembers(dbVars)

    def evaluate(self,dbVals,varName,outerIter):
        for ii in list(range(len(self.filters))):
            self.filters[ii].evaluate(
                dbVals,varName,outerIter)

    def apply(self,array,diagName,binVal):
        ibin = self.values.index(binVal)
        masked_array = deepcopy(array)
        for Filter in self.filters:
            masked_array = Filter.apply(
                masked_array,diagName,ibin)
        return masked_array


#=================================================
# binning configurations used for all bin methods
#=================================================

# Each binVarConfig member has the following properties
# key: binVar string describes the variable to bin over
# variables: list of NC variables used by the binning filters
# binMethod (e.g., defaultBinMethod, NAMED, etc.): used to distinguish
#  between multiple methods with the same binVar
#     filters: list of filters that will be applied for each method
#     for each filter:
#         where: logical function that determines locations that are blacklisted
#         variable (optional): variable that is used to initialize the IdFilterFunc class
#         function (optional): the where test is applied to function.values. defaults to None
#         bounds: a list of numerical values used in the where test (scalar or Iterable same length as values). At least one filter must have len(bounds)==len(values).
#         apply_to_diags (optional): list of diagnostics for which a particular filter applies
#     values (string): list of values associated with each bin
#
# Note: either variable or function must be provided to each filter
#
#TODO: add "diags" selection for each binMethod (some only need to be calculated for omb/oma)
#TODO: classify each binmethod as 1D, 2D, etc. to indicate which types of figures apply to it

# standard binMethod name
defaultBinMethod = 'default'

nullBinMethod = { 'filters': [], 'values': [] }
nullBinVarConfig = { defaultBinMethod:nullBinMethod }

binVarConfigs = {
    vu.obsVarQC:{
        defaultBinMethod:{
            'filters':[
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': goodFlagName,
        },
        'bad':{
            'filters':[
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': badFlags,
             'apply_to_diags': vu.nonObsDiags},
            {'where': equalBound,
             'variable': vu.selfQCValue,
             'bounds': badFlags,
             'apply_to_diags': vu.nonObsDiags,
             'mask_value': 0.0},
            ],
            'values': badFlagNames,
        },
    },
    vu.obsVarPrs:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': binLims[vu.obsVarPrs]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': binLims[vu.obsVarPrs]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarPrs]['values'],
        },
        PjetMethod:{
            'filters':[
# eliminate locations outside P_jet_min to P_jet_max
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_min},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': P_jet_val,
        },
    },
    vu.obsVarAlt:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': binLims[vu.obsVarAlt]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': binLims[vu.obsVarAlt]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarAlt]['values'],
        },
        altjetMethod:{
            'filters':[
# eliminate locations outside alt_jet_min to alt_jet_max
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_min},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': alt_jet_val,
        },
    },
    vu.obsVarLat:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        'NAMED':{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': namedLatBands['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': namedLatBands['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': namedLatBands['values'],
        },
        PjetMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_min},
            {'where': greatEqualBound,
             'variable': vu.prsMeta,
             'bounds': P_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        altjetMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_min},
            {'where': greatEqualBound,
             'variable': vu.altMeta,
             'bounds': alt_jet_max},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.latMeta,
             'bounds': binLims[vu.obsVarLat]['maxBounds']},
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
    },
    vu.obsVarSatZen:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
        clrskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': clrskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
        cldskyMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.satzenMeta,
             'bounds': binLims[vu.obsVarSatZen]['maxBounds']},
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': cldskyThresh},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarSatZen]['values'],
        },
    },
    vu.obsVarCldFrac:{
        defaultBinMethod:{
            'filters':[
            {'where': lessBound,
             'variable': vu.cldfracMeta,
             'bounds': binLims[vu.obsVarCldFrac]['minBounds']},
            {'where': greatEqualBound,
             'variable': vu.cldfracMeta,
             'bounds': binLims[vu.obsVarCldFrac]['maxBounds']},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': binLims[vu.obsVarCldFrac]['values'],
        },
    },
    vu.obsVarCa:{
        'Okamoto':{
            'filters':[
            {'where': lessBound,
             'function': CaOkamoto,
             'bounds': binLims[vu.obsVarCa]['minBounds']},
            {'where': greatEqualBound,
             'function': CaOkamoto,
             'bounds': binLims[vu.obsVarCa]['maxBounds']},
            ],
            'values': binLims[vu.obsVarCa]['values'],
        },
        'ScaledOkamoto':{
            'filters':[
            {'where': lessBound,
             'function': ScaledCaOkamoto,
             'bounds': binLims[vu.obsVarCa]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledCaOkamoto,
             'bounds': binLims[vu.obsVarCa]['maxBounds']},
            ],
            'values': binLims[vu.obsVarCa]['values'],
        },
        'ModHarnisch':{
            'filters':[
            {'where': lessBound,
             'function': CaModHarnisch,
             'bounds': binLims[vu.obsVarCa]['minBounds']},
            {'where': greatEqualBound,
             'function': CaModHarnisch,
             'bounds': binLims[vu.obsVarCa]['maxBounds']},
            ],
            'values': binLims[vu.obsVarCa]['values'],
        },
        'ScaledModHarnisch':{
            'filters':[
            {'where': lessBound,
             'function': ScaledCaModHarnisch,
             'bounds': binLims[vu.obsVarCa]['minBounds']},
            {'where': greatEqualBound,
             'function': ScaledCaModHarnisch,
             'bounds': binLims[vu.obsVarCa]['maxBounds']},
            ],
            'values': binLims[vu.obsVarCa]['values'],
        },
    },
    'ObsRegion':{
        'AFRICA':    nullBinMethod,
        'ATLANTIC':  nullBinMethod,
        'AUSTRALIA': nullBinMethod,
        'CONUS':{
            'filters':[
            {'where': lessBound,
             'variable': vu.lonMeta,
             'bounds': [234.0]},
            {'where': greatBound,
             'variable': vu.lonMeta,
             'bounds': [294.0]},
            {'where': lessBound,
             'variable': vu.latMeta,
             'bounds': [ 25.0]},
            {'where': greatBound,
             'variable': vu.latMeta,
             'bounds': [ 50.0]},
            {'where': notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'apply_to_diags': vu.nonObsDiags},
            ],
            'values': ['CONUS'],
        },
        'EUROPE':    nullBinMethod,
        'E_EUROPE':  nullBinMethod,
        'NAMERICA':  nullBinMethod,
        'PACIFIC':   nullBinMethod,
        'SAMERICA':  nullBinMethod,
        'SE_ASIA':   nullBinMethod,
        'S_ASIA':    nullBinMethod,
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#        'CONUS_POLYGON':{
#            'filters':[
#            {'where': notInBound,
#             'function': Regions,
#             'variable': [vu.lonMeta,vu.latMeta],
#             'bounds': ['CONUS']},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE_POLYGON':{
#            'filters':[
#            {'where': notInBound,
#             'function': Regions,
#             'variable': [vu.lonMeta,vu.latMeta],
#             'bounds': ['EUROPE']},
#            ],
#            'values': ['EUROPE'],
#        },
    },
#TODO: enable binning in MPAS Model space
#    vu.modVarPrs:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarPrs,
#             'bounds': binLims[vu.modVarPrs]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarPrs,
#             'bounds': binLims[vu.modVarPrs]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarPrs]['values'],
#        },
#    },
#    vu.modVarAlt:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarAlt,
#             'bounds': binLims[vu.modVarAlt]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarAlt,
#             'bounds': binLims[vu.modVarAlt]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarAlt]['values'],
#        },
#    },
#    vu.modVarLat:{
#        defaultBinMethod:{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': binLims[vu.modVarLat]['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarLat,
#             'bounds': binLims[vu.modVarLat]['maxBounds']},
#            ],
#            'values: binLims[vu.modVarLat]['values'],
#        },
#        'NAMED':{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': namedLatBands['minBounds']},
#            {'where': greatEqualBound,
#             'variable': vu.modVarLat,
#             'bounds': namedLatBands['maxBounds']},
#            ],
#            'values': namedLatBands['values'],
#        },
#    },
#    'ModelRegion':{
#        'AFRICA':    nullBinMethod,
#        'ATLANTIC':  nullBinMethod,
#        'AUSTRALIA': nullBinMethod,
#        'CONUS':{
#            'filters':[
#            {'where': lessBound,
#             'variable': vu.modVarLon,
#             'bounds': [234.0]},
#            {'where': greatBound,
#             'variable': vu.modVarLon,
#             'bounds': [294.0]},
#            {'where': lessBound,
#             'variable': vu.modVarLat,
#             'bounds': [ 25.0]},
#            {'where': greatBound,
#             'variable': vu.modVarLat,
#             'bounds': [ 50.0]},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE':    nullBinMethod,
#        'E_EUROPE':  nullBinMethod,
#        'NAMERICA':  nullBinMethod,
#        'PACIFIC':   nullBinMethod,
#        'SAMERICA':  nullBinMethod,
#        'SE_ASIA':   nullBinMethod,
#        'S_ASIA':    nullBinMethod,
#    },
}


#=============================
# Parameterized binVarConfigs
#=============================

# Add named latitude-band-specific pressure bins
for iband, latBand in enumerate(namedLatBands['values']):
    binVarConfigs[vu.obsVarPrs][latBand] = {
        'filters':[
        {'where': lessBound,
         'variable': vu.prsMeta,
         'bounds': binLims[vu.obsVarPrs]['minBounds']},
        {'where': greatEqualBound,
         'variable': vu.prsMeta,
         'bounds': binLims[vu.obsVarPrs]['maxBounds']},
        {'where': lessBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['minBounds'][iband]},
        {'where': greatEqualBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['maxBounds'][iband]},
        {'where': notEqualBound,
         'variable': vu.selfQCValue,
         'bounds': goodFlag,
         'apply_to_diags': vu.nonObsDiags},
        ],
        'values': binLims[vu.obsVarPrs]['values'],
    }


# Add named latitude-band-specific altitude bins
for iband, latBand in enumerate(namedLatBands['values']):
    binVarConfigs[vu.obsVarAlt][latBand] = {
        'filters':[
        {'where': lessBound,
         'variable': vu.altMeta,
         'bounds': binLims[vu.obsVarAlt]['minBounds']},
        {'where': greatEqualBound,
         'variable': vu.altMeta,
         'bounds': binLims[vu.obsVarAlt]['maxBounds']},
        {'where': lessBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['minBounds'][iband]},
        {'where': greatEqualBound,
         'variable': vu.latMeta,
         'bounds': namedLatBands['maxBounds'][iband]},
        {'where': notEqualBound,
         'variable': vu.selfQCValue,
         'bounds': goodFlag,
         'apply_to_diags': vu.nonObsDiags},
        ],
        'values': binLims[vu.obsVarAlt]['values'],
    }

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
