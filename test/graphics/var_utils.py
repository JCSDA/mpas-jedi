from copy import deepcopy
import numpy as np
import os

miss_f = -88888.8
miss_i = -88888
miss_s = 'null'
csvSEP = ';'

#=====================
# variable definitions
#=====================

## NC variable names for MPAS-JEDI
obsVarAlt     = 'altitude'
obsVarCa      = 'symmetric_cloud_impact'
obsVarCldFrac = 'cloud_area_fraction'
obsVarLon     = 'longitude'
obsVarLat     = 'latitude'
obsVarPrs     = 'air_pressure'
obsVarQC      = 'QCflag'
obsVarSatZen  = 'sensor_zenith_angle'

degree= u'\N{DEGREE SIGN}'

# columns: var_name            unit_used   abbr.
varDictObs = {
    'air_temperature':        [ 'K',     'T'       ]
  , 'bending_angle':          [ '%',     'Bnd'     ]
  , 'brightness_temperature': [ 'K',     'BT'      ]
  , 'eastward_wind':          [ 'm/s',   'U'       ]
  , 'northward_wind':         [ 'm/s',   'V'       ]
  , 'refractivity':           [ '%',     'Ref'     ]
  , 'specific_humidity':      [ 'kg/kg', 'qv'      ]
  , 'surface_pressure':       [ 'Pa',    'Ps'      ]
  , 'virtual_temperature':    [ 'K',     'Tv'      ]
  , obsVarAlt:                [ 'm',     'alt'     ]
  , obsVarCa:                 [ 'K',     'Ca'      ]
  , obsVarCldFrac:            [ miss_s,  'cldfrac' ]
  , obsVarLat:                [ degree,  'lat'     ]
  , obsVarLon:                [ degree,  'lon'     ]
  , obsVarPrs:                [ 'hPa',   'P'       ]
  , obsVarQC:                 [ miss_s,  obsVarQC  ]
  , obsVarSatZen:             [ degree,  'zenith'  ]
    }
#Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == obsVarAlt:

# NC observation variable name substitutions
vNameStr = 'varName'
vChanStr = 'varCHAN'

# MPAS-JEDI suffixes for observation-type variables
bakGroup    = 'HofX'
depbgGroup  = 'depbg'
depanGroup  = 'depan'
diagGroup   = 'ObsDiag'
geovalGroup = 'GeoVaLs'
obsGroup    = 'ObsValue'
metaGroup   = 'MetaData'
qcGroup     = 'EffectiveQC'

# Fixed NC variable names
selfObsValue = vNameStr+'@'+obsGroup
selfBakValue = vNameStr+'@'+bakGroup
selfQCValue  = vNameStr+'@'+qcGroup

altMeta      = obsVarAlt+'@'+metaGroup
cldfracMeta  = obsVarCldFrac+'@'+metaGroup
lonMeta      = obsVarLon+'@'+metaGroup
latMeta      = obsVarLat+'@'+metaGroup
prsMeta      = obsVarPrs+'@'+metaGroup
satzenMeta   = obsVarSatZen+'@'+metaGroup

clrskyBTDiag = 'brightness_temperature_assuming_clear_sky_'+vChanStr+'@'+diagGroup

# functions for extracting sub-parts of UFO variable names
def splitObsVarGrp(varATgroup):
    if "@" in varATgroup:
        var = ''.join(varATgroup.split("@")[:-1])
        grp = ''.join(varATgroup.split("@")[-1])
    else:
        var = varATgroup
        grp = miss_s
    return var, grp

def splitChan(var):
    varName, grpName = splitObsVarGrp(var)
    ch = ''.join(varName.split("_")[-1:])
    return ch

## NC variable names for MPAS-Model
#modVarAlt = 'zgrid' # --> needs to be interpolated to nVertLevels instead of nVertLevelsP1
modVarPrs = 'pressure_p'
modVarLat = 'latCell'
modVarLon = 'lonCell'

# columns: var_name            unit_used   abbr.
varDictModel = {
    'theta':                  [ 'K',    'Theta']
  , 'pressure':               [ 'Pa',    'P'   ]
  , modVarPrs:                [ 'Pa',    'PP'  ]
  , 'rho':                    [ 'kg/m\N{SUPERSCRIPT THREE}', 'rho' ]
  , 'qv':                     [ 'g/kg',  'Qv'  ]
  , 'uReconstructZonal':      [ 'm/s',   'U'   ]
  , 'uReconstructMeridional': [ 'm/s',   'V'   ]
  , 'u':                      [ 'm/s',   'uedge']
  , modVarLat:                [ degree,  'lat' ]
  , modVarLon:                [ degree,  'lon' ]
    }
#Note, qv unit is kg/kg in original mpas restart file. The unit is converted to g/kg when read qv.


#===========================
# obs diagnostic definitions
#===========================

allDiags = ['omb','oma','obs','bak','ana']
nonObsDiags = [diag for diag in allDiags if diag!='obs']


#==========================
# names and values of bins
#==========================

# heterogeneous/named bins
# latitude bands, north to south
allNamedLatBands = {}
allNamedLatBands['labels']    = ['NPol','NXTro','Tro','SXTro','SPol']
allNamedLatBands['minBounds'] = [60.0, 30.0, -30.0, -90.0, -90.0]
allNamedLatBands['maxBounds'] = [90.0, 90.0,  30.0, -30.0, -60.0]

namedLatBands = {}
namedLatBands['labels'] = ['NXTro','Tro','SXTro']
namedLatBands['nbins'] = len(namedLatBands['labels'])
namedLatBands['minBounds'] = []
namedLatBands['maxBounds'] = []

for latBand in namedLatBands['labels']:
    iband = allNamedLatBands['labels'].index(latBand)
    namedLatBands['minBounds'].append(allNamedLatBands['minBounds'][iband])
    namedLatBands['maxBounds'].append(allNamedLatBands['maxBounds'][iband])

# homogeneous bins
binLims = {}

binLims[obsVarPrs] = {}
binLims[obsVarPrs]['start']  = 1100.0
binLims[obsVarPrs]['finish'] = 0.0
binLims[obsVarPrs]['step']   = -100.0
binLims[obsVarPrs]['format'] = '{:.0f}'

binLims[obsVarAlt] = {}
binLims[obsVarAlt]['start']  = 1000.0
binLims[obsVarAlt]['finish'] = 50000.0
binLims[obsVarAlt]['step']   = 2000.0
binLims[obsVarAlt]['format'] = '{:.0f}'

binLims[obsVarLat] = {}
binLims[obsVarLat]['start']  = 90.0
binLims[obsVarLat]['finish'] = -90.0
binLims[obsVarLat]['step']   = -10.0
binLims[obsVarLat]['format'] = '{:.0f}'

binLims[obsVarSatZen] = {}
binLims[obsVarSatZen]['start']  = 0.0
binLims[obsVarSatZen]['finish'] = 70.0
binLims[obsVarSatZen]['step']   = 5.0
binLims[obsVarSatZen]['format'] = '{:.0f}'

binLims[obsVarCldFrac] = {}
binLims[obsVarCldFrac]['start']  = 0.0
binLims[obsVarCldFrac]['finish'] = 1.0
binLims[obsVarCldFrac]['step']   = 0.2
binLims[obsVarCldFrac]['format'] = '{:.2f}'

binLims[obsVarCa] = {}
binLims[obsVarCa]['start']  = 0.0
binLims[obsVarCa]['finish'] = 60.0
binLims[obsVarCa]['step']   = 1.0
binLims[obsVarCa]['format'] = '{:.0f}'

#binLims[modVarLat] = {}
#binLims[modVarLat]['start']  = 90.0
#binLims[modVarLat]['finish'] = -90.0
#binLims[modVarLat]['step']   = -30.0
#binLims[modVarLat]['format'] = '{:.0f}'

#binLims[modVarAlt] = {}
#binLims[modVarAlt]['start']  = 0.0
#binLims[modVarAlt]['finish'] = 50000.0
#binLims[modVarAlt]['step']   = 5000.0
#binLims[modVarAlt]['format'] = '{:.0f}'

for binType, param in binLims.items():
    binBounds = list(np.arange(
        param['start']-0.5*np.abs(param['step']),
        param['finish']+0.5*param['step'],
        param['step']))
    binLims[binType]['minBounds'] = []
    binLims[binType]['maxBounds'] = []
    binLims[binType]['labels'] = []
    for ibin in list(range(len(binBounds)-1)):
        if param['step'] > 0.0:
            binLims[binType]['minBounds'].append(binBounds[ibin])
            binLims[binType]['maxBounds'].append(binBounds[ibin+1])
        else:
            binLims[binType]['minBounds'].append(binBounds[ibin+1])
            binLims[binType]['maxBounds'].append(binBounds[ibin])

        binVal = 0.5 * (binBounds[ibin+1] + binBounds[ibin])
        binLims[binType]['labels'].append(param['format'].format(binVal))
    binLims[binType]['nbins'] = len(binLims[binType]['labels'])

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


#========================
# binning functions
#========================

# default is to take first key from argDict
def firstKey(argDict):
    return argDict[list(argDict.keys())[0]]

def CaOkamoto(argDict):
    # Okamoto, et al.
    BTobs = argDict[selfObsValue]
    BTbak = argDict[selfBakValue]
    BTclr = deepcopy(argDict[clrskyBTDiag])
    BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
    Ca = np.multiply( 0.5,
             np.add(np.abs(np.subtract(BTobs,BTclr)),
                    np.abs(np.subtract(BTbak,BTclr))) )
    return Ca

def ScaledCaOkamoto(argDict):
    # Okamoto, et al.
    BTobs = argDict[selfObsValue]
    BTbak = argDict[selfBakValue]
    BTclr = deepcopy(argDict[clrskyBTDiag])
    BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
    CldFrac = argDict[cldfracMeta]

    Ca = np.multiply( 0.5,
             np.multiply(CldFrac,
                 np.add(np.abs(np.subtract(BTobs,BTclr)),
                        np.abs(np.subtract(BTbak,BTclr))) ) )
    return Ca

def CaModHarnisch(argDict):
    # Modified Harnisch, et al.
    BTobs = argDict[selfObsValue]
    BTbak = argDict[selfBakValue]
    BTclr = deepcopy(argDict[clrskyBTDiag])
    BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
    zeros = np.full_like(BTbak,0.0)
    Ca = np.multiply( 0.5,
             np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                    np.maximum(zeros,np.subtract(BTclr,BTbak))) )
    return Ca

def ScaledCaModHarnisch(argDict):
    # Modified Harnisch, et al.
    BTobs = argDict[selfObsValue]
    BTbak = argDict[selfBakValue]
    BTclr = deepcopy(argDict[clrskyBTDiag])
    BTclr[BTclr < 1.0] = BTbak[BTclr < 1.0]
    CldFrac = argDict[cldfracMeta]

    zeros = np.full_like(BTbak,0.0)
    Ca = np.multiply( 0.5,
             np.multiply(CldFrac,
                 np.add(np.maximum(zeros,np.subtract(BTclr,BTobs)),
                        np.maximum(zeros,np.subtract(BTclr,BTbak))) ) )
    return Ca

#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(argDict,REGION_NAME):
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = argDict['longitdue']
#    lats = argDict['latitude']
#    nlocs = len(lons)
#    REGIONS = isinsideregions(lons,lats,shp)
#    return REGIONS

#========================================
# bin dictionary used for all bin groups
#========================================

# placeholder binMethod definition
nullBinMethod = { 'filters': [], 'labels': [miss_s] }

defaultBinMethod = 'default'

# Each binGrp has the following properties
# variables: list of NC variables used by the binning filters
# binVarName: binning variable name that will be used for plotting
# binMethod (e.g., defaultBinMethod, NAMED, etc.): used to distinguish
#  between multiple methods in the same group
#     filters: list of filters that will be applied for each method
#     for each filter:
#         where: logical function the determines locations that are blacklisted
#         binArgs: sub-set of variables that is used in this where test
#         binFunc (optional): function of binArgs.  The where test is applied to the
#                             output of binFunc.  defaults to firstKey.
#         bounds: a list of numerical values used in the where test
#         apply_to (optional): list of diagnostics for which a particular filter applies
#     labels: list of labels (same length as bounds).  can be strings, integers, or floats

binGrpDict = {
#    #name            binVars              binVarName ... followed by dictionary defitions of each binDesc
    'ObsQC':        { 'variables': [selfQCValue]
                    , 'binVarName': obsVarQC
                    , defaultBinMethod:{
                          'filters': [
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag],
                           'apply_to': nonObsDiags}],
                          'labels': [goodFlagName] }
                    , 'bad':{
                          'filters': [
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': badFlags,
                           'apply_to': nonObsDiags},
                          {'where': equalBound,
                           'binArgs': [selfQCValue],
                           'bounds': badFlags,
                           'apply_to': nonObsDiags,
                           'mask_value': 0.0}],
                          'labels': badFlagNames }
                    }
  , 'ObsPressure':  { 'variables': [selfQCValue,prsMeta]
                    , 'binVarName': obsVarPrs
                    , defaultBinMethod:{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [prsMeta],
                           'bounds': binLims[obsVarPrs]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [prsMeta],
                           'bounds': binLims[obsVarPrs]['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarPrs]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarPrs]['labels'] }
                    }
  , 'ObsAltitude':  { 'variables': [selfQCValue,altMeta]
                    , 'binVarName': obsVarAlt
                    , defaultBinMethod:{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [altMeta],
                           'bounds': binLims[obsVarAlt]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [altMeta],
                           'bounds': binLims[obsVarAlt]['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarAlt]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarAlt]['labels'] }
                    }
  , 'ObsLatBand':   { 'variables': [selfQCValue,latMeta,prsMeta,altMeta,cldfracMeta]
                    , 'binVarName': obsVarLat
                    , defaultBinMethod:{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarLat]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarLat]['labels'] }
                    , 'NAMED':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': namedLatBands['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': namedLatBands['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*namedLatBands['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': namedLatBands['labels'] }
                    , 'P=300hPa':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['maxBounds']},
                          {'where': lessBound,
                           'binArgs': [prsMeta],
                           'bounds': [250.0]*binLims[obsVarLat]['nbins']},
                          {'where': greatEqualBound,
                           'binArgs': [prsMeta],
                           'bounds': [350.0]*binLims[obsVarLat]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarLat]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarLat]['labels'] }
                    , 'alt=10km':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['maxBounds']},
                          {'where': lessBound,
                           'binArgs': [altMeta],
                           'bounds': [9500.0]*binLims[obsVarLat]['nbins']},
                          {'where': greatEqualBound,
                           'binArgs': [altMeta],
                           'bounds': [10500.0]*binLims[obsVarLat]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarLat]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarLat]['labels'] }
                    , 'clear-sky':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['maxBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [cldfracMeta],
                           'bounds': [0.2]*binLims[obsVarLat]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarLat]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarLat]['labels'] }
                    , 'cloud-sky':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': binLims[obsVarLat]['maxBounds']},
                          {'where': lessBound,
                           'binArgs': [cldfracMeta],
                           'bounds': [0.8]*binLims[obsVarLat]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarLat]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarLat]['labels'] }
                    }
  , 'ObsZenith':    { 'variables': [selfQCValue,satzenMeta]
                    , 'binVarName': obsVarSatZen
                    , defaultBinMethod:{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [satzenMeta],
                           'bounds': binLims[obsVarSatZen]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [satzenMeta],
                           'bounds': binLims[obsVarSatZen]['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarSatZen]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarSatZen]['labels'] }
                    }
  , 'ObsCldFrac':   { 'variables': [selfQCValue,cldfracMeta]
                    , 'binVarName': obsVarCldFrac
                    , defaultBinMethod:{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [cldfracMeta],
                           'bounds': binLims[obsVarCldFrac]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [cldfracMeta],
                           'bounds': binLims[obsVarCldFrac]['maxBounds']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarCldFrac]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarCldFrac]['labels'] }
                    }
  , 'ObsCa':        { 'variables': [selfObsValue,selfBakValue,clrskyBTDiag,cldfracMeta]
                    , 'binVarName': obsVarCa
                    , 'Okamoto':{
                          'filters': [
                          {'where': lessBound,
                           'binFunc': CaOkamoto,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag],
                           'bounds': binLims[obsVarCa]['minBounds']},
                          {'where': greatEqualBound,
                           'binFunc': CaOkamoto,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag],
                           'bounds': binLims[obsVarCa]['maxBounds']}],
                          'labels': binLims[obsVarCa]['labels'] }
                    , 'ScaledOkamoto':{
                          'filters': [
                          {'where': lessBound,
                           'binFunc': ScaledCaOkamoto,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag,cldfracMeta],
                           'bounds': binLims[obsVarCa]['minBounds']},
                          {'where': greatEqualBound,
                           'binFunc': ScaledCaOkamoto,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag,cldfracMeta],
                           'bounds': binLims[obsVarCa]['maxBounds']}],
                          'labels': binLims[obsVarCa]['labels'] }
                    , 'ModHarnisch':{
                          'filters': [
                          {'where': lessBound,
                           'binFunc': CaModHarnisch,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag],
                           'bounds': binLims[obsVarCa]['minBounds']},
                          {'where': greatEqualBound,
                           'binFunc': CaModHarnisch,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag],
                           'bounds': binLims[obsVarCa]['maxBounds']}],
                          'labels': binLims[obsVarCa]['labels'] }
                    , 'ScaledModHarnisch':{
                          'filters': [
                          {'where': lessBound,
                           'binFunc': ScaledCaModHarnisch,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag,cldfracMeta],
                           'bounds': binLims[obsVarCa]['minBounds']},
                          {'where': greatEqualBound,
                           'binFunc': ScaledCaModHarnisch,
                           'binArgs': [selfObsValue,selfBakValue,clrskyBTDiag,cldfracMeta],
                           'bounds': binLims[obsVarCa]['maxBounds']}],
                          'labels': binLims[obsVarCa]['labels'] }
                    }
  , 'ObsBox':       { 'variables': [selfQCValue,lonMeta,latMeta]
                    , 'binVarName': 'REGION'
                    , 'AFRICA':    nullBinMethod
                    , 'ATLANTIC':  nullBinMethod
                    , 'AUSTRALIA': nullBinMethod
                    , 'CONUS':{
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [lonMeta],
                           'bounds': [234.0]},
                          {'where': greatBound,
                           'binArgs': [lonMeta],
                           'bounds': [294.0]},
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': [ 25.0]},
                          {'where': greatBound,
                           'binArgs': [latMeta],
                           'bounds': [ 50.0]},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag],
                           'apply_to': nonObsDiags} ],
                          'labels': ['CONUS'] }
                    , 'EUROPE':    nullBinMethod
                    , 'E_EUROPE':  nullBinMethod
                    , 'NAMERICA':  nullBinMethod
                    , 'PACIFIC':   nullBinMethod
                    , 'SAMERICA':  nullBinMethod
                    , 'SE_ASIA':   nullBinMethod
                    , 'S_ASIA':    nullBinMethod
                    }
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#  , 'polygon':   { 'variables': [selfQCValue,lonMeta,latMeta]
#                 , 'binVarName': 'REGION'
#                 , 'CONUS':{
#                       'filters': [
#                       {'where': notInBound,
#                        'binFunc': Regions,
#                        'binArgs': [lonMeta,latMeta],
#                        'bounds': ['CONUS']}],
#                       'labels': ['CONUS'] }
#                 , 'EUROPE':{
#                       'filters': [
#                       {'where': notInBound,
#                        'binFunc': Regions,
#                        'binArgs': [lonMeta,latMeta],
#                        'bounds': ['EUROPE']}],
#                       'labels': ['EUROPE'] }
#                  }
#TODO: enable binning in MPAS Model space
#  , 'ModelPressure':  { 'variables': [selfQCValue,modVarPrs]
#                      , 'binVarName': modVarPrs
#                      , defaultBinMethod:{
#                            'filters': [
#                            {'where': lessBound,
#                             'binArgs': [modVarPrs],
#                             'bounds': binLims[modVarPrs]['minBounds']},
#                            {'where': greatEqualBound,
#                             'binArgs': [modVarPrs],
#                             'bounds': binLims[modVarPrs]['maxBounds']} ],
#                            'labels': binLims[modVarPrs]['labels'] }
#                      }
#  , 'ModelAltitude':  { 'variables': [modVarAlt]
#                      , 'binVarName': modVarAlt
#                      , defaultBinMethod:{
#                            'filters': [
#                            {'where': lessBound,
#                             'binArgs': [modVarAlt],
#                             'bounds': binLims[modVarAlt]['minBounds']},
#                            {'where': greatEqualBound,
#                             'binArgs': [modVarAlt],
#                             'bounds': binLims[modVarAlt]['maxBounds']} ],
#                            'labels': binLims[modVarAlt]['labels'] }
#                      }
#  , 'ModelLatBand':   { 'variables': [modVarLat]
#                      , 'binVarName': modVarLat
#                      , defaultBinMethod:{
#                            'filters': [
#                            {'where': lessBound,
#                             'binArgs': [modVarLat],
#                             'bounds': binLims[modVarLat]['minBounds']},
#                            {'where': greatEqualBound,
#                             'binArgs': [modVarLat],
#                             'bounds': binLims[modVarLat]['maxBounds']} ],
#                            'labels: binLims[modVarLat]['labels'] }
#                      , 'NAMED':      {
#                            'filters': [
#                            {'where': lessBound,
#                             'binArgs': [modVarLat],
#                             'bounds': namedLatBands['minBounds']},
#                            {'where': greatEqualBound,
#                             'binArgs': [modVarLat],
#                             'bounds': namedLatBands['maxBounds']},
#                            'labels': namedLatBands['labels'] }
#                      }
#  , 'ModelBox':       { 'variables': [modVarLon,modVarLat]
#                      , 'binVarName': 'REGION'
#                      , 'AFRICA':    nullBinMethod
#                      , 'ATLANTIC':  nullBinMethod
#                      , 'AUSTRALIA': nullBinMethod
#                      , 'CONUS':      {
#                            'filters': [
#                            {'where': lessBound,
#                             'binArgs': [modVarLon],
#                             'bounds': [234.0]},
#                            {'where': greatBound,
#                             'binArgs': [modVarLon],
#                             'bounds': [294.0]},
#                            {'where': lessBound,
#                             'binArgs': [modVarLat],
#                             'bounds': [ 25.0]},
#                            {'where': greatBound,
#                             'binArgs': [modVarLat],
#                             'bounds': [ 50.0]} ],
#                            'labels': ['CONUS'] }
#                      , 'EUROPE':    nullBinMethod
#                      , 'E_EUROPE':  nullBinMethod
#                      , 'NAMERICA':  nullBinMethod
#                      , 'PACIFIC':   nullBinMethod
#                      , 'SAMERICA':  nullBinMethod
#                      , 'SE_ASIA':   nullBinMethod
#                      , 'S_ASIA':    nullBinMethod
#                       }
}


#==================================================
# Sub-selections of binGrps for specific DiagSpaces
#==================================================

nullBinGrps = {miss_s: []}

## Generic binGrps that apply to all observation categories
obsBinGrps = { 'ObsQC':      [defaultBinMethod,'bad']
             , 'ObsLatBand': [defaultBinMethod,'NAMED']
             , 'ObsBox':     ['CONUS']
             }


## binGrps for surface obs
surfBinGrps = deepcopy(obsBinGrps)


## binGrps for profile obs w/ pressure vertical bins
profPressBinGrps = deepcopy(obsBinGrps)
profPressBinGrps['ObsPressure'] = [defaultBinMethod]
profPressBinGrps['ObsLatBand'].append('P=300hPa')

# Add named latitude-band-specific pressure bins
if latMeta not in binGrpDict['ObsPressure']['variables']:
    binGrpDict['ObsPressure']['variables'].append(latMeta)
for iband, latBand in enumerate(namedLatBands['labels']):
    profPressBinGrps['ObsPressure'].append(latBand)
    binGrpDict['ObsPressure'][latBand] = {
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [prsMeta],
                           'bounds': binLims[obsVarPrs]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [prsMeta],
                           'bounds': binLims[obsVarPrs]['maxBounds']},
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': [namedLatBands['minBounds'][iband]]*binLims[obsVarPrs]['nbins']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': [namedLatBands['maxBounds'][iband]]*binLims[obsVarPrs]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarPrs]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarPrs]['labels'] }


## binGrps for profile obs w/ altitude vertical bins
profAltBinGrps = deepcopy(obsBinGrps)
profAltBinGrps['ObsAltitude'] = [defaultBinMethod]
profAltBinGrps['ObsLatBand'].append('alt=10km')

# Add named latitude-band-specific altitude bins
if latMeta not in binGrpDict['ObsAltitude']['variables']:
    binGrpDict['ObsAltitude']['variables'].append(latMeta)
for iband, latBand in enumerate(namedLatBands['labels']):
    profAltBinGrps['ObsAltitude'].append(latBand)
    binGrpDict['ObsAltitude'][latBand] = {
                          'filters': [
                          {'where': lessBound,
                           'binArgs': [altMeta],
                           'bounds': binLims[obsVarAlt]['minBounds']},
                          {'where': greatEqualBound,
                           'binArgs': [altMeta],
                           'bounds': binLims[obsVarAlt]['maxBounds']},
                          {'where': lessBound,
                           'binArgs': [latMeta],
                           'bounds': [namedLatBands['minBounds'][iband]]*binLims[obsVarAlt]['nbins']},
                          {'where': greatEqualBound,
                           'binArgs': [latMeta],
                           'bounds': [namedLatBands['maxBounds'][iband]]*binLims[obsVarAlt]['nbins']},
                          {'where': notEqualBound,
                           'binArgs': [selfQCValue],
                           'bounds': [goodFlag]*binLims[obsVarAlt]['nbins'],
                           'apply_to': nonObsDiags} ],
                          'labels': binLims[obsVarAlt]['labels'] }


## binGrps for radiance obs
radianceBinGrps = deepcopy(obsBinGrps)
radianceBinGrps['ObsZenith'] = [defaultBinMethod]


## binGrps for GOES-ABI obs
abiBinGrps = deepcopy(radianceBinGrps)
abiBinGrps['ObsCldFrac'] = [defaultBinMethod]
abiBinGrps['ObsLatBand'].append('clear-sky')
abiBinGrps['ObsLatBand'].append('cloud-sky')
abiBinGrps['ObsCa'] = ['Okamoto','ScaledOkamoto']
#abiBinGrps['ObsCa'].append('ModHarnisch')
#abiBinGrps['ObsCa'].append('ScaledModHarnisch')

#modelBinGrps = { 'ModelLatBand':  ['NAMED',defaultBinMethod]
#               , 'ModelBox':      ['CONUS']
#               , 'ModelAltitude': [defaultBinMethod]
#               , 'ModelPressure': [defaultBinMethod]
#               }

#=======================
# DiagSpace definitions
# e.g. IODA ObsSpace
#      MPAS ModelSpace
#=======================

nullDiagSpaceInfo = {'DiagSpaceGrp': miss_s,    'process': False, 'binGrps': nullBinGrps }

profile_s = 'profile'
radiance_s = 'radiance'
model_s = 'MPAS'

# columns: DiagSpace name (YAML)    DiagSpaceGrp              process?                  binGrps
DiagSpaceDict = {
    'sondes':                {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps }
  , 'aircraft':              {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps }
  , 'satwind':               {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps }
  , 'gnssroref':             {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   }
  , 'gnssrobndropp1d':       {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   }
  , 'gnssro':                {'DiagSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   }
  , 'abi_g16':               {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': abiBinGrps,
                              'channels': [8,9,10,11,13,14,15,16] }
  , 'airs_aqua':             {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': [1,6,7] }
  , 'amsua_n15':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [5,6,7,8,9] }
  , 'amsua_n18':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [5,6,7,8,9] }
  , 'amsua_n19':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [5,6,7,9] }
  , 'amsua_metop-a':         {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [5,6,9] }
  , 'amsua_metop-b':         {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [] }
  , 'amsua_aqua':            {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [8,9] }
  , 'amsua_n19--hydro':      {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [1,2,3,15] }
  , 'amsua_n19--nohydro':    {'DiagSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps,
                              'channels': [4,5,6,7,9,10,11,12,13,14] }
  , 'cris-fsr_npp':          {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': [24,26,28,32,37,39] }
  , 'hirs4_metop-a':         {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,16) }
  , 'iasi_metop-a':          {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': [16,29,32,35,38,41,44] }
  , 'mhs_n19':               {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,6) }
  , 'seviri_m08':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': [5] }
  , 'sndrd1_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,16) }
  , 'sndrd2_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,16) }
  , 'sndrd3_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,16) }
  , 'sndrd4_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps,
                              'channels': range(1,16) }
#  , 'mpas-gfs':              {'DiagSpaceGrp': model_s,    'process': False, 'binGrps': modelBinGrps    }
    }

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
