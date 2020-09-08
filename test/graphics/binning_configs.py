#!/usr/bin/env python3

import binning_utils as bu
import diag_utils as du
import numpy as np
import var_utils as vu

#==========================
# names and values of bins
#==========================

## heterogeneous/named bins
# note: the order of subplots will follow what is specified here

# latbandsMethod bins, north to south
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


# cloudbandsMethod bins
clrskyThresh = 0.05
cldskyThresh = 1.0 - clrskyThresh

namedCldFracBands = {}
namedCldFracBands['values'] = [bu.clrskyMethod, bu.mixskyMethod, bu.cldskyMethod, bu.allskyMethod]
namedCldFracBands['minBounds'] = [0.0, clrskyThresh, cldskyThresh, 0.0]
namedCldFracBands['maxBounds'] = [clrskyThresh, cldskyThresh, 1.0, 1.0]


# surfbandsMethod bins
seasurfThresh = 0.05
landsurfThresh = 1.0 - seasurfThresh

namedLandFracBands = {}
namedLandFracBands['values'] = [bu.seasurfMethod, bu.mixsurfMethod, bu.landsurfMethod, bu.allsurfMethod]
namedLandFracBands['minBounds'] = [0.0, seasurfThresh, landsurfThresh, 0.0]
namedLandFracBands['maxBounds'] = [seasurfThresh, landsurfThresh, 1.0, 1.0]


## homogeneous bins
# note: the ordering described here does not make any difference
#       figure axes will be monotonically increasing except for pressure

# binLims is used to auto-generate a large number of binning
# configurations in the binVarConfigs dictionary

binLims = {}

binLims[vu.obsVarPrs] = {}
binLims[vu.obsVarPrs]['start']  =    0.0
binLims[vu.obsVarPrs]['finish'] = 1000.0
binLims[vu.obsVarPrs]['step']   =  100.0
binLims[vu.obsVarPrs]['format'] = '{:.0f}'

binLims[vu.obsVarAlt] = {}
binLims[vu.obsVarAlt]['start']  = 1000.0
binLims[vu.obsVarAlt]['finish'] = 50000.0
binLims[vu.obsVarAlt]['step']   = 2000.0
binLims[vu.obsVarAlt]['format'] = '{:.0f}'

binLims[vu.obsVarLat] = {}
binLims[vu.obsVarLat]['start']  = -90.0
binLims[vu.obsVarLat]['finish'] =  90.0
binLims[vu.obsVarLat]['step']   =  10.0
binLims[vu.obsVarLat]['format'] = '{:.0f}'

binLims[vu.obsVarLT] = {}
binLims[vu.obsVarLT]['start']  = bu.LH0
binLims[vu.obsVarLT]['finish'] = bu.LH1
binLims[vu.obsVarLT]['step']   = bu.LHDT
binLims[vu.obsVarLT]['format'] = '{:.0f}'

binLims[vu.obsVarSenZen] = {}
binLims[vu.obsVarSenZen]['start']  = 0.0
binLims[vu.obsVarSenZen]['finish'] = 70.0
binLims[vu.obsVarSenZen]['step']   = 5.0
binLims[vu.obsVarSenZen]['format'] = '{:.0f}'

binLims[vu.obsVarGlint] = {}
binLims[vu.obsVarGlint]['start']  = 0.0
binLims[vu.obsVarGlint]['finish'] = bu.maxGlint
binLims[vu.obsVarGlint]['step']   = 10.0
binLims[vu.obsVarGlint]['format'] = '{:.0f}'

binLims[vu.obsVarLandFrac] = {}
binLims[vu.obsVarLandFrac]['start']  = 0.0
binLims[vu.obsVarLandFrac]['finish'] = 1.0
binLims[vu.obsVarLandFrac]['step']   = seasurfThresh / 2.0
binLims[vu.obsVarLandFrac]['format'] = '{:.3f}'

binLims[vu.obsVarCldFrac] = {}
binLims[vu.obsVarCldFrac]['start']  = 0.0
binLims[vu.obsVarCldFrac]['finish'] = 1.0
binLims[vu.obsVarCldFrac]['step']   = clrskyThresh / 2.0
binLims[vu.obsVarCldFrac]['format'] = '{:.3f}'

binLims[vu.obsVarSCI] = {}
binLims[vu.obsVarSCI]['start']  = 0.0
binLims[vu.obsVarSCI]['finish'] = 60.0
binLims[vu.obsVarSCI]['step']   = 1.0
binLims[vu.obsVarSCI]['format'] = '{:.0f}'

binLims[vu.obsVarACI] = {}
binLims[vu.obsVarACI]['start']  = -20.0
binLims[vu.obsVarACI]['finish'] = 20.0
binLims[vu.obsVarACI]['step']   = 2.0
binLims[vu.obsVarACI]['format'] = '{:.0f}'

binLims[vu.obsVarNormErr] = {}
binLims[vu.obsVarNormErr]['start']  = -7.0
binLims[vu.obsVarNormErr]['finish'] =  7.0
binLims[vu.obsVarNormErr]['step']   =  0.25
binLims[vu.obsVarNormErr]['format'] = '{:.2f}'

#binLims[vu.modVarLat] = {}
#binLims[vu.modVarLat]['start']  = -90.0
#binLims[vu.modVarLat]['finish'] =  90.0
#binLims[vu.modVarLat]['step']   =  30.0
#binLims[vu.modVarLat]['format'] = '{:.0f}'

#binLims[vu.modVarAlt] = {}
#binLims[vu.modVarAlt]['start']  = 0.0
#binLims[vu.modVarAlt]['finish'] = 50000.0
#binLims[vu.modVarAlt]['step']   = 5000.0
#binLims[vu.modVarAlt]['format'] = '{:.0f}'

for binType, param in binLims.items():
    binBounds = list(np.arange(
        param['start']-0.5*np.abs(param['step']),
        param['finish']+1.5*param['step'],
        param['step']))
    binLims[binType]['minBounds'] = []
    binLims[binType]['maxBounds'] = []
    binLims[binType]['values'] = []
    for ibin in list(range(len(binBounds)-1)):
        binLims[binType]['minBounds'].append(binBounds[ibin])
        binLims[binType]['maxBounds'].append(binBounds[ibin+1])

        binVal = 0.5 * (binBounds[ibin+1] + binBounds[ibin])
        binLims[binType]['values'].append(param['format'].format(binVal))


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


#=================================================
# binning configurations used for all bin methods
#=================================================

# Each binVarConfig member has the following properties
# key: binVar string describes the variable to bin over
# binMethod: used to distinguish between multiple methods with the same binVar
#  (e.g., bu.identityBinMethod, bu.badQC, bu.latbandsMethod, etc.)
#     filters: list of filters that will blacklist locations for each method
#     for each filter:
#         where: logical function that determines locations that are blacklisted
#         variable: string or class that is used to initialize the IdObsFunction or ObsFunction class
#         bounds: numerical value(s) used in the where test (scalar or Iterable same length as values). At least one filter must have len(bounds)==len(values).
#         except_diags (optional): list of diagnostics for which a particular filter does not apply
#     values (string): list of values associated with each bin
#
#TODO: classify each binmethod as 1D, 2D, etc. to indicate which types of figures apply to it

nullBinMethod = { 'filters': [], 'values': [] }
nullBinVarConfig = { bu.identityBinMethod:nullBinMethod }

binVarConfigs = {
    vu.obsVarQC: {
        bu.goodQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': goodFlagName,
        },
        bu.badQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': badFlags,
                 'except_diags': du.nonQCedDiags},
                {'where': bu.equalBound,
                 'variable': vu.selfQCValue,
                 'bounds': badFlags,
                 'except_diags': du.nonQCedDiags,
                 'mask_value': 0.0},
            ],
            'values': badFlagNames,
        },
    },
    vu.obsVarPrs: {
        bu.PjetMethod: {
            'filters': [
# eliminate locations outside bu.P_jet_min to bu.P_jet_max
                {'where': bu.lessBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_max},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': bu.P_jet_val,
        },
    },
    vu.obsVarAlt: {
        bu.altjetMethod: {
            'filters': [
# eliminate locations outside bu.alt_jet_min to bu.alt_jet_max
                {'where': bu.lessBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_max},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': bu.alt_jet_val,
        },
    },
    vu.obsVarLat: {
        bu.latbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['maxBounds']},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': namedLatBands['values'],
        },
        bu.PjetMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': binLims[vu.obsVarLat]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': binLims[vu.obsVarLat]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_max},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
        bu.altjetMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': binLims[vu.obsVarLat]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': binLims[vu.obsVarLat]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_max},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[vu.obsVarLat]['values'],
        },
    },
    vu.obsVarLandFrac: {
        bu.surfbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['minBounds']},
                {'where': bu.greatBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['maxBounds']},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': namedLandFracBands['values'],
        },
    },
    vu.obsVarCldFrac: {
        bu.cloudbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['minBounds']},
                {'where': bu.greatBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['maxBounds']},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': namedCldFracBands['values'],
        },
    },
    vu.obsVarSCI: {
        bu.OkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.SCIOkamoto,
                 'bounds': binLims[vu.obsVarSCI]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.SCIOkamoto,
                 'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        bu.ScaleOkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ScaledSCIOkamoto,
                 'bounds': binLims[vu.obsVarSCI]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.ScaledSCIOkamoto,
                 'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        bu.ModHarnischMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.SCIModHarnisch,
                 'bounds': binLims[vu.obsVarSCI]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.SCIModHarnisch,
                 'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
        bu.ScaleModHarnischMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ScaledSCIModHarnisch,
                 'bounds': binLims[vu.obsVarSCI]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.ScaledSCIModHarnisch,
                 'bounds': binLims[vu.obsVarSCI]['maxBounds']},
            ],
            'values': binLims[vu.obsVarSCI]['values'],
        },
    },
    vu.obsVarNormErr: {
        bu.OkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.OkamotoNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.OkamotoNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        bu.ScaleOkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ScaledOkamotoNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.ScaledOkamotoNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        bu.ModHarnischMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ModHarnischNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.ModHarnischNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
        bu.ScaleModHarnischMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ScaledModHarnischNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': bu.ScaledModHarnischNormalizedError,
                 'bounds': binLims[vu.obsVarNormErr]['maxBounds']},
            ],
            'values': binLims[vu.obsVarNormErr]['values'],
        },
    },
    'ObsRegion': {
        'AFRICA': nullBinMethod,
        'ATLANTIC': nullBinMethod,
        'AUSTRALIA': nullBinMethod,
        'CONUS': {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.lonMeta,
                 'bounds': [234.0]},
                {'where': bu.greatBound,
                 'variable': vu.lonMeta,
                 'bounds': [294.0]},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': [ 25.0]},
                {'where': bu.greatBound,
                 'variable': vu.latMeta,
                 'bounds': [ 50.0]},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': ['CONUS'],
        },
        'EUROPE': nullBinMethod,
        'E_EUROPE': nullBinMethod,
        'NAMERICA': nullBinMethod,
        'PACIFIC': nullBinMethod,
        'SAMERICA': nullBinMethod,
        'SE_ASIA': nullBinMethod,
        'S_ASIA': nullBinMethod,
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#        'CONUS_POLYGON': {
#            'filters': [
#                {'where': bu.notInBound,
#                 'variable': bu.Regions,
#                 'variable': [vu.lonMeta, vu.latMeta],
#                 'bounds': ['CONUS']},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE_POLYGON': {
#            'filters': [
#                {'where': bu.notInBound,
#                 'variable': bu.Regions,
#                 'variable': [vu.lonMeta, vu.latMeta],
#                 'bounds': ['EUROPE']},
#            ],
#            'values': ['EUROPE'],
#        },
    },
#TODO: enable binning in MPAS Model space
#    vu.modVarPrs: {
#        bu.identityBinMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarPrs,
#                 'bounds': binLims[vu.modVarPrs]['minBounds']},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarPrs,
#                 'bounds': binLims[vu.modVarPrs]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarPrs]['values'],
#        },
#    },
#    vu.modVarAlt: {
#        bu.identityBinMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarAlt,
#                 'bounds': binLims[vu.modVarAlt]['minBounds']},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarAlt,
#                 'bounds': binLims[vu.modVarAlt]['maxBounds']},
#            ],
#            'values': binLims[vu.modVarAlt]['values'],
#        },
#    },
#    vu.modVarLat: {
#        bu.identityBinMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLat,
#                 'bounds': binLims[vu.modVarLat]['minBounds']},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarLat,
#                 'bounds': binLims[vu.modVarLat]['maxBounds']},
#            ],
#            'values: binLims[vu.modVarLat]['values'],
#        },
#        bu.latbandsMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLat,
#                 'bounds': namedLatBands['minBounds']},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarLat,
#                 'bounds': namedLatBands['maxBounds']},
#            ],
#            'values': namedLatBands['values'],
#        },
#    },
#    'ModelRegion': {
#        'AFRICA': nullBinMethod,
#        'ATLANTIC': nullBinMethod,
#        'AUSTRALIA': nullBinMethod,
#        'CONUS': {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLon,
#                 'bounds': [234.0]},
#                {'where': bu.greatBound,
#                 'variable': vu.modVarLon,
#                 'bounds': [294.0]},
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLat,
#                 'bounds': [ 25.0]},
#                {'where': bu.greatBound,
#                 'variable': vu.modVarLat,
#                 'bounds': [ 50.0]},
#            ],
#            'values': ['CONUS'],
#        },
#        'EUROPE': nullBinMethod,
#        'E_EUROPE': nullBinMethod,
#        'NAMERICA': nullBinMethod,
#        'PACIFIC': nullBinMethod,
#        'SAMERICA': nullBinMethod,
#        'SE_ASIA': nullBinMethod,
#        'S_ASIA': nullBinMethod,
#    },
}


#=============================
# Parameterized binVarConfigs
#=============================

# Add bu.identityBinMethod for identity ranged binning variables
identityRangeBinVars = {
    vu.obsVarAlt: [vu.altMeta, []],
    vu.obsVarACI: [bu.AsymmetricCloudImpact, ['obs','bak','ana','SCI']],
    vu.obsVarCldFrac: [vu.cldfracMeta, ['obs','bak','ana','SCI']],
    vu.obsVarGlint: [bu.GlintAngle, ['obs','bak','ana','SCI']],
    vu.obsVarLandFrac: [vu.landfracGeo, ['obs','bak','ana','SCI']],
    vu.obsVarLat: [vu.latMeta, ['obs','bak','ana','SCI']],
    vu.obsVarLT: [bu.LocalHour, ['obs','bak','ana','SCI']],
    vu.obsVarNormErr: [bu.NormalizedError, []],
    vu.obsVarPrs: [vu.prsMeta, []],
    vu.obsVarSenZen: [vu.senzenMeta, ['obs','bak','ana','SCI']],
}
for binVar, rangeVar in identityRangeBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    binVarConfigs[binVar][bu.identityBinMethod] = {
        'filters': [
            {'where': bu.lessBound,
             'variable': rangeVar[0],
             'bounds': binLims[binVar]['minBounds']},
            {'where': bu.greatEqualBound,
             'variable': rangeVar[0],
             'bounds': binLims[binVar]['maxBounds']},
            {'where': bu.notEqualBound,
             'variable': vu.selfQCValue,
             'bounds': goodFlag,
             'except_diags': du.nonQCedDiags},
        ],
        'values': binLims[binVar]['values'],
        'override_exclusiveDiags': rangeVar[1],
    }

# Add named latitude-band-specific bins for applicable ranged variables
latBinVars = {
    vu.obsVarAlt: [vu.altMeta, []],
    vu.obsVarPrs: [vu.prsMeta, []],
}
for binVar, rangeVar in latBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for iband, latBand in enumerate(namedLatBands['values']):
        binVarConfigs[binVar][latBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['minBounds'][iband]},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['maxBounds'][iband]},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': rangeVar[1],
        }

# Add named latitude-band-specific bins for clear-sky ranged variables
clrlatBinVars = {
#    vu.obsVarGlint: [bu.GlintAngle, []],
    vu.obsVarSenZen: [vu.senzenMeta, []],
}
clrband = namedCldFracBands['values'].index(bu.clrskyMethod)
clrlatMethods = {}
for binVar, rangeVar in clrlatBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for iband, latBand in enumerate(namedLatBands['values']):
        clrlatMethods[latBand] = latBand+'-'+bu.clrskyMethod
        binVarConfigs[binVar][clrlatMethods[latBand]] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['minBounds'][iband]},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['maxBounds'][iband]},
                {'where': bu.lessBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['minBounds'][clrband]},
                {'where': bu.greatBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['maxBounds'][clrband]},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': rangeVar[1],
        }

# Add named cloud fraction-band-specific bins for applicable ranged variables
cldfracBinVars = {
    vu.obsVarACI: [bu.AsymmetricCloudImpact, ['obs','bak','ana','SCI']],
    vu.obsVarGlint: [bu.GlintAngle, ['obs','bak','ana','SCI']],
    vu.obsVarLandFrac: [vu.landfracGeo, ['obs','bak','ana','SCI']],
    vu.obsVarLat: [vu.latMeta, ['obs','bak','ana','SCI']],
    vu.obsVarLT: [bu.LocalHour, ['obs','bak','ana','SCI']],
    vu.obsVarSenZen: [vu.senzenMeta, ['obs','bak','ana','SCI']],
}
for binVar, rangeVar in cldfracBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for iband, cldBand in enumerate(namedCldFracBands['values']):
        binVarConfigs[binVar][cldBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['minBounds'][iband]},
                {'where': bu.greatBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['maxBounds'][iband]},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': rangeVar[1],
    }

# Add named land fraction-band-specific bins for applicable ranged variables
landfracBinVars = {
    vu.obsVarSCI: [bu.SCIOkamoto, ['obs','bak','ana','SCI']],

}
for binVar, rangeVar in landfracBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for iband, surfBand in enumerate(namedLandFracBands['values']):
        binVarConfigs[binVar][surfBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['minBounds']},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binLims[binVar]['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['minBounds'][iband]},
                {'where': bu.greatBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['maxBounds'][iband]},
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlag,
                 'except_diags': du.nonQCedDiags},
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': rangeVar[1],
    }
