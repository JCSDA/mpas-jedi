#!/usr/bin/env python3

import binning_utils as bu
from collections import defaultdict
from copy import deepcopy
from jediApplicationArgs import jediAppName, nOuterIter
import numpy as np
import os
import var_utils as vu

#======================================================
# outer iteration settings for Variational applications
#======================================================

outerIter = str(nOuterIter)
anIter = str(nOuterIter)
appName = jediAppName


#=============================================================
# groups of diagnostics for which statistics can be calculated
#=============================================================

## model-space diagnostics
modelDiags = ['mmgfsan']

## difference diagnostics
diffDiagnostics_ = defaultdict(list)
diffDiagnostics_['variational'] += ['omb']
diffDiagnostics_['hofx'] += ['omf']

## relative difference diagnostics
relativeDiagnostics_ = defaultdict(list)
relativeDiagnostics_['variational'] += ['rltv_omb']
relativeDiagnostics_['hofx'] += ['rltv_omf']

## absolute diagnostics
absoluteDiagnostics_ = defaultdict(list)
absoluteDiagnostics_['variational'] += ['obs', 'bak']
absoluteDiagnostics_['hofx'] += ['obs', 'h(x)']

## cloudy radiance diagnostics
cloudyRadDiagnostics = ['SCI']

## STD diagnostics
sigmaDiagnostics_ = defaultdict(list)
sigmaDiagnostics_['variational'] = ['omb','sigmaob', 'sigmab', 'CRyb']
sigmaDiagnostics_['hofx'] = ['omf','sigmaof', 'sigmaf', 'CRyf']

# variational analysis diagnostics
if nOuterIter > 0:
    diffDiagnostics_['variational'] += ['oma']
    relativeDiagnostics_['variational'] += ['rltv_oma']
    absoluteDiagnostics_['variational'] += ['ana']
    sigmaDiagnostics_['variational'] += ['oma','sigmaoa','sigmaa', 'CRya']

diffDiagnostics = diffDiagnostics_[jediAppName]
relativeDiagnostics = relativeDiagnostics_[jediAppName]
absoluteDiagnostics = absoluteDiagnostics_[jediAppName]
sigmaDiagnostics = sigmaDiagnostics_[jediAppName]

defaultDiagnostics = diffDiagnostics

## diagnostics for which QC is irrelevant
nonQCedDiagnostics = ['obs']

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


# geoirlatlonboxMethod bins
geoirLonBands = deepcopy(bu.geoirlatlonBoxParams)
geoirLatBands = deepcopy(bu.geoirlatlonBoxParams)

# store minBounds/maxBounds
for centerLon in geoirLonBands['centerLon']:
  geoirLonBands['minBounds'] += \
    [centerLon - bu.geoirMaxZenith]
  geoirLonBands['maxBounds'] += \
    [centerLon + bu.geoirMaxZenith]
  geoirLatBands['minBounds'] += \
    [-bu.geoirMaxZenith]
  geoirLatBands['maxBounds'] += \
    [bu.geoirMaxZenith]

# ensure positive-definite longitude
for bound in ['minBounds', 'maxBounds']:
  for ii, lon in enumerate(geoirLonBands[bound]):
    while geoirLonBands[bound][ii] > 360.:
      geoirLonBands[bound][ii] -= 360.
    while geoirLonBands[bound][ii] < 0.:
      geoirLonBands[bound][ii] += 360.

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

binLims[vu.modVarLat] = {}
binLims[vu.modVarLat]['start']  = -90.0
binLims[vu.modVarLat]['finish'] =  90.0
binLims[vu.modVarLat]['step']   =  5.0
binLims[vu.modVarLat]['format'] = '{:.0f}'

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
# passive = 1;   // H(x) is computed (for monitoring, BC...) but obs not assimilated
# missing = 10;  // missing values prevent use of observation
# preQC   = 11;  // observation rejected by pre-processing
# bounds  = 12;  // observation value out of bounds
# domain  = 13;  // observation not within domain of use
# exclude = 14;  // observation excluded
# Hfailed = 15;  // H(x) computation failed
# thinned = 16;  // observation removed due to thinning
# diffref = 17;  // metadata too far from reference
# clw     = 18;  // observation removed due to cloud field
# fguess  = 19;  // observation too far from guess
# seaice  = 20;  // observation based sea ice detection, also flags land points
# track   = 21;  // observation removed as inconsistent with the rest of track
# buddy   = 22;  // observation rejected by the buddy check
# derivative = 23;  // observation removed due to metadata derivative value
# profile = 24;  // observation rejected by at least one profile QC check
# onedvar = 25;  // observation failed to converge in 1dvar check
# bayesianQC = 26;  // observation failed due to Bayesian background check
# modelobthresh = 27;  // observation failed modelob threshold check
# Static list above copied on 22 Apr 2021
# see ufo/src/ufo/filters/QCflags.h for up-to-date list

goodFlags = [0, 1]
goodFlagNames = ['pass', 'passive']

badFlags = [10, 11, 12, 13,
            14, 15, 16, 17,
            18, 19, 20, 21,
            22, 23, 24, 25, 26, 27,
]
badFlagNames = ['missing', 'preQC',      'bounds',  'domain',
                'exclude', 'Hfailed',    'thinned', 'diffref',
                'clw',     'fguess',     'seaice', 'track',
                'buddy',   'derivative', 'profile', 'onedvar', 'bayesianQC', 'modelobthresh',
]

#=================================================
# binning configurations used for all bin methods
#=================================================

# Each binVarConfig member has the following properties
# key: binVar string describes the variable to bin over
# binMethod: used to distinguish between multiple methods with the same binVar
#  (e.g., bu.identityBinMethod, bu.badQC, bu.latbandsMethod, etc.)
#     filters: list of filters that will exclude locations for each method
#     for each filter:
#         where: logical function that determines locations that are excluded
#         variable: string or class that is used to initialize the IdObsFunction or ObsFunction class;
#                   this may refer to a predefined key of vu.ObsGroups and vu.ObsVars
#         bounds: numerical value(s) used in the where test (scalar or Iterable same length as values). At least one filter must have len(bounds)==len(values).
#         except_diags (optional): list of diagnostics for which a particular filter does not apply
#     values (string): list of values associated with each bin
#
#TODO: classify each binmethod as 1D, 2D, etc. to indicate which types of figures apply to it

nullBinMethod = { 'filters': [], 'values': [] }
nullBinVarConfig = { bu.identityBinMethod:nullBinMethod }

# filter locations that do not have goodFlags
AnyBadQC = {
  'where': bu.notEqualAnyBound,
  'variable': vu.selfQCValue,
  'bounds': [goodFlags],
  'except_diags': nonQCedDiagnostics,
}

binVarConfigs = {
    vu.obsVarQC: {
        bu.goodQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlags,
                 'except_diags': nonQCedDiagnostics}
            ],
            'values': goodFlagNames,
        },
        bu.badQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': badFlags,
                 'except_diags': nonQCedDiagnostics},
                {'where': bu.equalBound,
                 'variable': vu.selfQCValue,
                 'bounds': badFlags,
                 'except_diags': nonQCedDiagnostics,
                 'mask_value': 0.0},
            ],
            'values': badFlagNames,
        },
        bu.allQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlags+badFlags,
                 'except_diags': nonQCedDiagnostics},
                {'where': bu.equalBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlags+badFlags,
                 'except_diags': nonQCedDiagnostics,
                 'mask_value': 0.0},
            ],
            'values': goodFlagNames+badFlagNames,
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
                AnyBadQC,
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
                AnyBadQC,
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
                AnyBadQC,
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
                AnyBadQC,
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
                AnyBadQC,
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
                AnyBadQC,
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
                AnyBadQC,
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
    vu.obsRegionBinVar: {
        'AFRICA': nullBinMethod,
        'ATLANTIC': nullBinMethod,
        'AUSTRALIA': nullBinMethod,
        'CONUS': {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.lonMeta,
                 'bounds': 234.0},
                {'where': bu.greatBound,
                 'variable': vu.lonMeta,
                 'bounds': 294.0},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': 25.0},
                {'where': bu.greatBound,
                 'variable': vu.latMeta,
                 'bounds': 50.0},
                AnyBadQC,
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
        bu.geoirlatlonboxMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.lonMeta,
                 'bounds': geoirLonBands['minBounds']},
                {'where': bu.greatBound,
                 'variable': vu.lonMeta,
                 'bounds': geoirLonBands['maxBounds']},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': geoirLatBands['minBounds']},
                {'where': bu.greatBound,
                 'variable': vu.latMeta,
                 'bounds': geoirLatBands['maxBounds']},
            ],
            'values': geoirLonBands['values'],
        },
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
#    vu.modelRegionBinVar: {
#        'AFRICA': nullBinMethod,
#        'ATLANTIC': nullBinMethod,
#        'AUSTRALIA': nullBinMethod,
#        'CONUS': {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLon,
#                 'bounds': 234.0},
#                {'where': bu.greatBound,
#                 'variable': vu.modVarLon,
#                 'bounds': 294.0},
#                {'where': bu.lessBound,
#                 'variable': vu.modVarLat,
#                 'bounds': 25.0},
#                {'where': bu.greatBound,
#                 'variable': vu.modVarLat,
#                 'bounds': 50.0},
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
#        bu.geoirlatlonboxMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modelVarLon,
#                 'bounds': geoirLonBands['minBounds']},
#                {'where': bu.greatBound,
#                 'variable': vu.modelVarLon,
#                 'bounds': geoirLonBands['maxBounds']},
#                {'where': bu.lessBound,
#                 'variable': vu.modelVarLat,
#                 'bounds': geoirLatBands['minBounds']},
#                {'where': bu.greatBound,
#                 'variable': vu.modelVarLat,
#                 'bounds': geoirLatBands['maxBounds']},
#            ],
#            'values': geoirLonBands['values'],
#        },
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
            AnyBadQC,
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
                AnyBadQC,
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': rangeVar[1],
        }

# Add named cloud fraction-band-specific bins for applicable ranged variables
cldfracBinVars = {
#    vu.obsVarACI: [bu.AsymmetricCloudImpact],
#    vu.obsVarGlint: [bu.GlintAngle],
#    vu.obsVarLandFrac: [vu.landfracGeo],
    vu.obsVarLat: [vu.latMeta],
#    vu.obsVarLT: [bu.LocalHour],
#    vu.obsVarSenZen: [vu.senzenMeta],
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
                AnyBadQC,
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': ['obs','bak','ana','SCI'],
    }

# Add named land fraction-band-specific bins for applicable ranged variables
landfracBinVars = {
    vu.obsVarSCI: [bu.SCIOkamoto],
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
                AnyBadQC,
            ],
            'values': binLims[binVar]['values'],
            'override_exclusiveDiags': ['obs','bak','ana','SCI'],
    }
