#!/usr/bin/env python3

import binning_utils as bu
from collections import defaultdict
from copy import deepcopy
from jediApplicationArgs import jediAppName, nOuterIter
import modelsp_utils as mu
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
modelDiags = [
  'mmgfsan',
  #'rltv_mmgfsan', 'log_mogfsan',
  'sigmaxb',
  #'sigmaxa', 'sigmaxf', 'sigmaxinf',
  #'SRx-eda',
  #'SRx-rtpp',
  'CRxb',
  #'CRxa', 'CRxf', 'CRxinf',
]

## difference diagnostics
absDiagnostics_ = defaultdict(list)
absDiagnostics_['variational'] += ['omb']
absDiagnostics_['hofx'] += ['omf']

## relative difference diagnostics
rltvDiagnostics_ = defaultdict(list)
rltvDiagnostics_['variational'] += ['rltv_omb']
rltvDiagnostics_['hofx'] += ['rltv_omf']

## absolute diagnostics
yDiagnostics_ = defaultdict(list)
yDiagnostics_['variational'] += ['obs', 'bak']
yDiagnostics_['hofx'] += ['obs', 'h(x)']

## non-bias-corrected diagnostics
nobcDiagnostics_ = defaultdict(list)
nobcDiagnostics_['variational'] += ['omb_nobc']

## cloudy radiance diagnostics
cloudyRadDiagnostics = set([
  'SCI-'+bu.OkamotoMethod,
#  'ACI-'+bu.MZ19Method,
#  'MCI',
#  'CFy',
])

## STD diagnostics
absSigmaDiagnostics_ = defaultdict(list)
absSigmaDiagnostics_['variational'] = [
  'omb',
  'sigmaob',
  'sigmab',
  #'ideal-sigmaob',
  'OENIb',
  #'CRyb',
]
absSigmaDiagnostics_['hofx'] = [
  'omf',
  'sigmaof',
  'sigmaf',
  #'ideal-sigmaof',
  #'ideal-rltv_sigmaof',
  #'OENIf',
  'CRyf',
]
spreadDiagnostics = [vu.DiagnosticVars[var] for var in [vu.EddT, vu.HBHT, vu.R, vu.HBHTplusR]]
absSigmaDiagnostics_['hofx'] += spreadDiagnostics

# variational analysis diagnostics
if nOuterIter > 0:
    absDiagnostics_['variational'] += ['oma', 'InnovationRatio']
    nobcDiagnostics_['variational'] += ['oma_nobc']
    rltvDiagnostics_['variational'] += [
      'rltv_oma',
      'rltv_sigmaob',
      'rltv_doadob',
      'rltv_dabdob',
      'ApproxRelativeDoaDob',
    ]
    yDiagnostics_['variational'] += ['ana']
    absSigmaDiagnostics_['variational'] += [
      'oma',
      'sigmaoa',
      'sigmaa',
      #'ideal-sigmaoa',
      'CRya',
      'OENIa',
      'SRy',
      'InnovationRatio',
      'doadob',
      'dabdob',
      'ApproxDoaDob',
    ]

absDiagnostics = set(absDiagnostics_[jediAppName])
rltvDiagnostics = set(rltvDiagnostics_[jediAppName])
yDiagnostics = set(yDiagnostics_[jediAppName])
absSigmaDiagnostics = set(absSigmaDiagnostics_[jediAppName])
nobcDiagnostics = set(nobcDiagnostics_[jediAppName])

## diagnostics for which QC is irrelevant
nonQCedDiagnostics = ['obs']

#==========================
# names and values of bins
#==========================

## heterogeneous/named bins
# note: the order of subplots will follow what is specified here

# latbandsMethod bins, north to south
allNamedLatBands = {}
allNamedLatBands['SXTro'] = (-90.0, -30.0)
allNamedLatBands['STro']  = (-30.0, -5.0)
allNamedLatBands['ITCZ']  = ( -5.0,  5.0)
allNamedLatBands['NTro']  = (  5.0, 30.0)
allNamedLatBands['NXTro'] = ( 30.0,  90.0)

allNamedLatBands['SPol']  = (-90.0, -60.0)
allNamedLatBands['SMid']  = (-60.0, -30.0)
allNamedLatBands['Tro']   = (-30.0,  30.0)
allNamedLatBands['NMid']  = ( 30.0,  60.0)
allNamedLatBands['NPol']  = ( 60.0,  90.0)

allNamedLatBands['SHem']  = (-90.0,  0.0)
allNamedLatBands['NHem']  = (  0.0, 90.0)

namedLatBands = {}
namedLatBands['values'] = ['NXTro','Tro','SXTro']
namedLatBands['starts'] = []
namedLatBands['stops'] = []

for latBand in namedLatBands['values']:
    namedLatBands['starts'].append(allNamedLatBands[latBand][0])
    namedLatBands['stops'].append(allNamedLatBands[latBand][1])


namedTropLatBands = {}
namedTropLatBands['values'] = ['NXTro','NTro','ITCZ','STro','SXTro','Tro']
namedTropLatBands['starts'] = []
namedTropLatBands['stops'] = []

for latBand in namedTropLatBands['values']:
    namedTropLatBands['starts'].append(allNamedLatBands[latBand][0])
    namedTropLatBands['stops'].append(allNamedLatBands[latBand][1])

namedPolarLatBands = {}
namedPolarLatBands['values'] = ['NPol','NMid','Tro','SMid','SPol']
namedPolarLatBands['starts'] = []
namedPolarLatBands['stops'] = []

for latBand in namedPolarLatBands['values']:
    namedPolarLatBands['starts'].append(allNamedLatBands[latBand][0])
    namedPolarLatBands['stops'].append(allNamedLatBands[latBand][1])


# cloudbandsMethod bins
namedCldFracBands = {}
namedCldFracBands['values'] = [bu.clrskyMethod, bu.mixskyMethod, bu.cldskyMethod, bu.allskyMethod]
namedCldFracBands['starts'] = [0.0, bu.clrskyCFThreshold, bu.cldskyCFThreshold, 0.0]
namedCldFracBands['stops'] = [bu.clrskyCFThreshold, bu.cldskyCFThreshold, 1.0, 1.0]


# surfbandsMethod bins
seasurfThresh = 0.05
landsurfThresh = 1.0 - seasurfThresh

namedLandFracBands = {}
namedLandFracBands['values'] = [bu.seasurfMethod, bu.mixsurfMethod, bu.landsurfMethod, bu.allsurfMethod]
namedLandFracBands['starts'] = [0.0, seasurfThresh, landsurfThresh, 0.0]
namedLandFracBands['stops'] = [seasurfThresh, landsurfThresh, 1.0, 1.0]


# geoirlatlonboxMethod bins
geoirLonBands = deepcopy(bu.geoirlatlonBoxParams)
geoirLatBands = deepcopy(bu.geoirlatlonBoxParams)

# store starts/stops
for centerLon in geoirLonBands['centerLon']:
  geoirLonBands['starts'] += \
    [centerLon - bu.geoirMaxZenith]
  geoirLonBands['stops'] += \
    [centerLon + bu.geoirMaxZenith]
  geoirLatBands['starts'] += \
    [-bu.geoirMaxZenith]
  geoirLatBands['stops'] += \
    [bu.geoirMaxZenith]

# ensure positive-definite longitude
for bound in ['starts', 'stops']:
  for ii, lon in enumerate(geoirLonBands[bound]):
    while geoirLonBands[bound][ii] >= 360.:
      geoirLonBands[bound][ii] -= 360.
    while geoirLonBands[bound][ii] < 0.:
      geoirLonBands[bound][ii] += 360.

# namedSatwindObsType bins (see GSI src/read_satwnd.f90)
GSISatwindObsTypes = {
  'GOES-SW': 240,
  'India': 241,
  'JMA-Vis': 242,
  'EUMETSAT-Vis': 243,
  'AVHRR': 244,
  'GOES-IR': 245,
  'GOES-WV-cloud-top': 246,
  'GOES-WV-deep-layer': 247,
  'JMA-WV-deep-layer': 250,
  'GOES-Vis': 251,
  'JMA-IR': 252,
  'EUMETSAT-IR': 253,
  'EUMETSAT-WV-deep-layer': 254,
  'LEOGEO': 255,
  'MODIS-IR': 257,
  'MODIS-WV-cloud-top': 258,
  'MODIS-WV-deep-layer': 259,
  'VIIRS-IR': 260,
}
namedSatWindObsTypes = {}
namedSatWindObsTypes['names'] = list(GSISatwindObsTypes.keys())
namedSatWindObsTypes['values'] = list(GSISatwindObsTypes.values())

## cloud impact coordinates
#notes:
# initial specification: ±np.floor(np.power(np.arange(0.4,3.6,0.33), np.e) * 10.) / 10.
# changed ±0.4 to ±0.5 to force the Count PDF to peak at 0.0 instead of -0.4
diagnosticCloudImpactValuesCoarse = np.array([-4.3, -2.4, -1.1, -0.5, 0., 0.5, 1.1, 2.4, 4.3, 7., 10.5, 15., 20.5, 27.1])

# initial specification: ±np.floor(np.power(np.arange(0.4,2.5,0.2), np.e) * 10.) / 10.
diagnosticCloudImpactValuesFine = np.array([ -4.9, -3.5, -2.4, -1.6, -1., -0.5, -0.2, 0., 0.2, 0.5, 1., 1.6, 2.4, 3.5, 4.9, 6.5, 8.5, 10.8])

## homogeneous bins
# note: the ordering described here does not make any difference
#       figure axes will be monotonically increasing except for pressure

# binLims1D and binLims2D are used to auto-generate a large portion of
# the binVarConfigs dict

############
# -- 1D -- #
############
binLims1D = {}

# base middle values on model-space diagnosticPressures, stop at 25.0 hPa for 30km model top
obsDiagnosticPressures = [25.]+list(mu.diagnosticPressures)+[1000.]

# add more bins as model top goes up
#obsDiagnosticPressures = [5.0, 10.0]+obsDiagnosticPressures

binLims1D[vu.obsVarPrs] = {}
binLims1D[vu.obsVarPrs]['format'] = '{:.1f}'
binLims1D[vu.obsVarPrs]['transforms'] = [np.log, np.exp]
binLims1D[vu.obsVarPrs]['values'] = obsDiagnosticPressures

binLims1D[vu.obsVarAlt] = {}
binLims1D[vu.obsVarAlt]['start']  = 1000.0
binLims1D[vu.obsVarAlt]['stop']   = 30000.0
binLims1D[vu.obsVarAlt]['step']   = 1000.0
binLims1D[vu.obsVarAlt]['format'] = '{:.0f}'

binLims1D[vu.obsVarImpact] = deepcopy(binLims1D[vu.obsVarAlt])

binLims1D[vu.obsVarLat] = {}
binLims1D[vu.obsVarLat]['start']  = -90.0
binLims1D[vu.obsVarLat]['stop']   =  90.0
binLims1D[vu.obsVarLat]['step']   =  10.0
binLims1D[vu.obsVarLat]['format'] = '{:.0f}'

binLims1D[vu.obsVarLT] = {}
binLims1D[vu.obsVarLT]['start']  = bu.LH0
binLims1D[vu.obsVarLT]['stop']   = bu.LH1
binLims1D[vu.obsVarLT]['step']   = bu.LHDT
binLims1D[vu.obsVarLT]['format'] = '{:.0f}'

binLims1D[vu.obsVarSenZen] = {}
binLims1D[vu.obsVarSenZen]['start']  = 0.0
binLims1D[vu.obsVarSenZen]['stop']   = 70.0
binLims1D[vu.obsVarSenZen]['step']   = 5.0
binLims1D[vu.obsVarSenZen]['format'] = '{:.0f}'

binLims1D[vu.obsVarGlint] = {}
binLims1D[vu.obsVarGlint]['start']  = 0.0
binLims1D[vu.obsVarGlint]['stop']   = bu.maxGlint
binLims1D[vu.obsVarGlint]['step']   = 10.0
binLims1D[vu.obsVarGlint]['format'] = '{:.0f}'

binLims1D[vu.obsVarLandFrac] = {}
binLims1D[vu.obsVarLandFrac]['start']  = 0.0
binLims1D[vu.obsVarLandFrac]['stop']   = 1.0
binLims1D[vu.obsVarLandFrac]['step']   = seasurfThresh / 2.0
binLims1D[vu.obsVarLandFrac]['format'] = '{:.3f}'

binLims1D[vu.obsVarCldFracY] = {}
binLims1D[vu.obsVarCldFracY]['start']  = 0.0
binLims1D[vu.obsVarCldFracY]['stop']   = 1.0
binLims1D[vu.obsVarCldFracY]['step']   = bu.clrskyCFThreshold / 2.0
binLims1D[vu.obsVarCldFracY]['format'] = '{:.3f}'
binLims1D[vu.obsVarCldFracX] = deepcopy(binLims1D[vu.obsVarCldFracY])

binLims1D[(vu.obsVarCldFracY, 'coarse')] = {}
binLims1D[(vu.obsVarCldFracY, 'coarse')]['start']  = 0.0
binLims1D[(vu.obsVarCldFracY, 'coarse')]['stop']   = 1.0
binLims1D[(vu.obsVarCldFracY, 'coarse')]['step']   = 0.1
binLims1D[(vu.obsVarCldFracY, 'coarse')]['format'] = '{:.1f}'
binLims1D[(vu.obsVarCldFracX, 'coarse')] = deepcopy(binLims1D[(vu.obsVarCldFracY, 'coarse')])

#binLims1D[(vu.obsVarCldFracY, 'coarse', 'positive')] = {}
#binLims1D[(vu.obsVarCldFracY, 'coarse', 'positive')]['start']  = 0.2
#binLims1D[(vu.obsVarCldFracY, 'coarse', 'positive')]['stop']   = 1.0
#binLims1D[(vu.obsVarCldFracY, 'coarse', 'positive')]['step']   = 0.2
#binLims1D[(vu.obsVarCldFracY, 'coarse', 'positive')]['format'] = '{:.1f}'

binLims1D[vu.obsVarCI] = {}
binLims1D[vu.obsVarCI]['start']  = 0.0
binLims1D[vu.obsVarCI]['stop']   = 60.0
binLims1D[vu.obsVarCI]['step']   = 1.0
binLims1D[vu.obsVarCI]['format'] = '{:.0f}'

binLims1D[(vu.obsVarCI, 'sqrt')] = {}
binLims1D[(vu.obsVarCI, 'sqrt')]['start']  = 0.0
binLims1D[(vu.obsVarCI, 'sqrt')]['stop']   = 7.0
binLims1D[(vu.obsVarCI, 'sqrt')]['nsteps'] = 50
binLims1D[(vu.obsVarCI, 'sqrt')]['format'] = '{:.2f}'

binLims1D[(vu.obsVarCI,'min0.0max0.3')] = {}
binLims1D[(vu.obsVarCI,'min0.0max0.3')]['start']  = 0.0
binLims1D[(vu.obsVarCI,'min0.0max0.3')]['stop']   = 0.3
binLims1D[(vu.obsVarCI,'min0.0max0.3')]['nsteps'] = 15
binLims1D[(vu.obsVarCI,'min0.0max0.3')]['format'] = '{:.3f}'

binLims1D[(vu.obsVarCI,'min1.0max1.3')] = {}
binLims1D[(vu.obsVarCI,'min1.0max1.3')]['start']  = 1.0
binLims1D[(vu.obsVarCI,'min1.0max1.3')]['stop']   = 1.3
binLims1D[(vu.obsVarCI,'min1.0max1.3')]['nsteps'] = 15
binLims1D[(vu.obsVarCI,'min1.0max1.3')]['format'] = '{:.3f}'

binLims1D[vu.obsVarLogCI] = {}
binLims1D[vu.obsVarLogCI]['start']  = 0.6
binLims1D[vu.obsVarLogCI]['stop']   = 60.0
binLims1D[vu.obsVarLogCI]['nsteps'] = 20
binLims1D[vu.obsVarLogCI]['format'] = '{:.2f}'
binLims1D[vu.obsVarLogCI]['transforms'] = [np.log, np.exp]

extendedCIValues = np.append(diagnosticCloudImpactValuesCoarse, [35., 44.1])

binLims1D[(vu.obsVarLogCI, 'coarse')] = {}
binLims1D[(vu.obsVarLogCI, 'coarse')]['format'] = '{:.1f}'
binLims1D[(vu.obsVarLogCI, 'coarse')]['values'] = extendedCIValues

binLims1D[(vu.obsVarLogCI, 'fine')] = {}
binLims1D[(vu.obsVarLogCI, 'fine')]['format'] = '{:.2f}'
#binLims1D[(vu.obsVarLogCI, 'fine')]['values'] = diagnosticCloudImpactValuesFine
binLims1D[(vu.obsVarLogCI, 'fine')]['values'] = np.sqrt(np.abs(extendedCIValues)) * np.sign(extendedCIValues)

binLims1D[vu.obsVarACI] = {}
binLims1D[vu.obsVarACI]['start']  = -20.0
binLims1D[vu.obsVarACI]['stop']   = 20.0
binLims1D[vu.obsVarACI]['step']   = 2.0
binLims1D[vu.obsVarACI]['format'] = '{:.0f}'

binLims1D[vu.obsVarNormDep] = {}
binLims1D[vu.obsVarNormDep]['start']  = -7.0
binLims1D[vu.obsVarNormDep]['stop']   =  7.0
binLims1D[vu.obsVarNormDep]['step']   =  0.25
binLims1D[vu.obsVarNormDep]['format'] = '{:.2f}'

binLims1D[(vu.obsVarDep, bu.allskyBTMethod)] = {}
binLims1D[(vu.obsVarDep, bu.allskyBTMethod)]['start']  = -40.0
binLims1D[(vu.obsVarDep, bu.allskyBTMethod)]['stop']   =  40.0
binLims1D[(vu.obsVarDep, bu.allskyBTMethod)]['step']   =  1.5
binLims1D[(vu.obsVarDep, bu.allskyBTMethod)]['format'] = '{:.1f}'

binLims1D[(vu.obsVarClearSkyDep, bu.allskyBTMethod)] = {}
binLims1D[(vu.obsVarClearSkyDep, bu.allskyBTMethod)]['start']  = -4.0
binLims1D[(vu.obsVarClearSkyDep, bu.allskyBTMethod)]['stop']   =  4.0
binLims1D[(vu.obsVarClearSkyDep, bu.allskyBTMethod)]['nsteps'] =  25
binLims1D[(vu.obsVarClearSkyDep, bu.allskyBTMethod)]['format'] = '{:.2f}'

#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)] = {}
#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)]['start']  = 0.6
#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)]['stop']   =  1.666667
#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)]['nsteps']   =  40.0
#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)]['format'] = '{:.2f}'
#binLims1D[(vu.obsVarDepRatio, bu.allskyBTMethod)]['transforms'] = [np.log, np.exp]

binLims1D[(vu.obsVarLogDepRatio, bu.allskyBTMethod)] = {}
binLims1D[(vu.obsVarLogDepRatio, bu.allskyBTMethod)]['start']  = -0.3
binLims1D[(vu.obsVarLogDepRatio, bu.allskyBTMethod)]['stop']   =  0.3
binLims1D[(vu.obsVarLogDepRatio, bu.allskyBTMethod)]['nsteps'] =  60
binLims1D[(vu.obsVarLogDepRatio, bu.allskyBTMethod)]['format'] = '{:.4f}'

binLims1D[vu.modVarLat] = {}
binLims1D[vu.modVarLat]['start']  = -90.0
binLims1D[vu.modVarLat]['stop']   =  90.0
binLims1D[vu.modVarLat]['step']   =  5.0
binLims1D[vu.modVarLat]['format'] = '{:.0f}'

binLims1D[vu.modVarLev] = {}
binLims1D[vu.modVarLev]['start']  = 1
binLims1D[vu.modVarLev]['stop']   = 55+1
binLims1D[vu.modVarLev]['step']   = 1
binLims1D[vu.modVarLev]['format'] = '{:.0f}'

binLims1D[vu.modVarDiagPrs] = {}
binLims1D[vu.modVarDiagPrs]['format'] = '{:.0f}'
binLims1D[vu.modVarDiagPrs]['values'] = mu.diagnosticPressures

binAxes1D = {}
for binVar, config in binLims1D.items():
  if len(binVar) > 1:
    binAxes1D[binVar] = bu.BinningAxes({(binVar[0], 0): config})
  else:
    binAxes1D[binVar] = bu.BinningAxes({(binVar, 0): config})

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

goodFlags = [0]
goodFlagNames = ['pass']
#Note: uncomment the following to also plot passive flagged data
#goodFlags += [1]
#goodFlagNames += ['passive']


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
#         variable: string or class that is used to initialize the IdBinFunction or BinFunction class;
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


# names for ClearCloudModeBins indicate the data that will be retained
#  I.e., the filters that are present describe the data that will be
#  removed, which is opposite the retained data
#1
OClearMClearBin = [
  {'where': bu.greatBound,
   'variable': vu.cldfracMeta,
   'bounds': 0.0},
  {'where': bu.greatBound,
   'variable': bu.MaximumCFx,
   'bounds': 0.0},
]

#2
OCloudMClearBin = [
  {'where': bu.lessEqualBound,
   'variable': vu.cldfracMeta,
   'bounds': 0.0},
  {'where': bu.greatBound,
   'variable': bu.MaximumCFx,
   'bounds': 0.0},
]

#3
OCloudMCloudBin = [
  {'where': bu.lessEqualBound,
   'variable': vu.cldfracMeta,
   'bounds': 0.0},
  {'where': bu.lessEqualBound,
   'variable': bu.MaximumCFx,
   'bounds': 0.0},
]

#4
OClearMCloudBin = [
  {'where': bu.greatBound,
   'variable': vu.cldfracMeta,
   'bounds': 0.0},
  {'where': bu.lessEqualBound,
   'variable': bu.MaximumCFx,
   'bounds': 0.0},
]

#5
OCloudBin = [
  {'where': bu.lessEqualBound,
   'variable': vu.cldfracMeta,
   'bounds': 0.0},
]

#6
MCloudBin = [
  {'where': bu.lessEqualBound,
   'variable': bu.MaximumCFx,
   'bounds': 0.0},
]

ClearCloudModeBins = {
'1': OClearMClearBin,
'2': OCloudMClearBin,
'3': OCloudMCloudBin,
'4': OClearMCloudBin,
'5': OCloudBin,
'6': MCloudBin,
}

ClearCloudSubgroupCases = [
  #Mode 2
  (vu.cldfracMeta, 'CFy', '2', (vu.obsVarCldFracY, 'coarse'), vu.obsVarCldFracY), # best so far, although STD vs. sqrt(CFy) is parabolic/exponential instead of linear
  (bu.UnscaledOCI, 'UnscaledOCI', '2', (vu.obsVarLogCI, 'coarse'), vu.obsVarCI), # keep for comparison to literature
  #Mode 3
  (bu.SCIOkamoto, 'SCIOkamoto', '3', (vu.obsVarLogCI, 'coarse'), vu.obsVarCI), # keep for comparison to literature, similar problems to using for all modes, flat, large spread, multi-modal
  #Mode 4
  (bu.MaximumCFx, 'CFx', '4', (vu.obsVarCldFracX, 'coarse'), vu.obsVarCldFracX),
  (bu.MCI, 'MCI', '4', (vu.obsVarLogCI, 'coarse'), vu.obsVarCI),
]

binVarConfigs = {
    vu.obsVarQC: {
        bu.passQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': [0],
                 'except_diags': nonQCedDiagnostics}
            ],
            'values': ['pass'],
        },
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
            ],
            'values': badFlagNames,
        },
        bu.allQCMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.selfQCValue,
                 'bounds': goodFlags+badFlags,
                 'except_diags': nonQCedDiagnostics},
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
    vu.obsVarImpact: {
        bu.altjetMethod: {
            'filters': [
# eliminate locations outside bu.alt_jet_min to bu.alt_jet_max
                {'where': bu.lessBound,
                 'variable': vu.impactMeta,
                 'bounds': bu.alt_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.impactMeta,
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
                 'bounds': namedLatBands['starts']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedLatBands['stops']},
                AnyBadQC,
            ],
            'values': namedLatBands['values'],
        },
        bu.troplatbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': namedTropLatBands['starts']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedTropLatBands['stops']},
                AnyBadQC,
            ],
            'values': namedTropLatBands['values'],
        },
        bu.polarlatbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': namedPolarLatBands['starts']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': namedPolarLatBands['stops']},
                AnyBadQC,
            ],
            'values': namedPolarLatBands['values'],
        },
        bu.PjetMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].starts()},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].stops()},
                {'where': bu.lessBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.prsMeta,
                 'bounds': bu.P_jet_max},
                AnyBadQC,
            ],
            'values': binAxes1D[vu.obsVarLat].values(),
        },
        bu.altjetMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].starts()},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].stops()},
                {'where': bu.lessBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.altMeta,
                 'bounds': bu.alt_jet_max},
                AnyBadQC,
            ],
            'values': binAxes1D[vu.obsVarLat].values(),
        },
        bu.impactjetMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].starts()},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': binAxes1D[vu.obsVarLat].stops()},
                {'where': bu.lessBound,
                 'variable': vu.impactMeta,
                 'bounds': bu.alt_jet_min},
                {'where': bu.greatEqualBound,
                 'variable': vu.impactMeta,
                 'bounds': bu.alt_jet_max},
                AnyBadQC,
            ],
            'values': binAxes1D[vu.obsVarLat].values(),
        },
    },
    vu.obsVarLandFrac: {
        bu.surfbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.landfracGeo,
                 'bounds': namedLandFracBands['stops']},
                AnyBadQC,
            ],
            'values': namedLandFracBands['values'],
        },
    },
    vu.obsVarCldFracY: {
        bu.cloudbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.cldfracMeta,
                 'bounds': namedCldFracBands['stops']},
                AnyBadQC,
            ],
            'values': namedCldFracBands['values'],
        },
    },
    vu.obsVarCldFracX: {
        bu.cloudbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.MaximumCFx,
                 'bounds': namedCldFracBands['starts']},
                {'where': bu.greatBound,
                 'variable': bu.MaximumCFx,
                 'bounds': namedCldFracBands['stops']},
                AnyBadQC,
            ],
            'values': namedCldFracBands['values'],
        },
    },
    vu.obsVarCI: {
        bu.OkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.SCIOkamoto,
                 'bounds': binAxes1D[vu.obsVarCI].starts()},
                {'where': bu.greatEqualBound,
                 'variable': bu.SCIOkamoto,
                 'bounds': binAxes1D[vu.obsVarCI].stops()},
            ],
            'values': binAxes1D[vu.obsVarCI].values(),
        },
        bu.CFQuadratureMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.CFQuadrature,
                 'bounds': binAxes1D[vu.obsVarCldFracX].starts()},
                {'where': bu.greatEqualBound,
                 'variable': bu.CFQuadrature,
                 'bounds': binAxes1D[vu.obsVarCldFracX].stops()},
            ],
            'values': binAxes1D[vu.obsVarCldFracX].values(),
        },
    },
    vu.obsVarACI: {
        bu.MZ19Method: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.ACIMZ19,
                 'bounds': binAxes1D[vu.obsVarACI].starts()},
                {'where': bu.greatEqualBound,
                 'variable': bu.ACIMZ19,
                 'bounds': binAxes1D[vu.obsVarACI].stops()},
            ],
            'values': binAxes1D[vu.obsVarACI].values(),
        },
    },
    vu.obsVarNormDep: {
        bu.OkamotoMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': bu.OkamotoNormalizedDeparture,
                 'bounds': binAxes1D[vu.obsVarNormDep].starts()},
                {'where': bu.greatEqualBound,
                 'variable': bu.OkamotoNormalizedDeparture,
                 'bounds': binAxes1D[vu.obsVarNormDep].stops()},
            ],
            'values': binAxes1D[vu.obsVarNormDep].values(),
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
        # TODO: this does not work for SEVIRI due to overlap with 0 deg longitude
        bu.geoirlatlonboxMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.lonMeta,
                 'bounds': geoirLonBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.lonMeta,
                 'bounds': geoirLonBands['stops']},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': geoirLatBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.latMeta,
                 'bounds': geoirLatBands['stops']},
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
#TODO: enable vertical binning by pressure/altitude in MPAS Model space
#    vu.modVarPrs: {
#        bu.identityBinMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarPrs,
#                 'bounds': binAxes1D[vu.modVarPrs].starts()},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarPrs,
#                 'bounds': binAxes1D[vu.modVarPrs].stops()},
#            ],
#            'values': binAxes1D[vu.modVarPrs].values(),
#        },
#    },
#    vu.modVarAlt: {
#        bu.identityBinMethod: {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.modVarAlt,
#                 'bounds': binAxes1D[vu.modVarAlt].starts()},
#                {'where': bu.greatEqualBound,
#                 'variable': vu.modVarAlt,
#                 'bounds': binAxes1D[vu.modVarAlt].stops()},
#            ],
#            'values': binAxes1D[vu.modVarAlt].values(),
#        },
#    },
    vu.modVarLev: {
        bu.identityBinMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.levModel,
                 'bounds': binAxes1D[vu.modVarLev].centrals(astype=int)},
            ],
            'include variables': vu.modVarNamesBase3d,
            'values': binAxes1D[vu.modVarLev].values(),
        },
    },
    vu.modVarDiagPrs: {
        bu.identityBinMethod: {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': vu.prsModelDiag,
                 'bounds': binAxes1D[vu.modVarDiagPrs].centrals},
            ],
            'include variables': vu.modDiagnosticVarNames,
            'values': binAxes1D[vu.modVarDiagPrs].values(),
        },
    },
    vu.modVarLat: {
        bu.identityBinMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latModel,
                 'bounds': binAxes1D[vu.modVarLat].starts()},
                {'where': bu.greatEqualBound,
                 'variable': vu.latModel,
                 'bounds': binAxes1D[vu.modVarLat].stops()},
            ],
            'exclude variables': vu.modDiagnosticVarNames,
            'values': binAxes1D[vu.modVarLat].values(),
        },
        bu.latbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latModel,
                 'bounds': namedLatBands['starts']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latModel,
                 'bounds': namedLatBands['stops']},
            ],
            'exclude variables': vu.modDiagnosticVarNames,
            'values': namedLatBands['values'],
        },
        bu.troplatbandsMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.latModel,
                 'bounds': namedTropLatBands['starts']},
                {'where': bu.greatEqualBound,
                 'variable': vu.latModel,
                 'bounds': namedTropLatBands['stops']},
            ],
            'exclude variables': vu.modDiagnosticVarNames,
            'values': namedTropLatBands['values'],
        },
    },
    vu.modelRegionBinVar: {
#        'AFRICA': nullBinMethod,
#        'ATLANTIC': nullBinMethod,
#        'AUSTRALIA': nullBinMethod,
#        'CONUS': {
#            'filters': [
#                {'where': bu.lessBound,
#                 'variable': vu.lonModel,
#                 'bounds': 234.0},
#                {'where': bu.greatBound,
#                 'variable': vu.lonModel,
#                 'bounds': 294.0},
#                {'where': bu.lessBound,
#                 'variable': vu.latModel,
#                 'bounds': 25.0},
#                {'where': bu.greatBound,
#                 'variable': vu.latModel,
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
        bu.geoirlatlonboxMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.lonModel,
                 'bounds': geoirLonBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.lonModel,
                 'bounds': geoirLonBands['stops']},
                {'where': bu.lessBound,
                 'variable': vu.latModel,
                 'bounds': geoirLatBands['starts']},
                {'where': bu.greatBound,
                 'variable': vu.latModel,
                 'bounds': geoirLatBands['stops']},
            ],
            'exclude variables': vu.modDiagnosticVarNames,
            'values': geoirLonBands['values'],
        },
    },
    vu.noBinVar: {
        bu.noBinMethod: {
            'filters': [
                {'where': bu.lessBound,
                 'variable': vu.levModel,
                 'bounds': [1]},
            ],
            'exclude variables': vu.modDiagnosticVarNames,
            'values': [vu.noBinVar],
        },
    },
}


#=============================
# Parameterized binVarConfigs
#=============================

# model level bins by geoir instrument
for ii, instrument in enumerate(geoirLonBands['values']):
  binVarConfigs[vu.modVarLev][instrument] = {
      'filters': [
          {'where': bu.notEqualBound,
           'variable': vu.levModel,
           'bounds': binAxes1D[vu.modVarLev].centrals(astype=int)},
          {'where': bu.lessBound,
           'variable': vu.lonModel,
           'bounds': geoirLonBands['starts'][ii]},
          {'where': bu.greatBound,
           'variable': vu.lonModel,
           'bounds': geoirLonBands['stops'][ii]},
          {'where': bu.lessBound,
           'variable': vu.latModel,
           'bounds': geoirLatBands['starts'][ii]},
          {'where': bu.greatBound,
           'variable': vu.latModel,
           'bounds': geoirLatBands['stops'][ii]},
      ],
      'include variables': vu.modVarNamesBase3d,
      'values': binAxes1D[vu.modVarLev].values(),
  }

# Add bu.identityBinMethod for identity ranged binning variables
identityRangeBinVars = {
    vu.obsVarAlt: [vu.altMeta, bu.BinMethod.commonDiags],
    vu.obsVarCldFracX: [bu.MaximumCFx, []],
    vu.obsVarCldFracY: [vu.cldfracMeta, []],
    vu.obsVarGlint: [bu.GlintAngle, []],
    vu.obsVarImpact: [vu.impactMeta, bu.BinMethod.commonDiags],
    vu.obsVarLandFrac: [vu.landfracGeo, []],
    vu.obsVarLat: [vu.latMeta, bu.BinMethod.commonDiags],
    vu.obsVarLT: [bu.LocalHour, []],
    vu.obsVarNormDep: [bu.NormalizedDeparture, []],
    vu.obsVarPrs: [vu.prsMeta, bu.BinMethod.commonDiags],
    vu.obsVarSenZen: [vu.senzenMeta, []],
}
for binVar, rangeVar in identityRangeBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    binVarConfigs[binVar][bu.identityBinMethod] = {
        'filters': [
            {'where': bu.lessBound,
             'variable': rangeVar[0],
             'bounds': binAxes1D[binVar].starts()},
            {'where': bu.greatEqualBound,
             'variable': rangeVar[0],
             'bounds': binAxes1D[binVar].stops()},
            AnyBadQC,
        ],
        'values': binAxes1D[binVar].values(),
        'override_exclusiveDiags': rangeVar[1],
    }


# Add Departure bins for subgroups of Cloud Impact, ClearCloudModeBins, CldFracY, etc...
for departureVar, departureFunction in [
  (vu.obsVarDep, bu.Departure),
  ]:

  binVarConfigs[departureVar] = {}
  binAxesDeparture = binAxes1D[(departureVar, bu.allskyBTMethod)]

  # Cloud Impact (CI) Bins
  # Add SCI-valued Departure bins for OkamotoMethod, QuadratureMethod, OCIMethod, and MCIMethod
  for (CIMethod, CIVariable, CILabel, axisKey) in list(zip(
    [bu.OkamotoMethod,
    ],
    [bu.SCIOkamoto,
    ],
    ['SCI',
    ],
    [(vu.obsVarLogCI, 'coarse'),
    ],
    )):
    for ii, CIBinVal in enumerate(binAxes1D[axisKey].values()):
      binVarConfigs[departureVar][CIMethod+','+CILabel+'='+CIBinVal] = {
        'filters': [
          {'where': bu.lessBound,
           'variable': departureFunction,
           'bounds': binAxesDeparture.starts()},
          {'where': bu.greatEqualBound,
           'variable': departureFunction,
           'bounds': binAxesDeparture.stops()},
          {'where': bu.lessBound,
           'variable': CIVariable,
           'bounds': binAxes1D[axisKey].starts()[ii]},
          {'where': bu.greatEqualBound,
           'variable': CIVariable,
           'bounds': binAxes1D[axisKey].stops()[ii]},
          AnyBadQC,
        ],
        'values': binAxesDeparture.values(),
      }

  # ClearCloudModeBins
  for mode, ClearCloudModeBin in ClearCloudModeBins.items():
    binVarConfigs[departureVar][bu.ClearCloudModeMethod+'='+mode] = {
      'filters': [
        {'where': bu.lessBound,
         'variable': departureFunction,
         'bounds': binAxesDeparture.starts()},
        {'where': bu.greatEqualBound,
         'variable': departureFunction,
         'bounds': binAxesDeparture.stops()},
        AnyBadQC,
      ]+ClearCloudModeBin,
      'values': binAxesDeparture.values(),
    }

  for (subgroupVariable, subgroupLabel, mode, pdfAxisKey, fittingAxisKey) in ClearCloudSubgroupCases:

    binAxes = binAxes1D[pdfAxisKey]

    for subgroupBinVal, subgroupStart, subgroupStop in zip(
      binAxes.values(), binAxes.starts(), binAxes.stops()):
      binVarConfigs[departureVar][bu.ClearCloudModeMethod+','+subgroupLabel+'='+mode+','+subgroupBinVal] = {
        'filters': [
          {'where': bu.lessBound,
           'variable': departureFunction,
           'bounds': binAxesDeparture.starts()},
          {'where': bu.greatEqualBound,
           'variable': departureFunction,
           'bounds': binAxesDeparture.stops()},
          {'where': bu.lessBound,
           'variable': subgroupVariable,
           'bounds': subgroupStart},
          {'where': bu.greatEqualBound,
           'variable': subgroupVariable,
           'bounds': subgroupStop},
          AnyBadQC,
        ]+ClearCloudModeBins[mode],
        'values': binAxesDeparture.values(),
      }

#Normalized departures binned into ClearCloudModeBins
departureVar = vu.obsVarNormDep
for departureFunction, obserrorMethod in [
  (bu.NormalizedDeparture, ''),
  (bu.OkamotoNormalizedDeparture, bu.OkamotoMethod+','),
  ]:

  binAxesDeparture = binAxes1D[departureVar]

  # ClearCloudModeBins
  for mode, ClearCloudModeBin in ClearCloudModeBins.items():
    binVarConfigs[departureVar][obserrorMethod+bu.ClearCloudModeMethod+'='+mode] = {
      'filters': [
        {'where': bu.lessBound,
         'variable': departureFunction,
         'bounds': binAxesDeparture.starts()},
        {'where': bu.greatEqualBound,
         'variable': departureFunction,
         'bounds': binAxesDeparture.stops()},
        AnyBadQC,
      ]+ClearCloudModeBin,
      'values': binAxesDeparture.values(),
    }


# Add named latitude-band-specific bins for observation-space ranged variables
obsLatBinVars = {
    vu.obsVarAlt: [vu.altMeta, [], namedPolarLatBands],
    vu.obsVarImpact: [vu.impactMeta, [], namedPolarLatBands],
    vu.obsVarPrs: [vu.prsMeta, [], namedLatBands],
}
for binVar, rangeVar in obsLatBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    latBands = deepcopy(rangeVar[2])
    for latBand, latStart, latStop in zip(
        latBands['values'],
        latBands['starts'],
        latBands['stops'],
        ):
        binVarConfigs[binVar][latBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].starts()},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].stops()},
                {'where': bu.lessBound,
                 'variable': vu.latMeta,
                 'bounds': latStart},
                {'where': bu.greatEqualBound,
                 'variable': vu.latMeta,
                 'bounds': latStop},
                AnyBadQC,
            ],
            'values': binAxes1D[binVar].values(),
            'override_exclusiveDiags': rangeVar[1],
        }

# Add named latitude-band-specific bins for model-space ranged variables
modelLatBinVars = {
    vu.modVarLev: [vu.levModel, int, namedTropLatBands, vu.modVarNamesBase3d],
    vu.modVarDiagPrs: [vu.prsModelDiag, float, namedLatBands, vu.modDiagnosticVarNames],
}
for binVar, rangeVar in modelLatBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    latBands = deepcopy(rangeVar[2])
    for latBand, latStart, latStop in zip(
        latBands['values'],
        latBands['starts'],
        latBands['stops'],
        ):
        binVarConfigs[binVar][latBand] = {
            'filters': [
                {'where': bu.notEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].centrals(astype=rangeVar[1])},
                {'where': bu.lessBound,
                 'variable': vu.latModel,
                 'bounds': latStart},
                {'where': bu.greatEqualBound,
                 'variable': vu.latModel,
                 'bounds': latStop},
            ],
            'values': binAxes1D[binVar].values(),
            'include variables': rangeVar[3]
        }

# pressure level bins by satellite ObsType and namedLatBand
for latBand, latStart, latStop in zip(
  namedLatBands['values'],
  namedLatBands['starts'],
  namedLatBands['stops'],
  ):

  for name, ID in zip(
    namedSatWindObsTypes['names'],
    namedSatWindObsTypes['values']):

    # add basic QC filter and composite binVal's (values)
    #binVarConfigs[vu.obsVarPrs][name] = {
    binVarConfigs[vu.obsVarPrs][latBand+'-ObsType'+str(ID)] = {
      'filters': [
        {'where': bu.notEqualBound, # eliminates locations where vu.selfObsType != ID
         'variable': vu.selfObsType,
         'bounds': ID},
        {'where': bu.lessBound,
         'variable': vu.prsMeta,
         'bounds': binAxes1D[vu.obsVarPrs].starts()},
        {'where': bu.greatEqualBound,
         'variable': vu.prsMeta,
         'bounds': binAxes1D[vu.obsVarPrs].stops()},
        {'where': bu.lessBound,
         'variable': vu.latMeta,
         'bounds': latStart},
        {'where': bu.greatEqualBound,
         'variable': vu.latMeta,
         'bounds': latStop},
        AnyBadQC,
      ],
      'values': binAxes1D[vu.obsVarPrs].values(),
    }


# Add named cloud fraction-band-specific bins for applicable ranged variables
cldfracBinVars = {
#    vu.obsVarACI: [bu.ACIMZ19],
#    vu.obsVarGlint: [bu.GlintAngle],
#    vu.obsVarLandFrac: [vu.landfracGeo],
    vu.obsVarLat: [vu.latMeta],
#    vu.obsVarLT: [bu.LocalHour],
#    vu.obsVarSenZen: [vu.senzenMeta],
}
for binVar, rangeVar in cldfracBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for cldBand, cldStart, cldStop in zip(
        namedCldFracBands['values'],
        namedCldFracBands['starts'],
        namedCldFracBands['stops'],
        ):
        binVarConfigs[binVar][cldBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].starts()},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].stops()},
                {'where': bu.lessBound,
                 'variable': vu.cldfracMeta,
                 'bounds': cldStart},
                {'where': bu.greatBound,
                 'variable': vu.cldfracMeta,
                 'bounds':cldStop},
                AnyBadQC,
            ],
            'values': binAxes1D[binVar].values(),
            #'override_exclusiveDiags': bu.BinMethod.commonDiags,
    }

# Add named land fraction-band-specific bins for applicable ranged variables
landfracBinVars = {
    vu.obsVarCI: [bu.SCIOkamoto],
}
for binVar, rangeVar in landfracBinVars.items():
    if binVar not in binVarConfigs: binVarConfigs[binVar] = {}
    for landBand, landStart, landStop in zip(
        namedLandFracBands['values'],
        namedLandFracBands['starts'],
        namedLandFracBands['stops'],
        ):
        binVarConfigs[binVar][landBand] = {
            'filters': [
                {'where': bu.lessBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].starts()},
                {'where': bu.greatEqualBound,
                 'variable': rangeVar[0],
                 'bounds': binAxes1D[binVar].stops()},
                {'where': bu.lessBound,
                 'variable': vu.landfracGeo,
                 'bounds':landStart},
                {'where': bu.greatBound,
                 'variable': vu.landfracGeo,
                 'bounds': landStop},
                AnyBadQC,
            ],
            'values': binAxes1D[binVar].values(),
            #'override_exclusiveDiags': bu.BinMethod.commonDiags,
    }


############
# -- 2D -- #
############
binLims2D = {}
binVars2D = {}
whereVars2D = {}
includeQC2D = {}
includeVars2D = {}
binAxes2D = {}
binMethods2D = defaultdict(list)

## vertical coordinates
binLims2D[(vu.obsVarPrs, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarPrs, bu.noBinMethod)]['format'] = '{:.1f}'
binLims2D[(vu.obsVarPrs, bu.noBinMethod)]['transforms'] = [np.log, np.exp]
binLims2D[(vu.obsVarPrs, bu.noBinMethod)]['values'] = obsDiagnosticPressures

binLims2D[(vu.obsVarAlt, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarAlt, bu.noBinMethod)]['start']  = 1000.0
binLims2D[(vu.obsVarAlt, bu.noBinMethod)]['stop']   = 30000.0
binLims2D[(vu.obsVarAlt, bu.noBinMethod)]['step']   = 2000.0
binLims2D[(vu.obsVarAlt, bu.noBinMethod)]['format'] = '{:.0f}'

binLims2D[(vu.obsVarImpact, bu.noBinMethod)] = deepcopy(binLims2D[(vu.obsVarAlt, bu.noBinMethod)])

binLims2D[(vu.modVarLev, bu.noBinMethod)] = {}
binLims2D[(vu.modVarLev, bu.noBinMethod)]['start']  = 1
binLims2D[(vu.modVarLev, bu.noBinMethod)]['stop']   = 56
binLims2D[(vu.modVarLev, bu.noBinMethod)]['nsteps'] = 10
binLims2D[(vu.modVarLev, bu.noBinMethod)]['format'] = '{:.0f}'

binLims2D[(vu.modVarDiagPrs, bu.noBinMethod)] = {}
binLims2D[(vu.modVarDiagPrs, bu.noBinMethod)]['format'] = '{:.0f}'
binLims2D[(vu.modVarDiagPrs, bu.noBinMethod)]['values'] = mu.diagnosticPressures

## horizontal coordinates
# globally distributed instruments/model
binLims2D[(vu.obsVarLat, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarLat, bu.noBinMethod)]['start']  = -88.0
binLims2D[(vu.obsVarLat, bu.noBinMethod)]['stop']   =  88.0
binLims2D[(vu.obsVarLat, bu.noBinMethod)]['step']   =  11.0
binLims2D[(vu.obsVarLat, bu.noBinMethod)]['format'] = '{:.1f}'

binLims2D[(vu.obsVarLon, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarLon, bu.noBinMethod)]['start']  = 0.0
binLims2D[(vu.obsVarLon, bu.noBinMethod)]['stop']   = 360.0
binLims2D[(vu.obsVarLon, bu.noBinMethod)]['step']   = 10.0
binLims2D[(vu.obsVarLon, bu.noBinMethod)]['period'] = 360.0
binLims2D[(vu.obsVarLon, bu.noBinMethod)]['format'] = '{:.1f}'

binLims2D[(vu.modVarLat, bu.noBinMethod)] = deepcopy(binLims2D[(vu.obsVarLat, bu.noBinMethod)])
binLims2D[(vu.modVarLon, bu.noBinMethod)] = deepcopy(binLims2D[(vu.obsVarLon, bu.noBinMethod)])

# geoir instruments
for ii, instrument in enumerate(geoirLonBands['values']):
  binLims2D[(vu.obsVarLat, instrument)] = {}
  binLims2D[(vu.obsVarLat, instrument)]['start']  = geoirLatBands['starts'][ii]
  binLims2D[(vu.obsVarLat, instrument)]['stop']   = geoirLatBands['stops'][ii]
  binLims2D[(vu.obsVarLat, instrument)]['step']   = 8.125
  binLims2D[(vu.obsVarLat, instrument)]['format'] = '{:.1f}'

  binLims2D[(vu.obsVarLon, instrument)] = {}
  binLims2D[(vu.obsVarLon, instrument)]['start']  = geoirLonBands['starts'][ii]
  binLims2D[(vu.obsVarLon, instrument)]['stop']   = geoirLonBands['stops'][ii]
  binLims2D[(vu.obsVarLon, instrument)]['step']   = 8.125
  binLims2D[(vu.obsVarLon, instrument)]['period'] = 360.0
  binLims2D[(vu.obsVarLon, instrument)]['format'] = '{:.1f}'

binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)]['start']  = 0.0
binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)]['stop']   = 1.0
binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)]['step']   = 0.05
binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)]['format'] = '{:.2f}'
binLims2D[(vu.obsVarCldFracX, bu.noBinMethod)] = deepcopy(binLims2D[(vu.obsVarCldFracY, bu.noBinMethod)])

binLims2D[(vu.obsVarMCI, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarMCI, bu.noBinMethod)]['values']   = diagnosticCloudImpactValuesCoarse
binLims2D[(vu.obsVarMCI, bu.noBinMethod)]['format'] = '{:.2f}'

binLims2D[(vu.obsVarOCI, bu.noBinMethod)] = {}
binLims2D[(vu.obsVarOCI, bu.noBinMethod)]['values']   = diagnosticCloudImpactValuesCoarse
binLims2D[(vu.obsVarOCI, bu.noBinMethod)]['format'] = '{:.2f}'

binLims2D[(vu.obsVarY, bu.allskyBTMethod)] = {}
binLims2D[(vu.obsVarY, bu.allskyBTMethod)]['start']  = 190.0
binLims2D[(vu.obsVarY, bu.allskyBTMethod)]['stop']   = 290.0
binLims2D[(vu.obsVarY, bu.allskyBTMethod)]['nsteps'] = 33
binLims2D[(vu.obsVarY, bu.allskyBTMethod)]['format'] = '{:.2f}'
#binLims2D[(vu.obsVarY, bu.allskyBTMethod)]['transforms'] = [np.log, np.exp]
binLims2D[(vu.obsVarH, bu.allskyBTMethod)] = deepcopy(binLims2D[(vu.obsVarY, bu.allskyBTMethod)])

## particular axis pairs
LonLat2D = 'LonLat2D'
binVars2D[LonLat2D] = [vu.obsVarLon, vu.obsVarLat]
whereVars2D[LonLat2D] = [vu.lonMeta, vu.latMeta]
binMethods2D[LonLat2D] += [bu.noBinMethod]
binMethods2D[LonLat2D] += geoirLonBands['values']

ModelLonLat2D = 'ModelLonLat2D'
binVars2D[ModelLonLat2D] = [vu.modVarLon, vu.modVarLat]
whereVars2D[ModelLonLat2D] = [vu.lonModel, vu.latModel]
includeQC2D[ModelLonLat2D] = False
includeVars2D[ModelLonLat2D] = vu.modVarNames2d+vu.modVarNamesBase3d
binMethods2D[ModelLonLat2D] += [bu.noBinMethod]

LatPrs2D = 'LatPrs2D'
binVars2D[LatPrs2D] = [vu.obsVarLat, vu.obsVarPrs]
whereVars2D[LatPrs2D] = [vu.latMeta, vu.prsMeta]
binMethods2D[LatPrs2D] += [bu.noBinMethod]

LatAlt2D = 'LatAlt2D'
binVars2D[LatAlt2D] = [vu.obsVarLat, vu.obsVarAlt]
whereVars2D[LatAlt2D] = [vu.latMeta, vu.altMeta]
binMethods2D[LatAlt2D] += [bu.noBinMethod]

LatImpact2D = 'LatImpact2D'
binVars2D[LatImpact2D] = [vu.obsVarLat, vu.obsVarImpact]
whereVars2D[LatImpact2D] = [vu.latMeta, vu.impactMeta]
binMethods2D[LatImpact2D] += [bu.noBinMethod]

ModelLatLev2D = 'ModelLatLev2D'
binVars2D[ModelLatLev2D] = [vu.modVarLat, vu.modVarLev]
whereVars2D[ModelLatLev2D] = [vu.latModel, vu.levModel]
includeQC2D[ModelLatLev2D] = False
includeVars2D[ModelLatLev2D] = vu.modVarNamesBase3d
binMethods2D[ModelLatLev2D] += [bu.noBinMethod]

ModelLatPrs2D = 'ModelLatPrs2D'
binVars2D[ModelLatPrs2D] = [vu.modVarLat, vu.modVarDiagPrs]
whereVars2D[ModelLatPrs2D] = [vu.latModel, vu.prsModelDiag]
includeQC2D[ModelLatPrs2D] = False
includeVars2D[ModelLatPrs2D] = vu.modDiagnosticVarNames
binMethods2D[ModelLatPrs2D] += [bu.noBinMethod]

ObsModel2D = 'ObsModel2D'
binVars2D[ObsModel2D] = [vu.obsVarY, vu.obsVarH]
whereVars2D[ObsModel2D] = [vu.selfObsValue, vu.selfHofXValue]
binMethods2D[ObsModel2D] += [bu.allskyBTMethod]

CldFrac2D = 'CldFrac2D'
binVars2D[CldFrac2D] = [vu.obsVarCldFracY, vu.obsVarCldFracX]
whereVars2D[CldFrac2D] = [vu.cldfracMeta, bu.MaximumCFx]
binMethods2D[CldFrac2D] += [bu.noBinMethod]

#CloudImpact2D = 'CloudImpact2D'
#binVars2D[CloudImpact2D] = [vu.obsVarMCI, vu.obsVarOCI]
#whereVars2D[CloudImpact2D] = [bu.MCI, bu.OCI]
#binMethods2D[CloudImpact2D] += [bu.noBinMethod]

#OkamotoCloudImpact2D = 'OkamotoCloudImpact2D'
#binVars2D[OkamotoCloudImpact2D] = [vu.obsVarMCI, vu.obsVarOCI]
#whereVars2D[OkamotoCloudImpact2D] = [bu.MCI, bu.UnscaledOCI]
#binMethods2D[OkamotoCloudImpact2D] += [bu.noBinMethod]

## construct all bu.BinningAxes objects
for binVar2D, binVars in binVars2D.items():
    for binMethod in binMethods2D[binVar2D]:
        axes = {}
        for ii, binVar in enumerate(binVars):
            axes[(binVar, ii)] = binLims2D[(binVar, binMethod)]
        binAxes2D[(binVar2D, binMethod)] = bu.BinningAxes(axes)


# parse all available binAxes2D objects
for (binVar2D, binMethod), binAxes in binAxes2D.items():

    if binVar2D not in binVarConfigs: binVarConfigs[binVar2D] = {}

    # add basic QC filter and composite binVal's (values)
    binVarConfigs[binVar2D][binMethod] = {
        'filters': [],
        'values': binAxes.values(),
        'override_exclusiveDiags': bu.BinMethod.commonDiags,
    }
    if includeQC2D.get(binVar2D, True):
      binVarConfigs[binVar2D][binMethod]['filters'].append(AnyBadQC)

    includeVars = includeVars2D.get(binVar2D, [])
    if len(includeVars) > 0:
      binVarConfigs[binVar2D][binMethod]['include variables'] = includeVars

    # add bounding filters for each BinningAxis
    for binVar, whereVar in zip(binVars2D[binVar2D], whereVars2D[binVar2D]):
        binVarConfigs[binVar2D][binMethod]['filters'].append(
            {'where': bu.lessBound,
             'variable': whereVar,
             'bounds': binAxes.starts(binVar)}
        )
        binVarConfigs[binVar2D][binMethod]['filters'].append(
            {'where': bu.greatEqualBound,
             'variable': whereVar,
             'bounds': binAxes.stops(binVar)}
        )

## ClearCloudMode-specific 2D binning, including for some bu.BinMethod.exclusiveDiags
#  in addition to all standard diagnostics (e.g., omf, sigmaof, sigmaf, etc...)

if False:

  for binVar2D, binMethod0 in [
    (ObsModel2D, bu.allskyBTMethod),
  #  (CldFrac2D, bu.noBinMethod),
  ]:
    binVars = binVars2D[binVar2D]
    binAxes = binAxes2D[(binVar2D, binMethod0)]

    for mode, ClearCloudModeBin in ClearCloudModeBins.items():

      binMethod = ''
      if binMethod0 != bu.noBinMethod:
        binMethod += binMethod0+'_'
      binMethod += bu.ClearCloudModeMethod+'='+mode

      # add basic QC filter and composite binVal's (values)
      binVarConfigs[binVar2D][binMethod] = {
        'filters': [
          AnyBadQC,
        ]+ClearCloudModeBin,
        'values': binAxes.values(),
      }

      # add bounding filters for each BinningAxis
      for binVar, whereVar in zip(binVars2D[binVar2D], whereVars2D[binVar2D]):
        binVarConfigs[binVar2D][binMethod]['filters'].append(
          {'where': bu.lessBound,
           'variable': whereVar,
           'bounds': binAxes.starts(binVar)}
        )
        binVarConfigs[binVar2D][binMethod]['filters'].append(
          {'where': bu.greatEqualBound,
           'variable': whereVar,
           'bounds': binAxes.stops(binVar)}
        )

if True:
  # latitude-band-specific 2D binning
  # #Poly2DLat: the CldFrac2D binning below is used to generate 5-band latitudinally distributed
  # polynomial fitting coefficients for STD(OMF) vs. (CFy, CFx)
  for binVar2D, binMethod0 in [
    #(ObsModel2D, bu.allskyBTMethod),
    (CldFrac2D, bu.noBinMethod),
  ]:
    binVars = binVars2D[binVar2D]
    binAxes = binAxes2D[(binVar2D, binMethod0)]

    for latBand, latStart, latStop in zip(
      namedTropLatBands['values'],
      namedTropLatBands['starts'],
      namedTropLatBands['stops']):

      binMethod = ''
      if binMethod0 != bu.noBinMethod:
        binMethod += binMethod0+'_'
      binMethod += latBand

      # add basic QC filter and composite binVal's (values)
      binVarConfigs[binVar2D][binMethod] = {
        'filters': [
          {'where': bu.lessBound,
           'variable': vu.latMeta,
           'bounds': latStart},
          {'where': bu.greatEqualBound,
           'variable': vu.latMeta,
           'bounds': latStop},
          AnyBadQC,
        ],
        'values': binAxes.values(),
      }

      # add bounding filters for each BinningAxis
      for binVar, whereVar in zip(binVars2D[binVar2D], whereVars2D[binVar2D]):
        binVarConfigs[binVar2D][binMethod]['filters'].append(
          {'where': bu.lessBound,
           'variable': whereVar,
           'bounds': binAxes.starts(binVar)}
        )
        binVarConfigs[binVar2D][binMethod]['filters'].append(
          {'where': bu.greatEqualBound,
           'variable': whereVar,
           'bounds': binAxes.stops(binVar)}
        )
