#!/usr/bin/env python3

import binning_utils as bu
import predefined_configs as pconf
from collections import defaultdict
from copy import deepcopy
import os
import var_utils as vu

#==============================================================
# Sub-selections of pconf.binVarConfigs for specific DiagSpaces
#==============================================================

#################################################################
## Generic binVarConfigs that apply to all observation categories
#################################################################
# priority 1
obsBinVars = defaultdict(list)
obsBinVars[vu.obsVarQC] += [bu.goodQCMethod, bu.badQCMethod, bu.allQCMethod]
obsBinVars[vu.obsVarLat] += [bu.identityBinMethod]
# TODO: re-enable obsVarLT binning after adapting to new
#       IODA dateTime format
#obsBinVars[vu.obsVarLT] += [bu.identityBinMethod]
obsBinVars[vu.obsVarNormDep] += [bu.identityBinMethod]
obsBinVars[vu.obsRegionBinVar] += ['CONUS']
obsBinVars[vu.obsRegionBinVar] += [bu.geoirlatlonboxMethod]


################################
## binVarConfigs for surface obs
################################
# priority 1
surfBinVars = deepcopy(obsBinVars)
surfBinVars[pconf.LonLat2D] += [bu.noBinMethod]
surfBinVars[vu.obsVarLat] += [bu.latbandsMethod]

##########################################################
## binVarConfigs for profile obs w/ pressure vertical bins
##########################################################
# priority 1
profilePressureBinVars = deepcopy(obsBinVars)
profilePressureBinVars[vu.obsVarPrs] += [bu.identityBinMethod, bu.PjetMethod]
profilePressureBinVars[vu.obsVarLat] += [bu.latbandsMethod]
# priority 2
profilePressureBinVars[vu.obsVarLat] += [bu.PjetMethod]

# pseudo-2D pressure bins with named latitude-band methods
# priority 1
for latBand in pconf.namedLatBands['values']:
    profilePressureBinVars[vu.obsVarPrs] += [latBand]

# 2D latitude-pressure bins
# priority 1b
profilePressureBinVars[pconf.LatPrs2D] += [bu.noBinMethod]

# 2D bins for satwind
# priority 1b
satwindBinVars = deepcopy(profilePressureBinVars)
satwindBinVars[pconf.LonLat2D] += [bu.noBinMethod]


##########################################################
## binVarConfigs for profile obs w/ altitude vertical bins
##########################################################
# priority 1
gnssroBinVars = deepcopy(obsBinVars)
gnssroBinVars[vu.obsVarAlt] += [bu.identityBinMethod, bu.altjetMethod]
gnssroBinVars[vu.obsVarLat] += [bu.latbandsMethod]
gnssroBinVars[vu.obsVarLat] += [bu.altjetMethod]

# pseudo-2D altitude bins with named latitude-band methods
# priority 1
for latBand in pconf.namedLatBands['values']:
    gnssroBinVars[vu.obsVarAlt] += [latBand]

# 2D bins
# priority 1b
gnssroBinVars[pconf.LonLat2D] += [bu.noBinMethod]
gnssroBinVars[pconf.LatAlt2D] += [bu.noBinMethod]


#################################
## binVarConfigs for radiance obs
#################################
# priority 1
radianceBinVars = deepcopy(obsBinVars)
radianceBinVars[vu.obsVarLat] += [bu.troplatbandsMethod]
# priority 2
radianceBinVars[vu.obsVarGlint] += [bu.identityBinMethod]
radianceBinVars[vu.obsVarSenZen] += [bu.identityBinMethod]


# binBYLandFrac controls whether to bin by land-ocean fractions and categories
# NOTE: requires geovals to be written for applicable DiagSpaces
# priority 4
binBYLandFrac = False
if binBYLandFrac:
    radianceBinVars[vu.obsVarLandFrac] += [bu.identityBinMethod,
                                           bu.surfbandsMethod]


############################################
## binVarConfigs for polar-orbitting obs,
## including cloud metrics
## e.g., AMSUA
############################################
polarBinVars = deepcopy(radianceBinVars)

# 2D bins
# priority 1b
polarBinVars[pconf.LonLat2D] += [bu.noBinMethod]

# cloud-specific metrics
polarcldBinVars = deepcopy(polarBinVars)

# symmetric and asymmetric cloud impact
# priority 4
#polarcldBinVars[vu.obsVarACI] += [bu.OkamotoMethod]
#polarcldBinVars[vu.obsVarCI] += [bu.OkamotoMethod]



##########################################
## binVarConfigs for geostationary IR obs
## including cloud metrics
## e.f., GOES-ABI, Himawari-AHI
##########################################
# cloud fraction
# priority 1
geoirBinVars = deepcopy(radianceBinVars)
geoirBinVars[vu.obsVarCldFracY] += [
  bu.identityBinMethod,
  bu.cloudbandsMethod,
]
geoirBinVars[vu.obsVarCldFracX] += [
  bu.identityBinMethod,
  bu.cloudbandsMethod,
]

# Binning variables with clr-/cld-sky methods (not typically needed)
#for var in pconf.cldfracBinVars.keys():
#    if var != vu.obsVarLandFrac or binBYLandFrac:
#        geoirBinVars[var] += [bu.clrskyMethod]
#        geoirBinVars[var] += [bu.cldskyMethod]

# Binning variables with land category methods
# priority 4
if binBYLandFrac:
    for var in pconf.landfracBinVars.keys():
        for surfBand in pconf.namedLandFracBands['values']:
            geoirBinVars[var] += [surfBand]

# symmetric cloud impact
# priority 2
geoirCIMethods = [
    bu.OkamotoMethod,
#    bu.CFQuadratureMethod,
]
for method in geoirCIMethods:
    geoirBinVars[vu.obsVarCI] += [method]
# uncomment to diagnose ENI = (HofX - ObsValue) / ObsError for SCI-parameterized ObsError
# note: requires SCIErrParams to be defined for all all DiagSpaces that use it
    geoirBinVars[vu.obsVarNormDep] += [method]

for mode, ClearCloudModeBin in pconf.ClearCloudModeBins.items():
  geoirBinVars[vu.obsVarNormDep] += [bu.ClearCloudModeMethod+'='+mode]
  geoirBinVars[vu.obsVarNormDep] += [bu.OkamotoMethod+','+bu.ClearCloudModeMethod+'='+mode]

# asymmetric cloud impact
# priority 3
geoirACIMethods = [
    bu.MZ19Method,
#    bu.QuadratureMethod,
]
for method in geoirACIMethods:
    geoirBinVars[vu.obsVarACI] += [method]

# Cloud Impact Departure bins (expensive, can be disabled unless investigating departure PDF's)
# priority 3
for method in pconf.binVarConfigs[vu.obsVarDep].keys():
  geoirBinVars[vu.obsVarDep] += [method]
#for method in pconf.binVarConfigs[vu.obsVarLogDepRatio].keys():
#  geoirBinVars[vu.obsVarLogDepRatio] += [method]

# 2D bins
# model vs. obs PDF
# priority 3
geoirBinVars[pconf.ObsModel2D] += [bu.allskyBTMethod]

# cloud-fraction
# priority 2
geoirBinVars[pconf.CldFrac2D] += [bu.noBinMethod]
# uncomment to enable Poly2DLat fitting
# priority 2b
#for latBand in pconf.namedTropLatBands['values']:
#  geoirBinVars[pconf.CldFrac2D] += [latBand]

# used for diagnosing 2D statistical patterns vs. model-CI and obs-CI
# priority 4
#geoirBinVars[pconf.CloudImpact2D] += [bu.noBinMethod]
#geoirBinVars[pconf.OkamotoCloudImpact2D] += [bu.noBinMethod]

# instrument-specific 2D bins
# priority 1b
abi_g16_binVars = deepcopy(geoirBinVars)
abi_g16_binVars[pconf.LonLat2D] += [bu.abi_g16]

ahi_himawari8_binVars = deepcopy(geoirBinVars)
ahi_himawari8_binVars[pconf.LonLat2D] += [bu.ahi_himawari8]

#seviri_m08_binVars = deepcopy(radianceBinVars)
#seviri_m08_binVars[pconf.LonLat2D] += [bu.seviri_m08]

#########################################
# binVarConfigs for model space variables
#########################################
# priority 1
modelBinVars = defaultdict(list)
modelBinVars[vu.noBinVar] += [bu.noBinMethod]
modelBinVars[vu.modVarLat] += [bu.identityBinMethod]
modelBinVars[vu.modVarLat] += [bu.troplatbandsMethod]
modelBinVars[vu.modVarLev] += [bu.identityBinMethod]
modelBinVars[vu.modelRegionBinVar] += [bu.geoirlatlonboxMethod]

# pseudo-2D model level bins with named latitude-band methods
# priority 1
for latBand in pconf.namedTropLatBands['values']:
  modelBinVars[vu.modVarLev] += [latBand]

# pseudo-2D model level bins with GEOIR lat-lon boxes
# priority 2
for instrument in bu.geoirlatlonBoxParams['values']:
  modelBinVars[vu.modVarLev] += [instrument]

# 2D bins
# priority 2
modelBinVars[pconf.ModelLonLat2D] += [bu.noBinMethod]
modelBinVars[pconf.ModelLatLev2D] += [bu.noBinMethod]

#=======================
# DiagSpace definitions
# e.g. IODA ObsSpace
#      MPAS ModelSpace
#=======================
nullBinVars = {vu.miss_s: []}

nullDiagSpaceInfo = {
    'DiagSpaceGrp': vu.miss_s,
    'process': False,
    'binVarConfigs': nullBinVars,
    'diagNames': [],
}

profile_s  = 'profile'
sfc_s      = 'surface'
radiance_s = 'radiance'
model_s    = 'model'

## analysis groups for configuring AnalyzeStats + Analyses
convGrp = 'conv'
satwindGrp = 'satwind'
abiGrp = 'abi'
ahiGrp = 'ahi'
amsuaGrp = 'amsua'
amsuacldGrp = 'amsuacld'
mhsGrp = 'mhs'
iasiGrp = 'iasi'
modelGrp = 'model'

anGroupConfig = {
    convGrp: {'npwork': 18, 'npread': 36, 'analyze walltime': '00:20:00'},
    satwindGrp: {'npwork': 16, 'npread': 36, 'analyze walltime': '00:40:00'},
    abiGrp: {'npwork': 1, 'npread': 36, 'analyze walltime': '03:00:00'},
    ahiGrp: {'npwork': 1, 'npread': 36, 'analyze walltime': '03:00:00'},
    amsuaGrp: {'npwork': 14, 'npread': 36, 'analyze walltime': '00:60:00'}, #bg
    amsuacldGrp: {'npwork': 14, 'npread': 36, 'analyze walltime': '00:25:00'},
    mhsGrp: {'npwork': 14, 'npread': 36, 'analyze walltime': '00:20:00'},
    iasiGrp: {'npwork': 14, 'npread': 36, 'analyze walltime': '00:50:00'},
    modelGrp: {'npwork': 10, 'npread': 36, 'analyze walltime': '03:30:00'},
}

# Each entry of DiagSpaceConfig is a key-value pair, with the following possible values:
#DiagSpaceName (fromYAML):{
#    'DiagSpaceGrp': any of profile_s, radiance_s, model_s
#    'process': True or False
#    'anGrp': used to configure statistical analysis jobs and scripts
#    'binVarConfigs': binning variable configurations used when calculating statistics on diagnostics
#    'diagNames': list of selected diagnostic names,
#                 e.g., pconf.diffDiagnostics[+pconf.absDiagnostics][+pconf.cloudyRadDiagnostics]
#                 new diagNames should be defined in diag_utils
#    'channels': channel selection for plot_obs_nc_loc
#}

DiagSpaceConfig = {
#conventional
    'aircraft': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profilePressureBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssro': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': gnssroBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssrobndropp1d': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': gnssroBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssroref': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': gnssroBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'satwind': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': satwindGrp,
        'binVarConfigs': satwindBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
    'satwnd': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': satwindBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
    'sfc': {
        'DiagSpaceGrp': sfc_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': surfBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
    'sondes': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profilePressureBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
#radiances
    'abi_g16': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': abiGrp,
        'binVarConfigs': abi_g16_binVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.cloudyRadDiagnostics,
        'channels': range(7,17),
        ### example for channel selection/ordering at plotting phase:
        'analyzed channels': [8, 9, 10, 13, 14, 15],
    },
    'ahi_himawari8': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': ahiGrp,
        'binVarConfigs': ahi_himawari8_binVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.cloudyRadDiagnostics,
        'channels': range(7,17),
        ### example for channel selection/ordering at plotting phase:
        'analyzed channels': [8, 9, 10, 13, 14, 15],
    },
    'abi-clr_g16': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': abiGrp,
        'binVarConfigs': abi_g16_binVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'ahi-clr_himawari8': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': ahiGrp,
        'binVarConfigs': ahi_himawari8_binVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'airs_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [1,6,7],
    },
    'amsua_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [8,9],
    },
    'amsua_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [5,6,9],
    },
    'amsua_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [],
    },
    'amsua_n15': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [5,6,7,8,9],
    },
    'amsua_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [5,6,7,8,9],
    },
    'amsua_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [5,6,7,9],
    },
    'amsua_n19--hydro': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [1,2,3,15],
    },
    'amsua_n19--nohydro': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [4,5,6,7,9,10,11,12,13,14],
    },
    'amsua-cld_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n15': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polarcldBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'cris-fsr_npp': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [24,26,28,32,37,39],
    },
    'hirs4_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
    'iasi_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': iasiGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [16,104,154,254,410,1585,2019,2993,3110],
        'analyzed channels': [16,104,154,254,410,1585,2019,2993,3110],
    },
    'iasi_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': iasiGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [16,104,154,254,410,1585,2019,2993,3110],
        'analyzed channels': [16,104,154,254,410,1585,2019,2993,3110],
    },
    'iasi_metop-c': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'anGrp': iasiGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': [16,104,154,254,410,1585,2019,2993,3110],
        'analyzed channels': [16,104,154,254,410,1585,2019,2993,3110],
    },
    'mhs_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': range(1,6),
    },
    'mhs_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': range(1,6),
    },
    'mhs_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': range(1,6),
    },
    'mhs_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': polarBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics+pconf.nobcDiagnostics,
        'channels': range(1,6),
    },
    'seviri_m08': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [5],
    },
    'sndrd1_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
    'sndrd2_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
    'sndrd3_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
    'sndrd4_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
#models
   'mpas': {
        'DiagSpaceGrp': model_s,
        'process': True,
        'anGrp': modelGrp,
        'binVarConfigs': modelBinVars,
        'diagNames': pconf.modelDiags,
        ### examples for variable selection/ordering at plotting phase:
        ## all standard variables
        #'analyzed variables': [
        #  'T2m', 'Q2m', 'U10m', 'V10m', 'Ps',
        #  'T', 'Theta', 'Qv', 'U', 'V', 'P', 'rho', 'W',
        #],
        ## all standard variables + Qv by levels
        'analyzed variables': [
          'T2m', 'Q2m', 'U10m', 'V10m', 'Ps',
          'T', 'Theta', 'Qv', 'U', 'V', 'P', 'rho', 'W',
          'Qv01to10',
          'Qv11to20',
          'Qv21to30',
          'Qv31to40',
          'Qv41to55',
        ],
        ## specific standard variables
        #'analyzed variables': [
        #  'T2m', 'Q2m', 'U10m', 'Ps',
        #  'T', 'Qv', 'U', 'V'
        #],
        ## all 2D
        #'analyzed variables': [
        #  'T2m', 'Q2m', 'U10m', 'V10m', 'Ps',
        #]
        ## all standard 3D
        #'analyzed variables': [
        #  'T', 'Theta', 'Qv', 'U', 'V', 'P', 'rho', 'W',
        #],
        ## special Qv
        #'analyzed variables': [
        #  'Qv01to10',
        #  'Qv11to20',
        #  'Qv21to30',
        #  'Qv31to40',
        #  'Qv41to55',
        #],
    },
}

