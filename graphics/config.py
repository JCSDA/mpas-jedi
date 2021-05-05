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
obsBinVars = defaultdict(list)
obsBinVars[vu.obsVarQC] += [bu.goodQCMethod, bu.badQCMethod, bu.allQCMethod]
obsBinVars[vu.obsVarLat] += [bu.identityBinMethod, bu.latbandsMethod]
obsBinVars[vu.obsVarLT] += [bu.identityBinMethod]
obsBinVars[vu.obsVarNormErr] += [bu.identityBinMethod]
obsBinVars[vu.obsRegionBinVar] += ['CONUS']
obsBinVars[vu.obsRegionBinVar] += [bu.geoirlatlonboxMethod]


################################
## binVarConfigs for surface obs
################################
surfBinVars = deepcopy(obsBinVars)


##########################################################
## binVarConfigs for profile obs w/ pressure vertical bins
##########################################################
profPressBinVars = deepcopy(obsBinVars)
profPressBinVars[vu.obsVarPrs] += [bu.identityBinMethod, bu.PjetMethod]
profPressBinVars[vu.obsVarLat] += [bu.PjetMethod]

# 2D pressure bins with named latitude-band methods
for latBand in pconf.namedLatBands['values']:
    profPressBinVars[vu.obsVarPrs] += [latBand]


##########################################################
## binVarConfigs for profile obs w/ altitude vertical bins
##########################################################
profAltBinVars = deepcopy(obsBinVars)
profAltBinVars[vu.obsVarAlt] += [bu.identityBinMethod, bu.altjetMethod]
profAltBinVars[vu.obsVarLat] += [bu.altjetMethod]

# 2D altitude bins with named latitude-band methods
for latBand in pconf.namedLatBands['values']:
    profAltBinVars[vu.obsVarAlt] += [latBand]


#################################
## binVarConfigs for radiance obs
#################################
radianceBinVars = deepcopy(obsBinVars)
radianceBinVars[vu.obsVarGlint] += [bu.identityBinMethod]
radianceBinVars[vu.obsVarSenZen] += [bu.identityBinMethod]


# binBYLandFrac controls whether to bin by land-ocean fractions and categories
# NOTE: requires geovals to be written for applicable DiagSpaces
binBYLandFrac = False
if binBYLandFrac:
    radianceBinVars[vu.obsVarLandFrac] += [bu.identityBinMethod,
                                           bu.surfbandsMethod]


############################################
## binVarConfigs for polar-orbitting MW obs,
## including cloud metrics
## e.g., AMSUA
############################################
polmwBinVars = deepcopy(radianceBinVars)
polmwBinVars[vu.obsVarACI] += [bu.identityBinMethod]
polmwBinVars[vu.obsVarSCI] += [bu.OkamotoMethod]


##########################################
## binVarConfigs for geostationary IR obs
## including cloud metrics
## e.f., GOES-ABI, Himawari-AHI
##########################################
geoirBinVars = deepcopy(radianceBinVars)
geoirBinVars[vu.obsVarACI] += [bu.identityBinMethod]
geoirBinVars[vu.obsVarCldFrac] += [bu.identityBinMethod,
                                  bu.cloudbandsMethod]

# Binning variables with clr-/cld-sky methods
for var in pconf.cldfracBinVars.keys():
    if var != vu.obsVarLandFrac or binBYLandFrac:
        geoirBinVars[var] += [bu.clrskyMethod]
        geoirBinVars[var] += [bu.cldskyMethod]

# Binning variables with land category methods
if binBYLandFrac:
    for var in pconf.landfracBinVars.keys():
        for surfBand in pconf.namedLandFracBands['values']:
            geoirBinVars[var] += [surfBand]

# symmetric cloud impact (expensive)
geoirSCIMethods = [
    bu.OkamotoMethod,
    bu.ScaleOkamotoMethod,
#    bu.ModHarnischMethod,
#    bu.ScaleModHarnischMethod,
]

for method in geoirSCIMethods:
    geoirBinVars[vu.obsVarSCI] += [method]
# uncomment to diagnose (HofX - ObsValue) / ObsError for SCI-parameterized ObsError
# note: requires SCIErrParams to be defined for all all DiagSpaces that use it
#    geoirBinVars[vu.obsVarNormErr] += [method]


#########################################
# binVarConfigs for model space variables
#########################################
modelBinVars = defaultdict(list)
modelBinVars[vu.noBinVar] += [bu.noBinMethod]
modelBinVars[vu.modVarLat] += [bu.identityBinMethod]
modelBinVars[vu.modVarLat] += [bu.latbandsMethod]
modelBinVars[vu.modVarLev] += [bu.identityBinMethod]
# 2D model level bins with named latitude-band methods
for latBand in pconf.namedLatBands['values']:
  modelBinVars[vu.modVarLev] += [latBand]
modelBinVars[vu.modelRegionBinVar] += [bu.geoirlatlonboxMethod]

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
abiGrp = 'abi'
ahiGrp = 'ahi'
amsuaGrp = 'amsua'
amsuacldGrp = 'amsuacld'
mhsGrp = 'mhs'
modelGrp = 'model'

anGroupConfig = {
    convGrp: {'npwork': 36, 'npread': 36, 'analyze_walltime': '00:10:00'},
    abiGrp: {'npwork': 12, 'npread': 36, 'analyze_walltime': '01:00:00'},
    ahiGrp: {'npwork': 12, 'npread': 36, 'analyze_walltime': '01:00:00'},
    amsuaGrp: {'npwork': 36, 'npread': 36, 'analyze_walltime': '00:12:00'},
    amsuacldGrp: {'npwork': 36, 'npread': 36, 'analyze_walltime': '00:12:00'},
    mhsGrp: {'npwork': 36, 'npread': 36, 'analyze_walltime': '00:08:00'},
    modelGrp: {'npwork': 10, 'npread': 36, 'analyze_walltime': '03:30:00'},
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
        'binVarConfigs': profPressBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssro': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profAltBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssrobndropp1d': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profAltBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'gnssroref': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profAltBinVars,
        'diagNames': pconf.relativeDiagnostics+pconf.sigmaDiagnostics,
    },
    'satwind': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'anGrp': convGrp,
        'binVarConfigs': profPressBinVars,
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
        'binVarConfigs': profPressBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
    },
#radiances
    'abi_g16': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': abiGrp,
        'binVarConfigs': geoirBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'ahi_himawari8': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': ahiGrp,
        'binVarConfigs': geoirBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'abi-clr_g16': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': abiGrp,
        'binVarConfigs': geoirBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'ahi-clr_himawari8': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': ahiGrp,
        'binVarConfigs': geoirBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'airs_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [1,6,7],
    },
    'amsua_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [8,9],
    },
    'amsua_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [5,6,9],
    },
    'amsua_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [],
    },
    'amsua_n15': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [5,6,7,8,9],
    },
    'amsua_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [5,6,7,8,9],
    },
    'amsua_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuaGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [5,6,7,9],
    },
    'amsua_n19--hydro': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [1,2,3,15],
    },
    'amsua_n19--nohydro': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [4,5,6,7,9,10,11,12,13,14],
    },
    'amsua-cld_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n15': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'amsua-cld_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': amsuacldGrp,
        'binVarConfigs': polmwBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': [1,2,3,4,15],
    },
    'cris-fsr_npp': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [24,26,28,32,37,39],
    },
    'hirs4_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': range(1,16),
    },
    'iasi_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics,
        'channels': [16,29,32,35,38,41,44],
    },
    'mhs_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': range(1,6),
    },
    'mhs_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': range(1,6),
    },
    'mhs_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
        'channels': range(1,6),
    },
    'mhs_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'anGrp': mhsGrp,
        'binVarConfigs': radianceBinVars,
        'diagNames': pconf.defaultDiagnostics+pconf.sigmaDiagnostics,
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
    },
}

