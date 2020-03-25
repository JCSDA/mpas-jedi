import binning_utils as bu
import binning_configs as bcs
from copy import deepcopy
import diag_utils as du
import os
import var_utils as vu


#========================================================================
# Sub-selections of binning_configs.binVarConfigs for specific DiagSpaces
#========================================================================

#################################################################
## Generic binVarConfigs that apply to all observation categories
#################################################################
obsBinVars = {
    vu.obsVarQC: [bu.goodQCMethod,bu.badQCMethod],
    vu.obsVarLat: [bu.identityBinMethod,bu.latbandsMethod],
    vu.obsVarLT: [bu.identityBinMethod],
    vu.obsVarNormErr: [bu.identityBinMethod],
    'ObsRegion': ['CONUS'],
}


################################
## binVarConfigs for surface obs
################################
surfBinVars = deepcopy(obsBinVars)


##########################################################
## binVarConfigs for profile obs w/ pressure vertical bins
##########################################################
profPressBinVars = deepcopy(obsBinVars)
profPressBinVars[vu.obsVarPrs] = [bu.identityBinMethod,bu.PjetMethod]
profPressBinVars[vu.obsVarLat].append(bu.PjetMethod)

# 2D pressure bins with named latitude-band methods
for latBand in bcs.namedLatBands['values']:
    profPressBinVars[vu.obsVarPrs].append(latBand)


##########################################################
## binVarConfigs for profile obs w/ altitude vertical bins
##########################################################
profAltBinVars = deepcopy(obsBinVars)
profAltBinVars[vu.obsVarAlt] = [bu.identityBinMethod,bu.altjetMethod]
profAltBinVars[vu.obsVarLat].append(bu.altjetMethod)

# 2D altitude bins with named latitude-band methods
for latBand in bcs.namedLatBands['values']:
    profAltBinVars[vu.obsVarAlt].append(latBand)


#################################
## binVarConfigs for radiance obs
#################################
radianceBinVars = deepcopy(obsBinVars)
radianceBinVars[vu.obsVarGlint] = [bu.identityBinMethod]
radianceBinVars[vu.obsVarSenZen] = [bu.identityBinMethod]


#########################################
## binVarConfigs for geostationary IR obs
## i.e., GOES-ABI, Himawari-AHI
#########################################
geoirBinVars = deepcopy(radianceBinVars)
geoirBinVars[vu.obsVarACI] = [bu.identityBinMethod]
geoirBinVars[vu.obsVarCldFrac] = [bu.identityBinMethod,
                                  bu.cloudbandsMethod]


# Binning variables with clr-/cld-sky methods
for var in bcs.cldfracBinVars.keys():
    geoirBinVars[var].append(bu.clrskyMethod)
    geoirBinVars[var].append(bu.cldskyMethod)


# 2D sensor zenith bins with named latitude-band methods and clear profiles
for var in bcs.clrlatBinVars.keys():
    for latBand in bcs.namedLatBands['values']:
        geoirBinVars[var].append(bcs.clrlatMethods[latBand])


# symmetric cloud impact (expensive)
selectSCIMethods = [
    bu.OkamotoMethod,
    bu.ScaleOkamotoMethod,
#    bu.ModHarnischMethod,
#    bu.ScaleModHarnischMethod,
]

geoirBinVars[vu.obsVarSCI] = []
for method in selectSCIMethods:
    geoirBinVars[vu.obsVarSCI].append(method)
    geoirBinVars[vu.obsVarNormErr].append(method)


#########################################
# binVarConfigs for model space variables
#########################################
#modelBinVars = { 'ModelLatBand':  [bu.latbandsMethod,bu.identityBinMethod]
#               , 'ModelBox':      ['CONUS']
#               , 'ModelAltitude': [bu.identityBinMethod]
#               , 'ModelPressure': [bu.identityBinMethod]
#               }

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
radiance_s = 'radiance'
model_s    = 'MPAS'

# Each entry of DiagSpaceConfig is a key-value pair, with the following possible values:
#DiagSpaceName (fromYAML):{
#    'DiagSpaceGrp': any of profile_s, radiance_s, model_s
#    'process': True or False
#    'binVarConfigs': binning variable configurations used when calculating statistics on diagnostics
#    'diagNames': list of selected diagnostic names,
#                 e.g., du.diffDiagNames[+du.absDiagNames][+du.cloudyRadDiagNames]
#                 new diagNames should be defined in diag_utils
#    'channels': channel selection for plot_obs_nc_loc
#}

defaultDiags = du.diffDiagNames

DiagSpaceConfig = {
#conventional
    'aircraft': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profPressBinVars,
        'diagNames': defaultDiags,
    },
    'gnssro': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profAltBinVars,
        'diagNames': du.relDiagNames,
    },
    'gnssrobndropp1d': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profAltBinVars,
        'diagNames': du.relDiagNames,
    },
    'gnssroref': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profAltBinVars,
        'diagNames': du.relDiagNames,
    },
    'satwind': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profPressBinVars,
        'diagNames': defaultDiags,
    },
    'sondes': {
        'DiagSpaceGrp': profile_s,
        'process': True,
        'binVarConfigs': profPressBinVars,
        'diagNames': defaultDiags,
    },
#radiances
    'abi_g16': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': geoirBinVars,
        'diagNames': defaultDiags,
#        'diagNames': du.diffDiagNames+du.absDiagNames+du.cloudyRadDiagNames,
        'channels': [8,9,10,11,13,14,15,16],
    },
    'airs_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [1,6,7],
    },
    'amsua_n15': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [5,6,7,8,9],
    },
    'amsua_n18': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [5,6,7,8,9],
    },
    'amsua_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [5,6,7,9],
    },
    'amsua_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [5,6,9],
    },
    'amsua_metop-b': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [],
    },
    'amsua_aqua': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [8,9],
    },
    'amsua_n19--hydro': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [1,2,3,15],
    },
    'amsua_n19--nohydro': {
        'DiagSpaceGrp': radiance_s,
        'process': True,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [4,5,6,7,9,10,11,12,13,14],
    },
    'cris-fsr_npp': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [24,26,28,32,37,39],
    },
    'hirs4_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,16),
    },
    'iasi_metop-a': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [16,29,32,35,38,41,44],
    },
    'mhs_n19': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,6),
    },
    'seviri_m08': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': [5],
    },
    'sndrd1_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,16),
    },
    'sndrd2_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,16),
    },
    'sndrd3_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,16),
    },
    'sndrd4_g15': {
        'DiagSpaceGrp': radiance_s,
        'process': False,
        'binVarConfigs': radianceBinVars,
        'diagNames': defaultDiags,
        'channels': range(1,16),
    },
#models
#   'mpas-gfs': {
#        'DiagSpaceGrp': model_s,
#        'process': False,
#        'binVarConfigs': modelBinVars,
#    },
}

