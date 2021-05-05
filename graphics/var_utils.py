#!/usr/bin/env python3

from copy import deepcopy
from jediApplicationArgs import depbgGroup, depanGroup
import numpy as np
import os
import re
import plot_utils as pu

miss_f = -88888.8
miss_i = -88888
miss_s = 'null'
csvSEP = ';'

#=====================
# variable definitions
#=====================

## NC variable names for MPAS-JEDI
obsVarAlt     = 'altitude'
obsVarACI     = 'asymmetric_cloud_impact'
obsVarCldFrac = 'cloud_area_fraction'
obsVarDT      = 'datetime'
obsVarLT      = 'LocalTime'
obsVarLandFrac= 'land_area_fraction'
obsVarLat     = 'latitude'
obsVarLon     = 'longitude'
obsVarNormErr = 'dσ\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
obsVarPrs     = 'air_pressure'
obsVarQC      = 'QCflag'
obsVarSCI     = 'symmetric_cloud_impact'
obsVarSenZen  = 'sensor_zenith_angle'
obsVarSenAzi  = 'sensor_azimuth_angle'
obsVarSolZen  = 'solar_zenith_angle'
obsVarSolAzi  = 'solar_azimuth_angle'
obsVarGlint   = 'glint'

degree= u'\N{DEGREE SIGN}'

# columns: var_name            unit_used   abbr.
varDictObs = {
    'air_temperature':        [ 'K',     'T'       ],
    'bending_angle':          [ '%',     'Bnd'     ],
    'brightness_temperature': [ 'K',     'BT'      ],
    'eastward_wind':          [ 'm/s',   'U'       ],
    'northward_wind':         [ 'm/s',   'V'       ],
    'refractivity':           [ '%',     'Ref'     ],
    'specific_humidity':      [ 'kg/kg', 'qv'      ],
    'surface_pressure':       [ 'Pa',    'Ps'      ],
    'virtual_temperature':    [ 'K',     'Tv'      ],
    obsVarAlt:                [ 'm',     'alt'     ],
    obsVarACI:                [ 'K',     'ACI'     ],
    obsVarCldFrac:            [ miss_s,  'cldfrac' ],
    obsVarLandFrac:           [ miss_s,  'landfrac'],
    obsVarLat:                [ degree,  'lat'     ],
    obsVarLon:                [ degree,  'lon'     ],
    obsVarLT:                 [ 'hr',    obsVarLT  ],
    obsVarNormErr:            [ miss_s,  obsVarNormErr ],
    obsVarPrs:                [ 'hPa',   'P'       ],
    obsVarQC:                 [ miss_s,  obsVarQC  ],
    obsVarSCI:                [ 'K',     'SCI'     ],
    obsVarSenZen:             [ degree,  'zenith'  ],
    obsVarGlint:              [ degree,  obsVarGlint  ],
}
#Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == obsVarAlt:

obsRegionBinVar = 'ObsRegion'
varDictObs[obsRegionBinVar] = [miss_s, obsRegionBinVar]

# NC observation variable name substitutions
vNameStr = 'varName'
vChanStr = 'varCHAN'

# MPAS-JEDI suffixes for observation-type variables
hofxGroup   = 'hofx'
depGroup    = 'depIter'
diagGroup   = 'ObsDiag'
geoGroup    = 'GeoVaLs'
obsGroup    = 'ObsValue'
metaGroup   = 'MetaData'
qcGroup     = 'EffectiveQC'
errorGroup  = 'EffectiveError'

# Fixed NC variable names
selfObsValue   = vNameStr+'@'+obsGroup
#selfDepBGValue = vNameStr+'@'+depbgGroup
#selfDepANValue = vNameStr+'@'+depanGroup
selfDepValue   = vNameStr+'@'+depGroup

selfHofXValue  = vNameStr+'@'+hofxGroup
selfQCValue    = vNameStr+'@'+qcGroup
selfErrorValue = vNameStr+'@'+errorGroup

altMeta      = obsVarAlt+'@'+metaGroup
cldfracMeta  = obsVarCldFrac+'@'+metaGroup
dtMeta       = obsVarDT+'@'+metaGroup
latMeta      = obsVarLat+'@'+metaGroup
lonMeta      = obsVarLon+'@'+metaGroup
prsMeta      = obsVarPrs+'@'+metaGroup
senzenMeta   = obsVarSenZen+'@'+metaGroup
senaziMeta   = obsVarSenAzi+'@'+metaGroup
solzenMeta   = obsVarSolZen+'@'+metaGroup
solaziMeta   = obsVarSolAzi+'@'+metaGroup

landfracGeo = obsVarLandFrac+'@'+geoGroup

clrskyBTDiag = 'brightness_temperature_assuming_clear_sky_'+vChanStr+'@'+diagGroup

mean = 'mean'
ensemble = 'ensemble'
ensSuffixBase = "&&&mem"
def ensSuffix(member):
    if member == 0:
        return ""
    else:
        return ensSuffixBase+str(member)


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


def splitIntSuffix(var):
    # separate integer suffixes (e.g., brightness_temperature_*)
    obsVarName, grpName = splitObsVarGrp(var)
    suf = ''.join(obsVarName.split("_")[-1:])
    if not pu.isint(suf):
        suf = ''
    else:
        obsVarName = '_'.join(obsVarName.split("_")[:-1])
    return obsVarName, suf


def varAttributes(var):
    # return short name and units
    dictName, suf = splitIntSuffix(var)
    varAtt = varDictObs.get(dictName,[miss_s,dictName])
    varShort = varAtt[1]+suf
    varUnits = varAtt[0]
    return varShort, varUnits


bgIter = '0'

def base2dbVar(baseVar, varName, outerIter = None):
    if outerIter is None:
        iterStr = ''
    else:
        iterStr = str(outerIter)
    dictName, suf = splitIntSuffix(varName)
    dbVar = re.sub(vNameStr,varName,baseVar)
    dbVar = re.sub(vChanStr,suf,dbVar)
    for group in [hofxGroup, errorGroup, qcGroup]:
        # append iterStr if one is not already appended
        if group in dbVar.split('@'):
            dbVar = re.sub(group,group+iterStr,dbVar)
    if iterStr != '':
        if iterStr == bgIter:
            dbVar = re.sub(depGroup,depbgGroup,dbVar)
        else:
            dbVar = re.sub(depGroup,depanGroup,dbVar)
    return dbVar


## NC variable names for MPAS-Model
#modVarAlt = 'zgrid' # --> needs to be interpolated to nVertLevels instead of nVertLevelsP1
modVarPrs = 'pressure_p'
modVarLat = 'latCell'
modVarLon = 'lonCell'
modVarLev = 'model_level'

kgm3 = 'kg/m\N{SUPERSCRIPT THREE}'

# columns: var_name            unit_used   abbr.
varDictModel = {
  modVarLev:                [ miss_s, 'ModLev'],
  modVarLat:                [ degree, 'lat'  ],
  modVarLon:                [ degree, 'lon'  ],
  modVarPrs:                [ 'Pa',   'PP'   ],
  'pressure':               [ 'Pa',   'P'    ],
  'q2':                     [ 'g/kg', 'Q2m'  ],
  'qv':                     [ 'g/kg', 'Qv'   ],
  'rho':                    [ kgm3,   'rho'  ],
  'surface_pressure':       [ 'Pa',   'Ps'   ],
  't2m':                    [ 'C',    'T2m'  ],
  'temperature':            [ 'C',    'T'    ],
  'theta':                  [ 'K',    'Theta'],
  'u':                      [ 'm/s',  'uedge'],
  'u10':                    [ 'm/s',  'U10m' ],
  'uReconstructZonal':      [ 'm/s',  'U'    ],
  'uReconstructMeridional': [ 'm/s',  'V'    ],
  'v10':                    [ 'm/s',  'V10m' ],
  'w':                      [ 'm/s',  'W'    ],
}
#Note, qv unit is kg/kg in original mpas restart file. The unit is converted to g/kg when read qv.

## add dummy variable for no binning
noBinVar = 'all'
varDictModel[noBinVar] = [miss_s, noBinVar]

modelRegionBinVar = 'ModelRegion'
varDictModel[modelRegionBinVar] = [miss_s, modelRegionBinVar]

modVarNames2d = ['t2m','surface_pressure','q2','u10','v10']
modVarNames3d = ['theta','temperature','rho','pressure','uReconstructZonal','uReconstructMeridional','qv','w']

def modelVarAttributes(var):
    # return short name and units
    dictName, suf = splitIntSuffix(var)
    varAtt = varDictModel.get(dictName,[miss_s,dictName])
    varShort = varAtt[1]+suf
    varUnits = varAtt[0]
    return varShort, varUnits

#misc. constants; TODO: collect into single script
deg2rad = np.pi / np.float(180.0)
rad2deg = np.float(180.0) / np.pi


# dictionary containing all analyzed variables
varDictAll = deepcopy(varDictObs)
for var, desc in varDictModel.items():
  if var not in varDictAll:
    varDictAll[var] = deepcopy(desc)
  else:
    assert desc[0] == varDictAll[var][0], var+' units differ between varDictObs and varDictModel'
    assert desc[1] == varDictAll[var][1], var+' abbreviation differs between varDictObs and varDictModel'
