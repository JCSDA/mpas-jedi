#!/usr/bin/env python3

from copy import deepcopy
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

# NC observation variable name substitutions
vNameStr = 'varName'
vChanStr = 'varCHAN'

# MPAS-JEDI suffixes for observation-type variables
hofxGroup   = 'hofx'
depbgGroup  = 'depbg'
depanGroup  = 'depan'
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

# functions for extracting sub-parts of UFO variable names
def splitObsVarGrp(varATgroup):
    if "@" in varATgroup:
        var = ''.join(varATgroup.split("@")[:-1])
        grp = ''.join(varATgroup.split("@")[-1])
    else:
        var = varATgroup
        grp = miss_s
    return var, grp


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
def base2dbVar(baseVar,varName,Iter):
    iterStr = str(Iter)
    dictName, suf = splitIntSuffix(varName)
    dbVar = re.sub(vNameStr,varName,baseVar)
    dbVar = re.sub(vChanStr,suf,dbVar)
    dbVar = re.sub(hofxGroup,hofxGroup+iterStr,dbVar)
    dbVar = re.sub(errorGroup,errorGroup+iterStr,dbVar)
    dbVar = re.sub(qcGroup,qcGroup+iterStr,dbVar)
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


#misc. constants; TODO: collect into single script
deg2rad = np.pi / np.float(180.0)
rad2deg = np.float(180.0) / np.pi

