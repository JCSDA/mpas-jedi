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

## spread diagnostic variable names
EddT = 'E[dd^T]'
HBHT = 'HBH^T'
R = 'R'
HBHTplusR = '[HBH^T+R]'

DiagnosticVars = {}
DiagnosticVars[EddT] = 'InnovSpread'
DiagnosticVars[HBHT] = 'EnsembleSpread'
DiagnosticVars[R] = 'ObsError'
DiagnosticVars[HBHTplusR] = 'TotalSpread'

## NC variable names for MPAS-JEDI
obsVarAlt     = 'altitude'
obsVarACI     = 'asymmetric_cloud_impact'
obsVarBT      = 'brightness_temperature'
obsVarY       = 'y'
obsVarH       = 'h'
obsVarBTClear = obsVarBT+'_assuming_clear_sky'
obsVarCldFracY = 'cloud_area_fraction'
obsVarCldFracX = 'cloud_area_fraction_in_atmosphere_layer'
obsVarDT      = 'datetime'
obsVarGlint   = 'glint'
obsVarLT      = 'LocalTime'
obsVarLandFrac= 'land_area_fraction'
obsVarLat     = 'latitude'
obsVarLon     = 'longitude'
obsVarDep     = 'y-h'
obsVarClearSkyDep = 'y-hclr'
obsVarLogDepRatio = 'logDeparture'
obsVarNormDep = 'dσ\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
obsVarNormDep = 'ENI' #error-normalized innovation
obsVarPrs     = 'air_pressure'
obsVarQC      = 'QCflag'
obsVarCI      = 'cloud_impact'
obsVarLogCI   = 'logarithm_cloud_impact'
obsVarMCI     = 'modeled_cloud_impact'
obsVarOCI     = 'observed_cloud_impact'
obsVarSenZen  = 'sensor_zenith_angle'
obsVarSenAzi  = 'sensor_azimuth_angle'
obsVarSolZen  = 'solar_zenith_angle'
obsVarSolAzi  = 'solar_azimuth_angle'

degree= u'\N{DEGREE SIGN}'

# columns: var_name            unit_used   abbr.
varDictObs = {
    'air_temperature':        [ 'K',     'T'       ],
    'bending_angle':          [ '%',     'Bnd'     ],
    obsVarBT:                 [ 'K',     'BT'      ],
    obsVarY:                  [ 'K',     obsVarY   ],
    obsVarH:                  [ 'K',     obsVarH   ],
    'eastward_wind':          [ 'm/s',   'U'       ],
    'northward_wind':         [ 'm/s',   'V'       ],
    'refractivity':           [ '%',     'Ref'     ],
    'specific_humidity':      [ 'kg/kg', 'qv'      ],
    'surface_pressure':       [ 'Pa',    'Ps'      ],
    'virtual_temperature':    [ 'K',     'Tv'      ],
    obsVarAlt:                [ 'm',     'alt'     ],
    obsVarACI:                [ 'K',     'ACI'     ],
    obsVarCldFracY:           [ miss_s,  'CFy'     ],
    obsVarCldFracX:           [ miss_s,  'CFx'     ],
    obsVarLandFrac:           [ miss_s,  'landfrac'],
    obsVarLat:                [ degree,  'lat'     ],
    obsVarLon:                [ degree,  'lon'     ],
    obsVarLT:                 [ 'hr',    obsVarLT  ],
    'Normalized departure':   [ miss_s,  obsVarNormDep ],
    obsVarPrs:                [ 'hPa',   'P'       ],
    obsVarQC:                 [ miss_s,  obsVarQC  ],
    obsVarCI:                 [ 'K',     'CI'     ],
    obsVarLogCI:              [ 'K',     'Log-CI'  ],
    obsVarMCI:                [ 'K',     'MCI'     ],
    obsVarOCI:                [ 'K',     'OCI'     ],
    obsVarSenZen:             [ degree,  'zenith'  ],
    obsVarGlint:              [ degree,  obsVarGlint  ],
}
#Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == obsVarAlt:

obsRegionBinVar = 'ObsRegion'
varDictObs[obsRegionBinVar] = [miss_s, obsRegionBinVar]

# IODA observation variable name substitutions
vNameStr = 'varName'
vChanStr = 'varCHAN'

# Dynamic Fixed IODA ObsGroups for observation-type variables
hofxGroup   = 'hofx'
errorGroup  = 'EffectiveError'
qcGroup     = 'EffectiveQC'
bcGroup     = 'ObsBias'
bgIter = '0'

# Fixed IODA ObsGroups for observation-type variables
obsGroup    = 'ObsValue'
metaGroup   = 'MetaData'

# Dynamic MPAS-Workflow ObsGroups for observation-type variables
depGroup    = 'depIter'

# GeoVaLs ObsGroups only used in post-processing, not in files
diagGroup   = 'ObsDiag'
geoGroup    = 'GeoVaLs'

# Generic variable names used as baseVar
selfObsValue = 'selfObsValue'
selfDepValue = 'selfDepValue'
selfHofXValue = 'selfHofXValue'
selfErrorValue = 'selfErrorValue'
selfQCValue = 'selfQCValue'
selfBCValue = 'selfBCValue'
bgHofXValue = 'bgHofXValue'
altMeta = 'altMeta'
cldfracMeta = 'cldfracMeta'
datetimeMeta = 'datetimeMeta'
latMeta = 'latMeta'
lonMeta = 'lonMeta'
prsMeta = 'prsMeta'
senzenMeta = 'senzenMeta'
senaziMeta = 'senaziMeta'
solzenMeta = 'solzenMeta'
solaziMeta = 'solaziMeta'
landfracGeo = 'landfracGeo'
clrskyBTDiag = 'clrskyBTDiag'
cldfracGeo = 'cldfracGeo'

# Context-dependent (dynamic) IODA variable names
ObsGroups = {}
ObsVars = {}
ObsGroups[selfObsValue] = obsGroup
ObsGroups[selfDepValue] = depGroup
ObsGroups[selfHofXValue] = hofxGroup
ObsGroups[selfErrorValue] = errorGroup
ObsGroups[selfQCValue] = qcGroup
ObsGroups[selfBCValue] = bcGroup

for key in ObsGroups.keys():
  ObsVars[key] = vNameStr

ObsGroups[bgHofXValue] = hofxGroup+bgIter
ObsVars[bgHofXValue] = vNameStr

# Fixed IODA variable names (MetaData)
ObsVars[altMeta]      = obsVarAlt
ObsVars[cldfracMeta]  = obsVarCldFracY
ObsVars[datetimeMeta] = obsVarDT
ObsVars[latMeta]      = obsVarLat
ObsVars[lonMeta]      = obsVarLon
ObsVars[prsMeta]      = obsVarPrs
ObsVars[senzenMeta]   = obsVarSenZen
ObsVars[senaziMeta]   = obsVarSenAzi
ObsVars[solzenMeta]   = obsVarSolZen
ObsVars[solaziMeta]   = obsVarSolAzi

for key in ObsVars.keys():
  if 'Meta' in key:
    ObsGroups[key] = metaGroup

# separator to be used for channel or other integer suffixes
intSufSeparator = '_'


# GeoVaLs variable names
ObsVars[landfracGeo] = obsVarLandFrac
ObsGroups[landfracGeo] = geoGroup

ObsVars[cldfracGeo] = obsVarCldFracX
ObsGroups[cldfracGeo] = geoGroup

ObsVars[clrskyBTDiag] = obsVarBTClear+intSufSeparator+vChanStr
ObsGroups[clrskyBTDiag] = diagGroup


# ensemble/mean classifiers
mean = 'mean'
ensemble = 'ensemble'
ensSuffixBase = "&&&mem"
def ensSuffix(member):
  if member == 0:
    return ""
  else:
    return ensSuffixBase+str(member)


# functions for extracting/combining sub-parts of UFO variable names
def splitObsVarGrp(WholeVarGrp):
  if "@" in WholeVarGrp:
    var = '@'.join(WholeVarGrp.split('@')[:-1])
    grp = WholeVarGrp.split('@')[-1]
  elif "/" in WholeVarGrp:
    grp = WholeVarGrp.split('/')[0]
    var = '/'.join(WholeVarGrp.split('/')[1:])
  else:
    var = WholeVarGrp
    grp = miss_s
  return var, grp


def splitIntSuffix(var):
  # separate integer suffixes (e.g., brightness_temperature_*)
  obsVarName, grpName = splitObsVarGrp(var)
  suf = obsVarName.split(intSufSeparator)[-1]
  if not pu.isint(suf):
    suf = ''
  else:
    obsVarName = intSufSeparator.join(obsVarName.split(intSufSeparator)[:-1])
  return obsVarName, suf


def appendSuffix(var, suf):
  return var+intSufSeparator+str(suf)


def obsVarAttributes(var):
  # return short name and units
  dictName, suf = splitIntSuffix(var)
  varAtt = varDictObs.get(dictName,[miss_s,dictName])
  varShort = varAtt[1]+suf
  varUnits = varAtt[0]
  return varShort, varUnits




## NC variable names for MPAS-Model
#modVarAlt = 'zgrid' # --> needs to be interpolated to nVertLevels instead of nVertLevelsP1
modVarPrs = 'pressure'
modVarLat = 'latCell'
modVarLon = 'lonCell'
modVarLev = 'modelLevel'

levModel = 'levModel'
latModel = 'latModel'
lonModel = 'lonModel'

AllVars = deepcopy(ObsVars)
AllGroups = deepcopy(ObsGroups)

AllVars[levModel] = modVarLev
AllVars[latModel] = modVarLat
AllVars[lonModel] = modVarLon

AllGroups[levModel] = None
AllGroups[latModel] = None
AllGroups[lonModel] = None

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
  'qv01to30':               [ 'g/kg', 'Qv01to30' ],
  'qv01to10':               [ 'g/kg', 'Qv01to10' ],
  'qv11to20':               [ 'g/kg', 'Qv11to20' ],
  'qv21to30':               [ 'g/kg', 'Qv21to30' ],
  'qv31to40':               [ 'g/kg', 'Qv31to40' ],
  'qv41to55':               [ 'g/kg', 'Qv41to55' ],
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

modVarNames2d = [
  'q2',
  'surface_pressure',
  't2m',
  'u10',
  'v10',
]

modVarNamesBase3d = [
  'pressure',
  'qv',
  'rho',
  'theta',
  'temperature',
  'uReconstructMeridional',
  'uReconstructZonal',
  'w',
]

modVarNames3d = modVarNamesBase3d+[
  #extra variables
  'qv01to30',
  'qv01to10',
  'qv11to20',
  'qv21to30',
  'qv31to40',
  'qv41to55',
]

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


def varAttributes(var):
  # return short name and units
  dictName, suf = splitIntSuffix(var)
  varAtt = varDictAll.get(dictName,[miss_s,dictName])
  varShort = varAtt[1]+suf
  varUnits = varAtt[0]
  return varShort, varUnits


# FileFormat-dependent variable name constructors
def groupSLASHvar(var, group):
  vv = var
  if group is not None:
    vv = group+'/'+vv
  return vv

def varATgroup(var, group):
  vv = var
  if group is not None:
    vv += '@'+group
  return vv

def rawVar(var, group=None):
  return var

ncFileFormat = 'nc'
hdfFileFormat = 'hdf'
modelFileFormat = 'model'
AllVarCtors = {
  ncFileFormat: varATgroup,
  hdfFileFormat: groupSLASHvar,
  modelFileFormat: rawVar,
}


#BaseVars describes the file-format-specific generic variable names

BaseVars = {}
for fileFormat, ctor in AllVarCtors.items():
  BaseVars[fileFormat] = {}
  for baseVar, ActualVar in AllVars.items():
    BaseVars[fileFormat][baseVar] = ctor(ActualVar, AllGroups[baseVar])

def base2dbVar(baseVar, varName, fileFormat, outerIter = None):
  # converts baseVar to a context-specific variable name to retrieve from a JediDB object
  dbVar = BaseVars[fileFormat][baseVar]

  dictName, suf = splitIntSuffix(varName)
  dbVar = re.sub(vNameStr,varName,dbVar)
  dbVar = re.sub(vChanStr,suf,dbVar)

  if outerIter is None:
    iterStr = ''
  else:
    iterStr = str(outerIter)
  for group in [hofxGroup, errorGroup, qcGroup, bcGroup]:
    # append iterStr if one is not already appended
    if group in dbVar.split('@') or group in dbVar.split('/'):
      dbVar = re.sub(group,group+iterStr,dbVar)
  if iterStr != '':
    if iterStr == bgIter:
      dbVar = re.sub(depGroup,depbgGroup,dbVar)
    else:
      dbVar = re.sub(depGroup,depanGroup,dbVar)
  return dbVar
