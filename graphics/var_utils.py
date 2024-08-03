#!/usr/bin/env python3

from copy import deepcopy
from jediApplicationArgs import ombgGroup, omanGroup, nOuterIter
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
DiagnosticVars[EddT] = 'Innov'
DiagnosticVars[HBHT] = 'EnsembleSpread'
DiagnosticVars[R] = 'ObsError'
DiagnosticVars[HBHTplusR] = 'TotalSpread'

## Variable names for MPAS-JEDI
# IODA/UFO variables (ObsSpace, GeoVaLs, ObsDiagnostics)
obsVarAlt     = 'height'
obsVarBT      = 'brightnessTemperature'
obsVarBTClear = 'brightness_temperature_assuming_clear_sky' # GeoVar
obsVarCldFracY = 'cloudAmount'
obsVarCldFracX = 'cloud_area_fraction_in_atmosphere_layer' # GeoVar
obsVarLandFrac= 'land_area_fraction'#'landAreaFraction'
obsVarLat     = 'latitude'
obsVarLon     = 'longitude'
obsVarPrs     = 'pressure'
obsVarSenZen  = 'sensorZenithAngle'
obsVarSenAzi  = 'sensorAzimuthAngle'
obsVarSolZen  = 'solarZenithAngle'
obsVarSolAzi  = 'solarAzimuthAngle'

# Post-processing derived variables
obsVarACI     = 'asymmetric_cloud_impact'
obsVarY       = 'y'
obsVarH       = 'h'
obsVarDT      = 'datetime'
obsVarGlint   = 'glint'
obsVarImpact  = 'impactHeightRO'
obsVarLT      = 'LocalTime'
obsVarDep     = 'y-h'
obsVarClearSkyDep = 'y-hclr'
obsVarLogDepRatio = 'logDeparture'
obsVarNormDep = 'dσ\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
obsVarNormDep = 'ENI' #error-normalized innovation
obsVarQC      = 'QCflag'
obsVarCI      = 'cloud_impact'
obsVarLogCI   = 'logarithm_cloud_impact'
obsVarMCI     = 'modeled_cloud_impact'
obsVarOCI     = 'observed_cloud_impact'

degree= u'\N{DEGREE SIGN}'

# key = variable name
# value = [unit used, abbreviation]
varDictObs = {
    'airTemperature': ['K', 'T'],
    'bendingAngle': ['%', 'Bnd'],
    obsVarBT: ['K', 'BT'],
    obsVarY: ['K', obsVarY],
    obsVarH: ['K', obsVarH],
    'windEastward': ['m/s', 'U'],
    'windNorthward': ['m/s', 'V'],
    'atmosphericRefractivity': ['%', 'Ref'],
    'specificHumidity': ['kg/kg', 'qv'],
    'stationPressure': ['Pa', 'Ps'],
    'virtualTemperature': ['K', 'Tv'],
    obsVarAlt: ['m', 'alt'],
    obsVarACI: ['K', 'ACI'],
    obsVarCldFracY: [miss_s, 'CFy'],
    obsVarCldFracX: [miss_s, 'CFx'],
    obsVarImpact: ['m', 'impact_height'],
    obsVarLandFrac: [miss_s, 'landfrac'],
    obsVarLat: [degree, 'lat'],
    obsVarLon: [degree, 'lon'],
    obsVarLT: ['hr', obsVarLT],
    'Normalized departure': [miss_s, obsVarNormDep],
    obsVarPrs: ['hPa', 'P'],
    obsVarQC: [miss_s, obsVarQC],
    obsVarCI: ['K', 'CI'],
    obsVarLogCI: ['K', 'Log-CI'],
    obsVarMCI: ['K', 'MCI'],
    obsVarOCI: ['K', 'OCI'],
    obsVarSenZen: [degree, 'zenith'],
    obsVarGlint: [degree, obsVarGlint],
}
#Note, atmosphericRefractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bendingAngle: we plot RMSE of OMB/O and OMA/O; bending angle binVar == obsVarAlt:

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
anIter = str(nOuterIter)

# Fixed IODA ObsGroups for observation-type variables
obsGroup    = 'ObsValue'
metaGroup   = 'MetaData'
typeGroup   = 'ObsType'

# Dynamic MPAS-Workflow ObsGroups for observation-type variables
ommGroup    = 'ommIter'

# GeoVaLs ObsGroups only used in post-processing, not in files
diagGroup   = 'ObsDiag'
geoGroup    = 'GeoVaLs'

# Generic variable names used as baseVar
selfObsValue = 'selfObsValue'
selfObsType = 'selfObsType'
selfOMMValue = 'selfOMMValue'
selfOMBValue = 'selfOMBValue'
selfOMAValue = 'selfOMAValue'
selfHofXValue = 'selfHofXValue'
selfErrorValue = 'selfErrorValue'
selfQCValue = 'selfQCValue'
selfBCValue = 'selfBCValue'
bgHofXValue = 'bgHofXValue'
anHofXValue = 'anHofXValue'
altMeta = 'altMeta'
cldfracMeta = 'cldfracMeta'
datetimeMeta = 'datetimeMeta'
impactMeta = 'impactMeta'
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
ObsGroups[selfObsType] = typeGroup

ObsGroups[selfOMMValue] = ommGroup
ObsGroups[selfOMBValue] = ombgGroup
ObsGroups[selfOMAValue] = omanGroup

ObsGroups[selfHofXValue] = hofxGroup
ObsGroups[bgHofXValue] = hofxGroup+bgIter
ObsGroups[anHofXValue] = hofxGroup+anIter

ObsGroups[selfErrorValue] = errorGroup
ObsGroups[selfQCValue] = qcGroup
ObsGroups[selfBCValue] = bcGroup

for key in ObsGroups.keys():
  ObsVars[key] = vNameStr

# Fixed IODA variable names (MetaData)
ObsVars[altMeta]      = obsVarAlt
ObsVars[cldfracMeta]  = obsVarCldFracY
ObsVars[datetimeMeta] = obsVarDT
ObsVars[impactMeta]   = obsVarImpact
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
  # separate integer suffixes (e.g., brightnessTemperature_*)
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




## Variable names for MPAS-Model
#modVarAlt = 'zgrid' # --> needs to be interpolated to nVertLevels instead of nVertLevelsP1
modVarLev = 'modelLevel'
modVarLat = 'latCell'
modVarLon = 'lonCell'
modVarPrs = 'P3D'
modVarDiagPrs = 'diagnostic_pressure'

levModel = 'levModel'
latModel = 'latModel'
lonModel = 'lonModel'
prsModelDiag = 'prsModelDiag'

AllVars = deepcopy(ObsVars)
AllGroups = deepcopy(ObsGroups)

AllVars[levModel] = modVarLev
AllVars[latModel] = modVarLat
AllVars[lonModel] = modVarLon
AllVars[prsModelDiag] = modVarDiagPrs

AllGroups[levModel] = None
AllGroups[latModel] = None
AllGroups[lonModel] = None
AllGroups[prsModelDiag] = None

modDiagnosticIndependentVars = [prsModelDiag]

kgm3 = 'kg/m\N{SUPERSCRIPT THREE}'

# key = variable name
# value = [unit used, abbreviation]
varDictModel = {
  modVarLev: [miss_s, 'ModLev'],
  modVarLat: [degree, 'lat'],
  modVarLon: [degree, 'lon'],
  modVarPrs: ['Pa', 'P'],
  'q2': ['g/kg', 'Q2m'],
  'qv': ['g/kg', 'Qv'],
  'qv01to30': ['g/kg', 'Qv01to30'],
  'qv01to10': ['g/kg', 'Qv01to10'],
  'qv11to20': ['g/kg', 'Qv11to20'],
  'qv21to30': ['g/kg', 'Qv21to30'],
  'qv31to40': ['g/kg', 'Qv31to40'],
  'qv41to55': ['g/kg', 'Qv41to55'],
  'rho': [kgm3, 'rho'],
  'surface_pressure': ['Pa', 'Ps'],
  't2m': ['C', 'T2m'],
  'temperature': ['C', 'T'],
  'theta': ['K', 'Theta'],
  'u': ['m/s', 'uedge'],
  'u10': ['m/s', 'U10m'],
  'uReconstructZonal': ['m/s', 'U'],
  'uReconstructMeridional': ['m/s', 'V'],
  'v10': ['m/s', 'V10m'],
  'w': ['m/s', 'W'],
  # diagnostic variables at fixed pressure levels
  modVarDiagPrs: ['Pa', 'DiagnosticP'],
  'diagnostic_relhum': ['%', 'RH'],
  'diagnostic_dewpoint': ['K', 'Tdp'],
  'diagnostic_temperature': ['K', 'T'],
  'diagnostic_height': ['m', 'height'],
  'diagnostic_uzonal': ['m/s', 'U'],
  'diagnostic_umeridional': ['m/s', 'V'],
  'diagnostic_w': ['m/s', 'W'],
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
  modVarPrs,
  'qv',
  'rho',
  'temperature',
  'theta',
  'uReconstructMeridional',
  'uReconstructZonal',
  'w',
]

modDiagnosticVarNames = [
  'diagnostic_relhum',
  'diagnostic_dewpoint',
  'diagnostic_temperature',
  'diagnostic_height',
  'diagnostic_uzonal',
  'diagnostic_umeridional',
  'diagnostic_w',
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

modVarNames3d += modDiagnosticVarNames

def modelVarAttributes(var):
    # return short name and units
    dictName, suf = splitIntSuffix(var)
    varAtt = varDictModel.get(dictName,[miss_s,dictName])
    varShort = varAtt[1]+suf
    varUnits = varAtt[0]
    return varShort, varUnits

#misc. constants; TODO: collect into single script
deg2rad = np.pi / float(180.0)
rad2deg = float(180.0) / np.pi


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
      dbVar = re.sub(ommGroup,ombgGroup,dbVar)
    else:
      dbVar = re.sub(ommGroup,omanGroup,dbVar)
  return dbVar
