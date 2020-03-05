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
obsVarCa      = 'symmetric_cloud_impact'
obsVarCldFrac = 'cloud_area_fraction'
obsVarLon     = 'longitude'
obsVarLat     = 'latitude'
obsVarPrs     = 'air_pressure'
obsVarQC      = 'QCflag'
obsVarSatZen  = 'sensor_zenith_angle'

degree= u'\N{DEGREE SIGN}'

# columns: var_name            unit_used   abbr.
varDictObs = {
    'air_temperature':        [ 'K',     'T'       ]
  , 'bending_angle':          [ '%',     'Bnd'     ]
  , 'brightness_temperature': [ 'K',     'BT'      ]
  , 'eastward_wind':          [ 'm/s',   'U'       ]
  , 'northward_wind':         [ 'm/s',   'V'       ]
  , 'refractivity':           [ '%',     'Ref'     ]
  , 'specific_humidity':      [ 'kg/kg', 'qv'      ]
  , 'surface_pressure':       [ 'Pa',    'Ps'      ]
  , 'virtual_temperature':    [ 'K',     'Tv'      ]
  , obsVarAlt:                [ 'm',     'alt'     ]
  , obsVarCa:                 [ 'K',     'Ca'      ]
  , obsVarCldFrac:            [ miss_s,  'cldfrac' ]
  , obsVarLat:                [ degree,  'lat'     ]
  , obsVarLon:                [ degree,  'lon'     ]
  , obsVarPrs:                [ 'hPa',   'P'       ]
  , obsVarQC:                 [ miss_s,  obsVarQC  ]
  , obsVarSatZen:             [ degree,  'zenith'  ]
    }
#Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == obsVarAlt:

# NC observation variable name substitutions
vNameStr = 'varName'
vChanStr = 'varCHAN'

# MPAS-JEDI suffixes for observation-type variables
bakGroup    = 'HofX'
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

#selfBakValue  = vNameStr+'@'+bakGroup
selfQCValue    = vNameStr+'@'+qcGroup
selfErrorValue = vNameStr+'@'+errorGroup

altMeta      = obsVarAlt+'@'+metaGroup
cldfracMeta  = obsVarCldFrac+'@'+metaGroup
lonMeta      = obsVarLon+'@'+metaGroup
latMeta      = obsVarLat+'@'+metaGroup
prsMeta      = obsVarPrs+'@'+metaGroup
satzenMeta   = obsVarSatZen+'@'+metaGroup

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
        obsVarName = ''.join(obsVarName.split("_")[:-1])
    return obsVarName, suf


def varAttributes(var):
    # return short name and units
    dictName, suf = splitIntSuffix(var)
    varVal = varDictObs.get(dictName,[miss_s,dictName])
    if suf != '':
        varShort = varVal[1]+'_'+suf
    else:
        varShort = varVal[1]
    varUnits = varVal[0]
    return varShort, varUnits


bgIter = '0'
def base2dbVar(baseVar,varName,Iter):
    dictName, suf = splitIntSuffix(varName)
    dbVar = re.sub(vNameStr,varName,baseVar)
    dbVar = re.sub(vChanStr,suf,dbVar)
    dbVar = re.sub(qcGroup,qcGroup+Iter,dbVar)
    dbVar = re.sub(errorGroup,errorGroup+Iter,dbVar)
    if Iter == bgIter:
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


#===========================
# obs diagnostic definitions
#===========================

allDiags = ['omb','oma','obs','bak','ana']
nonObsDiags = [diag for diag in allDiags if diag!='obs']


def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
