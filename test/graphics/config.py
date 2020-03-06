import binning_utils as bu
from copy import deepcopy
import os
import var_utils as vu


#=========================================================
# Sub-selections of binVarConfigs for specific DiagSpaces
#=========================================================

nullBinVars = {vu.miss_s: []}

## Generic binVarConfigs that apply to all observation categories
obsBinVars = {
    vu.obsVarQC:  [bu.defaultBinMethod,'bad'],
    vu.obsVarLat: [bu.defaultBinMethod,'NAMED'],
    vu.obsVarLT: [bu.defaultBinMethod],
    vu.obsVarNormErr: [bu.defaultBinMethod],
    'ObsRegion':  ['CONUS'],
}


## binVarConfigs for surface obs
surfBinVars = deepcopy(obsBinVars)

##########################################################
## binVarConfigs for profile obs w/ pressure vertical bins
##########################################################
profPressBinVars = deepcopy(obsBinVars)
profPressBinVars[vu.obsVarPrs] = [bu.defaultBinMethod,bu.PjetMethod]
profPressBinVars[vu.obsVarLat].append(bu.PjetMethod)

# 2D pressure bins with named latitude-band methods
for iband, latBand in enumerate(bu.namedLatBands['values']):
    profPressBinVars[vu.obsVarPrs].append(latBand)


##########################################################
## binVarConfigs for profile obs w/ altitude vertical bins
##########################################################
profAltBinVars = deepcopy(obsBinVars)
profAltBinVars[vu.obsVarAlt] = [bu.defaultBinMethod,bu.altjetMethod]
profAltBinVars[vu.obsVarLat].append(bu.altjetMethod)

# 2D altitude bins with named latitude-band methods
for iband, latBand in enumerate(bu.namedLatBands['values']):
    profAltBinVars[vu.obsVarAlt].append(latBand)


#################################
## binVarConfigs for radiance obs
#################################
radianceBinVars = deepcopy(obsBinVars)
radianceBinVars[vu.obsVarSatZen] = [bu.defaultBinMethod]


#################################
## binVarConfigs for GOES-ABI obs
#################################
abiBinVars = deepcopy(radianceBinVars)
abiBinVars[vu.obsVarCldFrac] = [bu.defaultBinMethod]

# Binning variables with clr-/cld-sky methods
clrcldVars = [
   vu.obsVarLat,
   vu.obsVarLT,
   vu.obsVarCldFrac,
   vu.obsVarSatZen
]
for var in clrcldVars:
    abiBinVars[var].append(bu.clrskyMethod)
    abiBinVars[var].append(bu.cldskyMethod)

# symmetric cloud impact (expensive)
selectSCIMethods = [
    bu.OkamotoMethod,
    bu.ScaleOkamotoMethod,
#    bu.ModHarnischMethod,
#    bu.ScaleModHarnischMethod,
]
abiBinVars[vu.obsVarSCI] = []
for method in selectSCIMethods:
    abiBinVars[vu.obsVarSCI].append(method)
    abiBinVars[vu.obsVarNormErr].append(method)


#########################################
# binVarConfigs for model space variables
#########################################
#modelBinVars = { 'ModelLatBand':  ['NAMED',bu.defaultBinMethod]
#               , 'ModelBox':      ['CONUS']
#               , 'ModelAltitude': [bu.defaultBinMethod]
#               , 'ModelPressure': [bu.defaultBinMethod]
#               }

#=======================
# DiagSpace definitions
# e.g. IODA ObsSpace
#      MPAS ModelSpace
#=======================

nullDiagSpaceInfo = {'DiagSpaceGrp': vu.miss_s,    'process': False, 'binVarConfigs': nullBinVars }

profile_s  = 'profile'
radiance_s = 'radiance'
model_s    = 'MPAS'

# columns: DiagSpace name (YAML)    DiagSpaceGrp              process?                  binVarConfigs
DiagSpaceConfig = {
    'sondes':                {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profPressBinVars }
  , 'aircraft':              {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profPressBinVars }
  , 'satwind':               {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profPressBinVars }
  , 'gnssroref':             {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profAltBinVars   }
  , 'gnssrobndropp1d':       {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profAltBinVars   }
  , 'gnssro':                {'DiagSpaceGrp': profile_s,  'process': True, 'binVarConfigs': profAltBinVars   }
  , 'abi_g16':               {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': abiBinVars,
                              'channels': [8,9,10,11,13,14,15,16] }
  , 'airs_aqua':             {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': [1,6,7] }
  , 'amsua_n15':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [5,6,7,8,9] }
  , 'amsua_n18':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [5,6,7,8,9] }
  , 'amsua_n19':             {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [5,6,7,9] }
  , 'amsua_metop-a':         {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [5,6,9] }
  , 'amsua_metop-b':         {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [] }
  , 'amsua_aqua':            {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [8,9] }
  , 'amsua_n19--hydro':      {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [1,2,3,15] }
  , 'amsua_n19--nohydro':    {'DiagSpaceGrp': radiance_s, 'process': True, 'binVarConfigs': radianceBinVars,
                              'channels': [4,5,6,7,9,10,11,12,13,14] }
  , 'cris-fsr_npp':          {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': [24,26,28,32,37,39] }
  , 'hirs4_metop-a':         {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,16) }
  , 'iasi_metop-a':          {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': [16,29,32,35,38,41,44] }
  , 'mhs_n19':               {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,6) }
  , 'seviri_m08':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': [5] }
  , 'sndrd1_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,16) }
  , 'sndrd2_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,16) }
  , 'sndrd3_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,16) }
  , 'sndrd4_g15':            {'DiagSpaceGrp': radiance_s, 'process': False, 'binVarConfigs': radianceBinVars,
                              'channels': range(1,16) }
#  , 'mpas-gfs':              {'DiagSpaceGrp': model_s,    'process': False, 'binVarConfigs': modelBinVars    }
    }

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
