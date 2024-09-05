import binning_utils as bu
from copy import deepcopy
import datetime as dt
from datetime import datetime, timedelta
import multiprocessing as mp
from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
import os
import sys
import var_utils as vu

ncWriteFormat = 'NETCDF3_64BIT_OFFSET'

fcHours = os.getenv('fcHours', '0')
intervalHours = os.getenv('intervalHours', '6')
fcNums = os.getenv('fcNums', '1')
initDate = os.getenv('start_init', '2018041500')
endDate  = os.getenv('end_init', '2018051400')
diff2exp = os.getenv('diff2exp', 'False')
expDirectory = os.getenv('TOP_DIR','/glade/scratch/$user/pandac/')

GFSANA_DIR = os.getenv('GFSANA_DIR', 'Please link GFSANA_DIR')
expLongNames = os.getenv('expLongNames', 'please set expLongNames')
expNames = os.getenv('expNames','please set expNames')
ncells = os.getenv('ncells', '40962')
#for ci:
EXP_DIR1 = os.getenv('FCDIAG_WORK_DIR1','FC1DIAG DIR OR FC2DIAG DIR FOR CONTROL')
EXP_DIR2 = os.getenv('FCDIAG_WORK_DIR2','FC1DIAG DIR OR FC2DIAG DIR FOR CURRENT/target exp')
exp1Name = os.getenv('exp1Name','name for control expt')
exp2Name = os.getenv('exp2Name','name for current/target expt')
#
aggregatableFileStats = ['RMS','Mean'] #,'STD','MS', 'Min','Max']
allFileStats = aggregatableFileStats
expStats = 'expmgfs'
varNames2d = ['t2m','surface_pressure','q2','u10','v10']
varNames3d = ['theta','temperature','rho',vu.modVarPrs,'uReconstructZonal','uReconstructMeridional','qv','w']
varNames = varNames2d + varNames3d
latBands = ['NXTro','Tro','SXTro']
latBandsBounds = [90.0, 30.0, -30.0, -90.0]

inflationVariables = [
  'temperature',
  'uReconstructZonal',
  'uReconstructMeridional',
  'spechum',
  'surface_pressure',
  'qc',
  'qi',
  'qr',
  'qs',
  'qg',
#non-increment-variables in output stream
  'qv',
]


variableTraits = {
  'temperature': {
    'templateVar': '2D-c-c',
    'units': 'K',
    'long_name': 'temperature',
  },
  'uReconstructZonal': {
    'templateVar': '2D-c-c',
    'units': 'm s^{-1}',
    'long_name': 'Zonal component of reconstructed horizontal velocity at cell centers',
  },
  'uReconstructMeridional': {
    'templateVar': '2D-c-c',
    'units': 'm s^{-1}',
    'long_name': 'Meridional component of reconstructed horizontal velocity at cell centers',
  },
  'spechum': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Specific humidity',
  },
  'surface_pressure': {
    'templateVar': '1D-c',
    'units': 'Pa',
    'long_name': 'Diagnosed surface pressure',
  },
  'qc': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Cloud water mixing ratio',
  },
  'qi': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Ice mixing ratio',
  },
  'qr': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Rain mixing ratio',
  },
  'qs': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Snow mixing ratio',
  },
  'qg': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Graupel mixing ratio',
  },
  'qv': {
    'templateVar': '2D-c-c',
    'units': 'kg kg^{-1}',
    'long_name': 'Water vapor mixing ratio',
  },
#  'pressure_base': '2D-c-c',
#  'pressure_p': '2D-c-c',
#  'theta': '2D-c-c',
#  'rho': '2D-c-c',
##  'u': '2D-e-c',
#  'xice': '1D-c',
#  'snowc': '1D-c',
#  'skintemp': '1D-c',
#  'snowh': '1D-c',
#  'vegfra': '1D-c',
#  'u10': '1D-c',
#  'v10': '1D-c',
#  'lai': '1D-c',
#  'smois': '2D-c-s',
#  'tslb': '2D-c-s',
##  'w': '2D-c-w',
#  're_cloud': '2D-c-c',
#  're_ice': '2D-c-c',
#  're_snow': '2D-c-c',
#  'cldfrac': '2D-c-c',
# TODO: replace template variables with actual dimensions?
#  'temperature': [nCells, nVertLevels],
#  'uReconstructZonal': [nCells, nVertLevels],
#  'uReconstructMeridional': [nCells, nVertLevels],
#  'spechum': [nCells, nVertLevels],
#  'surface_pressure': [nCells],
}
templateVariables = {
  '1D-c': 'surface_pressure',
  '2D-c-c': 'theta',
#  '2D-c-w': 'w',
#  '2D-e-c': 'u',
#  '2D-c-s': 'smois',
}

#setting for 120km res.:
nlevelSurface = 1
nlevels = 55
nlevelsP1 =56

#
fcRange = int(fcHours)/24.
interval = int(intervalHours)/24
SDATE = datetime.strptime(initDate,'%Y%m%d%H')
EDATE = datetime.strptime(endDate,'%Y%m%d%H')

DATEINC = dt.timedelta(days=interval)
expLongNames = expLongNames.split()
expNames = expNames.split()
nExp  = len(expNames)

def getGridFile(date = initDate, gfsAnaDir = GFSANA_DIR, nCells = ncells):
  date = initDate
  print(date)
  d = datetime.strptime(date,'%Y%m%d%H')
  filedate= d.strftime('%Y-%m-%d_%H.%M.%S')

  return gfsAnaDir+'/x1.'+str(nCells)+'.init.'+filedate+'.nc'

def readGrid(date=initDate, gridFile=None):
  if gridFile is None: gridFile = getGridFile(date)
  ncData = Dataset(gridFile, 'r')
  grid = {}
  grid['latitude'] = np.array( ncData.variables['latCell'][:] ) * 180.0 / np.pi
  grid['longitude'] = np.array( ncData.variables['lonCell'][:] ) * 180.0 / np.pi
  grid['area'] = np.array( ncData.variables['areaCell'][:] )
  grid['R'] = ncData.__dict__['sphere_radius']
  ncData.close()

  return grid

def hasVar(varName, ncData):
  return (varName in ncData.variables)

def varDims(varName, ncData):
  return ncData.variables[varName].dimensions

def varAttrs(varName, ncData):
  return ncData[varName].__dict__

def varDatatype(varName, ncData):
  return ncData.variables[varName].datatype

def getPressure(ncData):
  if 'pressure' in ncData.variables:
    pressure = np.array( ncData.variables['pressure'][0,:,:] )
  elif 'pressure_p' in ncData.variables and 'pressure_base' in ncData.variables:
    pressure_p = np.array( ncData.variables['pressure_p'][0,:,:] )
    pressure_base = np.array( ncData.variables['pressure_base'][0,:,:] )
    pressure = pressure_p + pressure_base
  else:
    print("Error: 'pressure', 'pressure_p', or 'pressure_base' variables are missing from the NetCDF data.")
  return pressure

def getTemperature(ncData):
  rgas = 287.0
  cp = 1004.5
  pressure = getPressure(ncData)
  theta = np.array(ncData.variables['theta'][0,:,:])
  tmp1 = (1.e5 / pressure)**(rgas / cp)
  temperature = np.subtract(np.divide(theta, tmp1), 273.15)

  return temperature

def getNCData(ncFile, mode='r'):
  return Dataset(ncFile, mode)

# diagnostic pressures (hPa) for diagnostics
diagnosticPressures = np.array([50., 100., 200., 250., 500., 700., 850., 925.])

def varRead(varName, ncData, nprocs=1):

  if varName in vu.modDiagnosticVarNames:
    nLocs = getTemperature(ncData['state']).shape[0]
    nLevs = len(diagnosticPressures)
    varVals = np.full((nLocs, nLevs), np.NaN)

    # allow for missing values for ensembles (None)
    if ncData['diagnostics'] is not None:
      # read all 1D nc arrays into single 2D numpy array
      for iLev, p in enumerate(diagnosticPressures):
        layerVarName = (varName+'_{:d}hPa'.format(int(p))).replace('diagnostic_','')
        if layerVarName in ncData['diagnostics'].variables:
          varVals[:, iLev] = np.array( ncData['diagnostics'].variables[layerVarName][0,:] )

  else:
    if (varName == 'temperature'):
      varVals = getTemperature(ncData['state'])
    elif (varName == vu.modVarPrs):
      varVals = getPressure(ncData['state'])
    else:
      if varName in varNames3d:
        varVals = np.array( ncData['state'].variables[varName][0,:,:] )
      else:
        varVals = np.array( ncData['state'].variables[varName][0,:] )

      if (varName in ['qv', 'rho', 'q2']):
        varVals = varVals * 1000.

  return varVals

def varWrite(varName, varVals, dst,
             varAttrs,
             dims, datatype):

  #TODO: Only appending to existing file is enabled in varWrite
  #      need to make sure dst has 'a' mode

  #assert os.path.exists(ncFile), 'Only appending to existing file is enabled in varWrite!'

  nDims = len(varVals.shape)

  if not hasVar(varName, dst):
    dst.createVariable(varName, datatype, dims)
    dst[varName].setncatts(deepcopy(varAttrs))

  if 'Time' in dims:
    if nDims == 1:
      dst[varName][0,:] = np.asarray(varVals[:], dtype=datatype)
    elif nDims == 2:
      dst[varName][0,:,:] = np.asarray(varVals[:,:], dtype=datatype)
  else:
    if nDims == 1:
      dst[varName][:] = np.asarray(varVals[:], dtype=datatype)
    elif nDims == 2:
      dst[varName][:,:] = np.asarray(varVals[:,:], dtype=datatype)

def headerOnlyFileFromTemplate(template, outfile, date):
  dst = Dataset(outfile, 'w', format=ncWriteFormat)

  # copy global attributes all at once via dictionary
  dstatts = deepcopy(template.__dict__)

  d = datetime.strptime(date, '%Y%m%d%H')
  confdate = d.strftime('%Y-%m-%d_%H:%M:%S')
  dstatts['config_start_time'] = confdate
  dstatts['config_stop_time'] = confdate
  dst.setncatts(dstatts)

  # copy dimensions
  for name, dimension in template.dimensions.items():
    dst.createDimension(
      name, (len(dimension) if not dimension.isunlimited() else None))

  requiredMetaVars = ['xtime']
  for varname in requiredMetaVars:
    srcvar = template.variables[varname]
    # create variable in destination
    dst.createVariable(
       varname,
       srcvar.datatype,
       srcvar.dimensions)

    # copy variable attributes all at once via dictionary
    dst[varname].setncatts(template[varname].__dict__)

    # finally copy variable data
    if varname == 'xtime':
      xtime_ = template[varname][0][:]
      xtime = nc.stringtoarr(confdate, len(template[varname][0][:]))
      dst[varname][0] = xtime
    else:
      if 'Time' in srcvar.dimensions:
        dst[varname][0,:] = template[varname][0,:]
      else:
        dst[varname][:] = template[varname][:]

  dst.close()


## field functions for a single varName
#class fieldFunction():
#  def evaluate(self, varName, fieldsDB, nprocs=1):
#    '''
#    virtual method
#    '''
#    raise NotImplementedError()


class fieldSpread():
  def __init__(self, stateType: str):
    self.stateType = stateType

  def evaluate(self, varName: str, fieldsDB, nprocs=1):
    exp = varRead(varName, fieldsDB['mean'], nprocs)
    nLocs = exp.size
    ensemble = fieldsDB[self.stateType+'Ens']

    nEns = len(ensemble)
    if nEns > 1:
      exps = np.full((nEns, nLocs), np.NaN)
      for ii, member in enumerate(ensemble):
        exps[ii,:] = varRead(varName, member, nprocs).flatten()

      std = np.nanstd(exps, axis=0, ddof=1).reshape(exp.shape)

      del exps

    else:
      std = np.full(exp.shape, np.NaN)

    return np.asarray(std)


class fieldDiff():
  def __init__(self):
    pass

  def evaluate(self, varName: str, fieldsDB, nprocs=1):
    ref = varRead(varName, fieldsDB['reference'], nprocs)
    exp = varRead(varName, fieldsDB['mean'], nprocs)
    return np.asarray(exp - ref)


class fieldRelativeDiff():
  def __init__(self):
    pass

  def evaluate(self, varName: str, fieldsDB, nprocs=1):
    # only valid for positive-definite reference state
    ref = varRead(varName, fieldsDB['reference'], nprocs)
    exp = varRead(varName, fieldsDB['mean'], nprocs)

    d = np.empty_like(ref)
    valid = bu.greatBound(np.abs(ref), 0.)
    d[valid] = (exp[valid] - ref[valid]) / ref[valid]

    return np.asarray(np.multiply(d, 100.0))


class fieldLogRatio():
  def __init__(self):
    pass

  def evaluate(self, varName: str, fieldsDB, nprocs=1):
    # only valid for positive-definite reference and mean
    ref = varRead(varName, fieldsDB['reference'], nprocs)
    exp = varRead(varName, fieldsDB['mean'], nprocs)

    d = np.empty_like(ref)
    valid = bu.greatBound(ref, 0.)
    valid = np.logical_and(bu.greatBound(exp, 0.), valid)

    d[valid] = np.log10(np.divide(exp[valid], ref[valid]))

    return np.asarray(d)


diagnosticFunctions = {
  'mmgfsan': fieldDiff(),
  'rltv_mmgfsan': fieldRelativeDiff(),
  'log_mogfsan': fieldLogRatio(),
  'sigmaxb': fieldSpread('bg'),
  'sigmaxa': fieldSpread('an'),
  'sigmaxinf': fieldSpread('inf'),
  'sigmaxf': fieldSpread('fc'),
}

variableSpecificDiagnosticConfigs = {
  # default is 'mmgfsan' below
  vu.modVarPrs: ['mmgfsan', 'log_mogfsan'],
  'q2': ['mmgfsan', 'log_mogfsan'],
  'qv': ['mmgfsan', 'log_mogfsan'],
  'qv01to30': ['mmgfsan', 'log_mogfsan'],
  'qv01to10': ['mmgfsan', 'log_mogfsan'],
  'qv11to20': ['mmgfsan', 'log_mogfsan'],
  'qv21to30': ['mmgfsan', 'log_mogfsan'],
  'qv31to40': ['mmgfsan', 'log_mogfsan'],
  'qv41to55': ['mmgfsan', 'log_mogfsan'],
}

def variableSpecificDiagnostics(varName: str, nEns: int):
  diags = variableSpecificDiagnosticConfigs.get(varName, ['mmgfsan'])
  if nEns > 1:
    diags.append('sigmaxb')
    diags.append('sigmaxa')
    diags.append('sigmaxinf')
    #diags.append('sigmaxf')

  return diags

aggregatedVariableConfig = {
  'qv01to30': {
    'model variable': 'qv',
    'max level': 30,
  },
  'qv01to10': {
    'model variable': 'qv',
    'min level': 1,
    'max level': 10,
  },
  'qv11to20': {
    'model variable': 'qv',
    'min level': 11,
    'max level': 20,
  },
  'qv21to30': {
    'model variable': 'qv',
    'min level': 21,
    'max level': 30,
  },
  'qv31to40': {
    'model variable': 'qv',
    'min level': 31,
    'max level': 40,
  },
  'qv41to55': {
    'model variable': 'qv',
    'min level': 41,
    'max level': 55,
  },
}

def aggVariableConfig(varName: str):
  return aggregatedVariableConfig.get(
    varName,
    {'model variable': varName}
  )

def aggModelVariable(varName: str):
  return aggVariableConfig(varName)['model variable']

def aggMinLevel(varName: str):
  return np.max([aggVariableConfig(varName).get('min level', 1), 1])

def aggMaxLevel(varName: str, nLevels: int):
  return np.min([aggVariableConfig(varName).get('max level', nLevels), nLevels])

def main():
  print ('This is not a runnable program.')
  os._exit(0)

if __name__ == '__main__': main()

