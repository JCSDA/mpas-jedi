from netCDF4 import Dataset
import numpy as np
import os
import sys
import datetime as dt
from datetime import datetime, timedelta

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
varNames3d = ['theta','temperature','rho','pressure','uReconstructZonal','uReconstructMeridional','qv','w']
varNames = varNames2d + varNames3d
latBands = ['NXTro','Tro','SXTro']
latBandsBounds = [90.0, 30.0, -30.0, -90.0]

#setting for 120km res.:
nlevelSurface = 1
nlevels = 55
nlevelsP1 =56
ncells  = 40962

#
fcRange = int(fcHours)/24.
interval = int(intervalHours)/24
SDATE = datetime.strptime(initDate,'%Y%m%d%H')
EDATE = datetime.strptime(endDate,'%Y%m%d%H')

DATEINC = dt.timedelta(days=interval)
expLongNames = expLongNames.split()
expNames = expNames.split()
nExp  = len(expNames)

def getGridFile(date = initDate, gfsAnaDir = GFSANA_DIR, nCells = 40962):
  date = initDate
  print(date)
  d = datetime.strptime(date,'%Y%m%d%H')
  filedate= d.strftime('%Y-%m-%d_%H.%M.%S')

  return gfsAnaDir+'/x1.'+str(nCells)+'.init.'+filedate+'.nc'

def readGrid(date=initDate, gridFile=None):
  if gridFile is None: gridFile = getGridFile(date)
  ncData = Dataset(gridFile, 'r')
  lats = np.array( ncData.variables['latCell'][:] ) * 180.0 / np.pi
  lons = np.array( ncData.variables['lonCell'][:] ) * 180.0 / np.pi
  ncData.close()

  return (lats,lons)

def getPressure(ncData):
  pressure_p = np.array( ncData.variables['pressure_p'][0,:,:] )
  pressure_base = np.array( ncData.variables['pressure_base'][0,:,:] )
  pressure = pressure_p + pressure_base
  return pressure

def getTemperature(ncFile):
  ncData = Dataset(ncFile, 'r')
  rgas = 287.0
  cp = 1004.5
  pressure = getPressure(ncData)
  theta = np.array(ncData.variables['theta'][0,:,:])
  tmp1 = (1.e5 / pressure)**(rgas / cp)
  temperature = np.subtract(np.divide(theta, tmp1), 273.15)

  ncData.close()

  return temperature

def varRead(varName, ncFile):
  if (varName == 'temperature'):
    varVals = getTemperature(ncFile)
  else:
    ncData = Dataset(ncFile, 'r')
    if (varName == 'pressure'):
      varVals = getPressure(ncData)
    elif varName in varNames3d:
      varVals = np.array( ncData.variables[varName][0,:,:] )
    else:
      varVals = np.array( ncData.variables[varName][0,:] )

    if (varName == 'qv' or varName == 'rho' or varName == 'q2'):
      varVals = varVals * 1000.

    ncData.close()

  return varVals

def varDiff(varName, ncFile1, ncFile2):
  var1 = varRead(varName, ncFile1)
  var2 = varRead(varName, ncFile2)
  return var2 - var1

def main():
  print ('This is not a runnable program.')
  os._exit(0)

if __name__ == '__main__': main()

