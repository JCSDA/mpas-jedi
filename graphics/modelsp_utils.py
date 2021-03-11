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
SDATE = datetime.strptime(initDate,"%Y%m%d%H")
EDATE = datetime.strptime(endDate,"%Y%m%d%H")

DATEINC = dt.timedelta(days=interval)
expLongNames = expLongNames.split()
expNames = expNames.split()
nExp  = len(expNames)

def readGrid():
    date = initDate
    print(date)
    d = datetime.strptime(date,"%Y%m%d%H")
    filedate= d.strftime("%Y-%m-%d_%H.%M.%S")
    diagdir    = '../'
    ncGFS = GFSANA_DIR+'/x1.40962.init.'+filedate+'.nc'
    nc_fid2 = Dataset(ncGFS, "r", format="NETCDF4")
    lats = np.array( nc_fid2.variables['latCell'][:] ) * 180.0 / np.pi
    lons = np.array( nc_fid2.variables['lonCell'][:] ) * 180.0 / np.pi
    return (lats,lons)

def varDiff(varNames,ncFile1,ncFile2):
    nc_fid1 = Dataset(ncFile1, "r", format="NETCDF4")
    nc_fid2 = Dataset(ncFile2, "r", format="NETCDF4")
    if varNames in varNames3d:
        varsVals = np.array( nc_fid2.variables[varNames][0,:,:] ) - np.array( nc_fid1.variables[varNames][0,:,:] )
    else:
        varsVals = np.array( nc_fid2.variables[varNames][0,:] ) - np.array( nc_fid1.variables[varNames][0,:] )
    if (varNames == 'qv' or varNames == 'rho' or varNames == 'q2'):
        varsVals = varsVals * 1000.
    return varsVals

def varRead(varNames,ncFile1):
    nc_fid1 = Dataset(ncFile1, "r", format="NETCDF4")
    if varNames in varNames3d:
        varsVals = np.array( nc_fid1.variables[varNames][0,:,:] )
    else:
        varsVals = np.array( nc_fid1.variables[varNames][0,:] )

    if (varNames == 'qv' or varNames == 'rho' or varNames == 'q2'):
        varsVals = varsVals * 1000.
    return varsVals

def getTemperature(varNames,ncFile1):
    nc_fid1 = Dataset(ncFile1, "r", format="NETCDF4")
    rgas = 287.0
    cp = 1004.5
    pressure_p = np.array( nc_fid1.variables['pressure_p'][0,:,:] )
    pressure_base = np.array( nc_fid1.variables['pressure_base'][0,:,:] )
    pressure = np.divide (pressure_p + pressure_base, 100.0)
    theta = np.array( nc_fid1.variables['theta'][0,:,:] )
    tmp1 = (1000/pressure)**(rgas/cp)
    temperature = np.divide(theta, tmp1)
    varsVals = temperature - 273.15
    return varsVals

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()

