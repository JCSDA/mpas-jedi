import datetime as dt
import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import deepcopy
from datetime import datetime, timedelta

fcHours = os.getenv('fcHours', '0')
intervalHours = os.getenv('intervalHours', '6')
GFSANA_DIR = os.getenv('GFSANA_DIR', '/gpfs/fs1/scratch/jban/pandac/120km_GFSANA')
fcNums = os.getenv('fcNums', '1')

aggregatableFileStats = ['Mean','RMS','STD','MS'] # [ 'Min','Max']
allFileStats = aggregatableFileStats

expStats = 'expmgfs'
varNames = ['pressure','theta','uReconstructZonal','uReconstructMeridional','qv']

latBands = ['NXTro','Tro','SXTro']
latBandsBounds = [90.0, 30.0, -30.0, -90.0]

fcRange = int(fcHours)/24.
interval = int(intervalHours)/24.

def write_diag_stats():
    if os.path.exists(expStats+'.nc'):
        os.remove(expStats+'.nc')    
    path = os.getcwd()
    date = path.split('/')[-2]
    initDate = datetime.strptime(date,"%Y%m%d%H")
    lats, lons = readGrid()

    for varName in varNames:
        for latBand in range(0, len(latBands)):
            tmp = []
            meanncs = []
            rmsncs  = []
            msncs  = []

            for fcTDelta in np.arange(0, fcRange+interval,interval):
                fcDate = initDate + timedelta(days=fcTDelta) 
                fileDate= fcDate.strftime("%Y-%m-%d_%H.%M.%S")
                ncFile1 = GFSANA_DIR+'/x1.40962.init.'+fileDate+'.nc'
                if (fcTDelta ==0):
                    ncFile2 = '../analysis.'+fileDate+'.nc'
                else:
                    ncFile2 = '../restart.'+fileDate+'.nc'

                if (varName == 'pressure'):
                    if (fcTDelta == 0):
                        pressure_model = varRead(varName,ncFile2)
                    else:
                        pressure_model =  varRead('pressure_base',ncFile2) + varRead('pressure_p',ncFile2)
                    pressure_gfs = varRead('pressure_base',ncFile1) + varRead('pressure_p',ncFile1)
                    tmp = pressure_model - pressure_gfs
                else:
                    tmp = varDiff(varName,ncFile1,ncFile2) 

                #bin for regions
                tmpbin = []
                tmpbin = deepcopy(tmp)
                tmpbin[np.logical_or(lats < latBandsBounds [latBand+1], lats > latBandsBounds [latBand])] = np.NaN
                #save every level stat
                newfile = write_stats(tmpbin,varName,latBands[latBand],str(fcTDelta))

                meannc = np.nanmean(tmpbin.flatten(),axis=0)
                rmsnc = np.sqrt(np.nanmean(tmpbin.flatten()**2,axis=0))
                msnc = np.nanmean(tmpbin.flatten()**2,axis=0)

                meanncs = np.append(meanncs,meannc)
                rmsncs = np.append(rmsncs,rmsnc)
                msncs = np.append(msncs,msnc)                
            #save all levels stat
            newfile = write_stats_levall(meanncs,varName,latBands[latBand],'Mean')
            newfile = write_stats_levall(rmsncs,varName,latBands[latBand],'RMS')
            newfile = write_stats_levall(msncs,varName,latBands[latBand],'MS')
#TODO: move the following functions to plot_utils.py or basic_plot_functions.py
def readGrid():
    path = os.getcwd()
    date = path.split('/')[-2]
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
    varsVals = np.array( nc_fid2.variables[varNames][0,:,:] ) - np.array( nc_fid1.variables[varNames][0,:,:] )
    if (varNames == 'qv'):
        varsVals = varsVals * 1000.
    return varsVals

def varRead(varNames,ncFile1):
    nc_fid1 = Dataset(ncFile1, "r", format="NETCDF4")
    varsVals = np.array( nc_fid1.variables[varNames][0,:,:] )
    if (varNames == 'qv'):
        varsVals = varsVals * 1000.
    return varsVals

def write_stats(array_f,varNames,band,fcTDelta):
    STATS = {}
    STATS['Mean'] = np.nanmean(array_f,axis=0)
    STATS['MS']  = np.nanmean(array_f**2,axis=0)
    STATS['RMS']  = np.sqrt(np.nanmean(array_f**2,axis=0))
    STATS['STD']  = np.nanstd(array_f,axis=0)
    STATS['Min']  = np.nanmin(array_f,axis=0)
    STATS['Max']  = np.nanmax(array_f,axis=0)

    if os.path.exists(expStats+'.nc'):
        w_nc_fid = Dataset(expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics" 
        w_nc_fid.createDimension('level', 55)
        w_nc_fid.createDimension('fcnums', int(fcNums))
    for statName in allFileStats:
        bias2exp = w_nc_fid.createVariable(expStats+"_day"+fcTDelta+"_"+band+"_"+varNames+"_"+statName,'f4', "level")
        w_nc_fid.variables[expStats+"_day"+fcTDelta+'_'+band+'_'+varNames+'_'+statName][:] = STATS[statName]

    w_nc_fid.close()

def write_stats_levall(array_f,varNames,band,statName):
    if os.path.exists(expStats+'.nc'):
        w_nc_fid = Dataset(expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics"

    bias2exp = w_nc_fid.createVariable(expStats+"_"+band+"_"+varNames+"_"+statName,'f4', "fcnums")
    w_nc_fid.variables[expStats+'_'+band+'_'+varNames+'_'+statName][:] = array_f

    w_nc_fid.close()

def main():
    write_diag_stats()

if __name__ == '__main__': main()
