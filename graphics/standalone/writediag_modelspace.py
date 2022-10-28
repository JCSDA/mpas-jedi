import datetime as dt
import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from copy import deepcopy
from datetime import datetime, timedelta
import modelsp_utils as mu

def write_diag_stats():
    if os.path.exists(mu.expStats+'.nc'):
        os.remove(mu.expStats+'.nc')
    path = os.getcwd()
    date = path.split('/')[-2]
    initDate = datetime.strptime(date,"%Y%m%d%H")
    grid = mu.readGrid()
    lats = grid['latitude']
    lons = grid['longitude']

    for varName in mu.varNames:
        for latBand in range(0, len(mu.latBands)):
            tmp = []
            meanncs = []
            rmsncs  = []
            msncs  = []

            for fcTDelta in np.arange(0, mu.fcRange+mu.interval,mu.interval):
                fcDate = initDate + timedelta(days=fcTDelta) 
                fileDate= fcDate.strftime("%Y-%m-%d_%H.%M.%S")
                ncFile1 = mu.GFSANA_DIR+'/x1.'+mu.ncells+'.init.'+fileDate+'.nc'
                ncFile2 = '../restart.'+fileDate+'.nc'
                tmp = mu.varDiff(varName,ncFile1,ncFile2)

                #bin for regions
                tmpbin = []
                tmpbin = deepcopy(tmp)
                tmpbin[np.logical_or(lats < mu.latBandsBounds [latBand+1], lats > mu.latBandsBounds [latBand])] = np.NaN
                #save every level stat
                newfile = write_stats(tmpbin,varName,mu.latBands[latBand],str(fcTDelta))

                meannc = np.nanmean(tmpbin.flatten(),axis=0)
                rmsnc = np.sqrt(np.nanmean(tmpbin.flatten()**2,axis=0))
                msnc = np.nanmean(tmpbin.flatten()**2,axis=0)

                meanncs = np.append(meanncs,meannc)
                rmsncs = np.append(rmsncs,rmsnc)
                msncs = np.append(msncs,msnc)                
            #save all levels stat
            newfile = write_stats_levall(meanncs,varName,mu.latBands[latBand],'Mean')
            newfile = write_stats_levall(rmsncs,varName,mu.latBands[latBand],'RMS')
            newfile = write_stats_levall(msncs,varName,mu.latBands[latBand],'MS')
#TODO: move the following functions to plot_utils.py or basic_plot_functions.py
def write_stats(array_f,varNames,band,fcTDelta):
    STATS = {}
    STATS['Mean'] = np.nanmean(array_f,axis=0)
    STATS['MS']  = np.nanmean(array_f**2,axis=0)
    STATS['RMS']  = np.sqrt(np.nanmean(array_f**2,axis=0))
    STATS['STD']  = np.nanstd(array_f,axis=0)
    STATS['Min']  = np.nanmin(array_f,axis=0)
    STATS['Max']  = np.nanmax(array_f,axis=0)
    if os.path.exists(mu.expStats+'.nc'):
        w_nc_fid = Dataset(mu.expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(mu.expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics" 
        w_nc_fid.createDimension('level', 55)
        w_nc_fid.createDimension('fcnums', int(mu.fcNums))
        w_nc_fid.createDimension('levelP1', 56)
        w_nc_fid.createDimension('levelSurface', 1)
    for statName in mu.allFileStats:
        if varNames in mu.varNames3d:
            if (varNames == 'w'):
                bias2exp = w_nc_fid.createVariable(mu.expStats+"_day"+fcTDelta+"_"+band+"_"+varNames+"_"+statName,'f4', "levelP1")
            else:
                bias2exp = w_nc_fid.createVariable(mu.expStats+"_day"+fcTDelta+"_"+band+"_"+varNames+"_"+statName,'f4', "level")
        else:
            bias2exp = w_nc_fid.createVariable(mu.expStats+"_day"+fcTDelta+"_"+band+"_"+varNames+"_"+statName,'f4', "levelSurface")
        w_nc_fid.variables[mu.expStats+"_day"+fcTDelta+'_'+band+'_'+varNames+'_'+statName][:] = STATS[statName]

    w_nc_fid.close()

def write_stats_levall(array_f,varNames,band,statName):
    if os.path.exists(mu.expStats+'.nc'):
        w_nc_fid = Dataset(mu.expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(mu.expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics"

    bias2exp = w_nc_fid.createVariable(mu.expStats+"_"+band+"_"+varNames+"_"+statName,'f4', "fcnums")
    w_nc_fid.variables[mu.expStats+'_'+band+'_'+varNames+'_'+statName][:] = array_f

    w_nc_fid.close()

def main():
    write_diag_stats()

if __name__ == '__main__': main()
