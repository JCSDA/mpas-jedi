import datetime as dt
import glob
import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from copy import deepcopy
from datetime import datetime, timedelta
import modelsp_utils as mu

def write_diag_stats():
    if os.path.exists(mu.expStats+'.nc'):
        os.remove(mu.expStats+'.nc')
    path = os.getcwd()
    fcdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    allfiledate = []
    grid = mu.readGrid()
    lats = grid['latitude']
    lons = grid['longitude']

    for varName in mu.varNames:
        for latBand in range(0, len(mu.latBands)):
            for fcTDelta in np.arange(0,mu.fcRange+mu.interval,mu.interval):
                alltmp = [[]]
                tmp = []
                for fcdir in range(0,len(fcdirs)):
                    d = datetime.strptime(fcdirs[fcdir], "%Y%m%d%H") + timedelta(days=fcTDelta)
                    fileDate= d.strftime("%Y-%m-%d_%H.%M.%S")
 
                    ncFile1 = mu.GFSANA_DIR+'/x1.'+mu.ncells+'.init.'+fileDate+'.nc'
                    ncFile2 = './'+fcdirs[fcdir]+'/restart.'+fileDate+'.nc'
                    tmp = mu.varDiff(varName,ncFile1,ncFile2)

                    tmpbin = deepcopy(tmp)
                    tmpbin[np.logical_or(lats < mu.latBandsBounds [latBand+1], lats > mu.latBandsBounds [latBand])] = np.NaN

                    if (fcdir == 0):
                        alltmp = np.append(alltmp, tmpbin)
                        if varName in mu.varNames3d:
                            if (varName == 'w'):
                                alltmp = alltmp.reshape(int(mu.ncells),mu.nlevelsP1)
                            else:
                                alltmp = alltmp.reshape(int(mu.ncells),mu.nlevels)
                        else:
                            alltmp = alltmp.reshape(int(mu.ncells),mu.nlevelSurface)

                    if (fcdir != 0):
                        if varName in mu.varNames3d:
                            alltmp = np.append(alltmp, tmpbin, axis=0)
                        else:
                            alltmp = np.append(alltmp, tmpbin)
                    if (fcdir == len(fcdirs)-1):
                        newfile = write_stats(alltmp,varName,mu.latBands[latBand],str(fcTDelta))
#TODO: move the following functions to plot_utils.py or basic_plot_functions.py
def write_stats(array_f,varNames,band,fcTDelta):
    STATS = {}
    STATS['Mean'] = np.nanmean(array_f,axis=0)
    STATS['RMS']  = np.sqrt(np.nanmean(array_f**2,axis=0))
    STATS['STD']  = np.nanstd(array_f,axis=0)
    STATS['Min']  = np.nanmin(array_f,axis=0)
    STATS['Max']  = np.nanmax(array_f,axis=0)

    if os.path.exists(mu.expStats+'.nc'):
        w_nc_fid = Dataset(mu.expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(mu.expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics" 
        w_nc_fid.createDimension('level', mu.nlevels)
        w_nc_fid.createDimension('levelP1',mu.nlevelsP1)
        w_nc_fid.createDimension('levelSurface', mu.nlevelSurface)
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

def main():
    write_diag_stats()

if __name__ == '__main__': main()
