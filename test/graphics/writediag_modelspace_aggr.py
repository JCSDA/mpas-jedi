import datetime as dt
import glob
import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import deepcopy
from datetime import datetime, timedelta

aggregatableFileStats = ['Mean','RMS','STD','Min','Max']
allFileStats = aggregatableFileStats

expStats = 'expmgfs'
varNames = ['pressure','theta','uReconstructZonal','uReconstructMeridional','qv']

#latitude, north to south
#'NAMED': N extrtropics, tropics, S extrtropics
latBands = ['NXTro','Tro','SXTro']
latBandsBounds = [90.0, 30.0, -30.0, -90.0]

fcRange = 10 # 10day fc

#setting for 120km res.:
nlevels = 55
ncells  = 40962

EXP_DIR = '/gpfs/fs1/scratch/jban/pandac/'
EXP_NAME1 = 'gfsanal_x1.'+str(ncells)
def write_diag_stats():
    if os.path.exists(expStats+'.nc'):
        os.remove(expStats+'.nc')    
    path = os.getcwd()
    fcdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    allfiledate = []
    lats, lons = readGrid()
    for varName in varNames:
        for latBand in range(0, len(latBands)): 
            for fcTDelta in range(0, fcRange+1): 
                alltmp = [[]]
                tmp = []
                for fcdir in range(0,len(fcdirs)):
                    d = datetime.strptime(fcdirs[fcdir], "%Y%m%d%H") + timedelta(days=fcTDelta)
                    fileDate= d.strftime("%Y-%m-%d_%H.%M.%S")
 
                    ncFile1 = EXP_DIR+EXP_NAME1+'/x1.'+str(ncells)+'.init.'+fileDate+'.nc'
                    ncFile2 = './'+fcdirs[fcdir]+'/restart.'+fileDate+'.nc'

                    if (varName == 'pressure'):
                        pressure_model =  varRead('pressure_base',ncFile2) + varRead('pressure_p',ncFile2)
                        pressure_gfs = varRead('pressure_base',ncFile1) + varRead('pressure_p',ncFile1)
                        tmp = pressure_model - pressure_gfs
                    else:
                        tmp = varDiff(varName,ncFile1,ncFile2)

                    tmpbin = deepcopy(tmp)
                    tmpbin[np.logical_or(lats < latBandsBounds [latBand+1], lats > latBandsBounds [latBand])] = np.NaN

                    if (fcdir == 0):
                        alltmp = np.append(alltmp, tmpbin)
                        alltmp = alltmp.reshape(ncells,nlevels)
                    else:
                        alltmp = np.append(alltmp, tmpbin, axis=0)
                    if (fcdir == len(fcdirs)-1):
                        newfile = write_stats(alltmp,varName,latBands[latBand],str(fcTDelta))
#TODO: move the following functions to plot_utils.py or basic_plot_functions.py
def readGrid():
    path = os.getcwd()
    date = '2018041500'
    d = datetime.strptime(date,"%Y%m%d%H")
    filedate= d.strftime("%Y-%m-%d_%H.%M.%S")
    diagdir    = '../'
    ncGFS = EXP_DIR+EXP_NAME1+'/x1.'+str(ncells)+'.init.'+filedate+'.nc'
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
    STATS['RMS']  = np.sqrt(np.nanmean(array_f**2,axis=0))
    STATS['STD']  = np.nanstd(array_f,axis=0)
    STATS['Min']  = np.nanmin(array_f,axis=0)
    STATS['Max']  = np.nanmax(array_f,axis=0)

    if os.path.exists(expStats+'.nc'):
        w_nc_fid = Dataset(expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics" 
        w_nc_fid.createDimension('level', nlevels)

    for statName in allFileStats:
        bias2exp = w_nc_fid.createVariable(expStats+"_day"+fcTDelta+"_"+band+"_"+varNames+"_"+statName,'f4', "level")
        w_nc_fid.variables[expStats+"_day"+fcTDelta+'_'+band+'_'+varNames+'_'+statName][:] = STATS[statName]

    w_nc_fid.close()

def main():
    write_diag_stats()

if __name__ == '__main__': main()
