import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import deepcopy
import datetime as dt
import plot_utils as pu
import var_utils as vu
import modelsp_utils as mu

def readdata():
   for metrics in mu.allFileStats:
       for varName in mu.varNames2d:
           for latBand in mu.latBands:
             arraylist = []
             for iexp, expName in enumerate(mu.expNames):
               xlabeltime  = []
               alldata = []
               alldiffdata = []
               alldiffdata_rmsdiv = [] 
               nc_file =  mu.expDirectory+'/'+mu.expLongNames[iexp]+'/FC2DIAG/expmgfs.nc'
               nc_fid = Dataset(nc_file, "r", format="NETCDF4")
               for fcTDelta in np.arange(0,mu.fcRange+mu.interval,mu.interval):
                   varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                   # data1: exp1-GFSANA 
                   data = np.array( nc_fid.variables[''.join(varNamesList)][:]) 
                   alldata = np.append(alldata, data)
                   xlabeltime = np.append(xlabeltime,fcTDelta)
               varNamesListUse = 'expmgfs_fc_'+ latBand +'_'+ varName + '_' + metrics
               if (iexp == 0):
                   arraylist = [alldata]
               else:
                   arraylist= arraylist + [alldata]
             plotTimeSerial(arraylist,xlabeltime,varNamesListUse) 
#TODO: move this part to basic_plot_functions.py
def plotTimeSerial(linesVals,xlabeltime,VarName):

    fig,ax1 = plt.subplots(1,sharex=True)
    plt.grid(True)
    xarray = range(len(xlabeltime))
    major_ticks = np.arange(0, len(xlabeltime), 1)

    if (VarName == 'qv' or VarName == 'rho' or VarName == 'q2'):
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    nx = mu.nExp
    for iexp in range(0,mu.nExp):
       ax1.plot(xarray,linesVals[iexp],pu.plotSpecs[iexp],markersize=5)

    ax1.set_xticks(major_ticks[::4])
    # upper right
    ax1.legend(mu.expNames, loc='upper left',fontsize=12,frameon=False)

    FCDay = VarName.split("_")[1]
    if (FCDay == 'day0.0'):
        ax1.set_xlabel('Analysis Time',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    elif (FCDay == 'day0.25'):
        ax1.set_xlabel( '6h Forecast',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    else:
        ax1.set_xlabel( 'Lead Time (day)',fontsize=15)
        ax1.set_xticks(major_ticks[::1])
        ax1.set_xticks(xarray[::1])
        ax1.set_xticklabels(xlabeltime[::1],rotation=0)

    ax1.grid()
    #ax1.set_xticklabels(xlabeltime[::4])
    plt.xticks(rotation=90)
    plt.grid(True)
    region = VarName.split("_")[2]
    var    = '_'.join(VarName.split("_")[3:][:-1])  # surface_pressure
    stats  = ''.join(VarName.split("_")[-1:])
    ax1.set_ylabel(stats,fontsize=15)
    plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
    plt.savefig(VarName+'.png',dpi=300,bbox_inches='tight')
    plt.close()
def main():
    readdata()

if __name__ == '__main__': main()
