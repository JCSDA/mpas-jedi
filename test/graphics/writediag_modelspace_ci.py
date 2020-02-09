import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import deepcopy
import datetime as dt
from datetime import datetime, timedelta
import numpy.random as npr

EXP_DIR1 = os.getenv('FCDIAG_WORK_DIR1','FC1DIAG DIR OR FC2DIAG DIR FOR EXP1')
EXP_DIR2 = os.getenv('FCDIAG_WORK_DIR2','FC1DIAG DIR OR FC2DIAG DIR FOR CONTROL')

initDate = os.getenv('start_init', '2018041500')
endDate  = os.getenv('end_init', '2018051400')
diff2exp = os.getenv('diff2exp', 'True')

fcHours = os.getenv('fcHours', 6)
intervalHours = os.getenv('intervalHours', 6)
interval = int(intervalHours)/24.
fcNums = os.getenv('fcNums', '1')

SDATE = datetime.strptime(initDate,"%Y%m%d%H")
EDATE = datetime.strptime(endDate,"%Y%m%d%H")

DATEINC = dt.timedelta(days=interval)
expStats = 'expmgfsci'
varNames = ['pressure','theta','uReconstructZonal','uReconstructMeridional','qv']
allmetrics = ['RMS'] #,'MS']
latBands = ['NXTro','Tro','SXTro']

fcRange = int(fcHours)/24.
interval = int(intervalHours)/24.

def readdata():

   if os.path.exists(expStats+'.nc'):
       os.remove(expStats+'.nc')

   varNamesListAll = []
   for metrics in allmetrics:
       for varName in varNames:
           for latBand in latBands:
               for fcTDelta in np.arange(0, fcRange+interval,interval):
                   varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                   varNamesListAll.append(varNamesList)
                   #print('check',varNamesListAll)
   for metrics in allmetrics:
       for varName in varNames:
           for latBand in latBands:
               #for all levels:
               varNamesList2 = ['expmgfs_'+ latBand +'_'+ varName + '_' + metrics]
               varNamesListAll.append(varNamesList2)
   print('check varNamesListAll=',(varNamesListAll))
   for i in range(0,len(varNamesListAll)): 
       alldata = []
       alldatacon = []
       alldiffdata = []
       alldiffdata_rmsdiv = []
       xlabeltime  = []
       TDATE = SDATE 
       while TDATE <= EDATE:
           print('TDATE=',TDATE)
           data1 = []
           diffdata = []
           diffdata_rmsdiv = []
           date = TDATE.strftime('%Y%m%d%H')
           xlabeltime = np.append(xlabeltime,date[4:][:-2])
           nc_file1 = EXP_DIR1+'/'+date+'/diagnostic_stats/expmgfs.nc'
           nc_fid1 = Dataset(nc_file1, "r", format="NETCDF4")
           # data1: exp1-GFSANA 
           data1 = np.array( nc_fid1.variables[''.join(varNamesListAll[i])][:] )
           alldata = np.append(alldata, data1)
           if (diff2exp == 'True'): # and ''.join(varNamesListAll[i])[-3:] == 'rms'):
               nc_file2 = EXP_DIR2+'/'+date+'/diagnostic_stats/expmgfs.nc'
               nc_fid2 = Dataset(nc_file2, "r", format="NETCDF4")
               data2 = np.array( nc_fid2.variables[''.join(varNamesListAll[i])][:] )
               #diffdata, e.g. rms(exp1-GFSANA)-rms(exp2-GFSANA) 
               diffdata = data1 - data2
               alldiffdata = np.append(alldiffdata, diffdata)

           TDATE += DATEINC

       alldiffdata = alldiffdata.reshape(len(xlabeltime),len(diffdata)).T
       newfile = write_stats(alldiffdata,''.join(varNamesListAll[i]))
#TODO: move the following functions to plot_utils.py or basic_plot_functions.py
def write_stats(array_f,varNames):
    lows =[]
    highs = []
    mean = []
    rms = []

    #get dimesion for array_f (function of fcnums or levels)
    x, y = array_f.shape

    #loop for every fcnums or levels:
    for i in range(0,x): 

        STATS = {}     
        STATS['Mean'] = np.nanmean(array_f,axis=1)

        low, high = bootstrap(array_f[i,:], 10000, np.mean, 0.05)
        lows = np.append(lows, low)
        highs = np.append(highs,high) 
    stat = np.column_stack((STATS['Mean'], lows, highs))
 
    if os.path.exists(expStats+'.nc'):
        w_nc_fid = Dataset(expStats+'.nc', 'a', format='NETCDF4')
    else:
        w_nc_fid = Dataset(expStats+'.nc', 'w', format='NETCDF4')
        w_nc_fid.description = "MPAS diagnostics/statistics"
        w_nc_fid.createDimension('levels',55)
        w_nc_fid.createDimension('fcnums',int(fcNums))
        w_nc_fid.createDimension('stat',3)
    if x == 55:
        w_nc_var = w_nc_fid.createVariable(varNames,'f4', ('levels','stat')) 
    else:
        w_nc_var = w_nc_fid.createVariable(varNames,'f4', ('fcnums','stat')) 
    w_nc_fid.variables[varNames][:] = stat
    w_nc_fid.close()

def bootstrap(data, num_samples, statistic, alpha):
    #Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic.
    n = data.size
    idx = npr.randint(0, n, (num_samples, n))
    samples_with_replacement = data[idx]
    stat = np.sort(statistic(samples_with_replacement, 1))
    return (stat[int((alpha/2.0)*num_samples)], stat[int((1-alpha/2.0)*num_samples)])

def main():
    readdata()

if __name__ == '__main__': main()
