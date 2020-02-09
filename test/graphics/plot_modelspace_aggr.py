import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from copy import deepcopy
import basic_plot_functions as BasicPF

EXP_DIR1 = os.getenv('FCDIAG_WORK_DIR1','FC1DIAG DIR OR FC2DIAG DIR FOR EXP1')
EXP_DIR2 = os.getenv('FCDIAG_WORK_DIR2','FC1DIAG DIR OR FC2DIAG DIR FOR CONTROL')

diff2exp = os.getenv('diff2exp', 'False')

varNames = ['pressure','theta','uReconstructZonal','uReconstructMeridional','qv']

allmetrics = ['Mean','RMS','STD'] #,'Min','Max']
latBands = ['NXTro','Tro','SXTro']

def readdata():
   fcRange = 10

   nc_file1 = EXP_DIR1+'/expmgfs.nc'
   nc_fid1 = Dataset(nc_file1, "r", format="NETCDF4")
   if (diff2exp == 'True'):
      nc_file2 = EXP_DIR2+'/expmgfs.nc'
      nc_fid2 = Dataset(nc_file2, "r", format="NETCDF4")

   for metrics in allmetrics:
       for varName in varNames:
           for latBand in latBands:
               xlabeltime  = []
               alldata = []
               alldiffdata = []
               alldiffdata_rmsdiv = []
               for fcTDelta in range(0, fcRange+1):
                   varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                   # data1: exp1-GFSANA
                   data1 = np.array( nc_fid1.variables[''.join(varNamesList)][:])
                   alldata = np.append(alldata, data1)
                   xlabeltime = np.append(xlabeltime,fcTDelta)
                   if (diff2exp == 'True' and ''.join(varNamesList)[:][-3:] == 'RMS'):
                       data2 = np.array( nc_fid2.variables[''.join(varNamesList)][:])
                       diffdata = data1 - data2
                       alldiffdata = np.append(alldiffdata, diffdata)
                       print('(expt-control)/control for RMSE,')
                       diffdata_rmsdiv = (data1-data2)/data2 *100
                       alldiffdata_rmsdiv = np.append(alldiffdata_rmsdiv, diffdata_rmsdiv)

               alldata = alldata.reshape(11,len(data1)).T
               ylevels = list(np.arange(0,len(data1),1))
               varNamesListUse = 'expmgfs_fc_'+ latBand +'_'+ varName + '_' + metrics 
               BasicPF.plotTimeserial2D(alldata,xlabeltime,ylevels,varNamesListUse)
                
               if (diff2exp == 'True' and ''.join(varNamesList)[:][-3:] == 'RMS'):
                   alldiffdata = alldiffdata.reshape(11,len(diffdata)).T
                   varNamesListUse = 'expmgfs_fc_'+ latBand +'_'+ varName + '_' + metrics+'diff2exp'
                   BasicPF.plotTimeserial2D(alldiffdata,xlabeltime,ylevels,varNamesListUse)

                   alldiffdata_rmsdiv = alldiffdata_rmsdiv.reshape(11,len(diffdata_rmsdiv)).T
                   varNamesListUse = 'expmgfs_fc_'+ latBand +'_'+ varName + '_' + metrics+'dividediff2exp'
                   BasicPF.plotTimeserial2D(alldiffdata_rmsdiv,xlabeltime,ylevels,varNamesListUse)

def main():
    readdata()

if __name__ == '__main__': main()
