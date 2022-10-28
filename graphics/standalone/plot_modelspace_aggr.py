import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import basic_plot_functions as BasicPF
import plot_utils as pu
import var_utils as vu
import modelsp_utils as mu

def readdata():

   for metrics in mu.allFileStats:
       for varName in mu.varNames3d:
           for latBand in mu.latBands:
             for iexp, expName in enumerate(mu.expNames):
               xlabeltime  = []
               alldata = []
               alldiffdata = []
               alldiffdata_rmsdiv = []
               nc_file =  mu.expDirectory+'/'+mu.expLongNames[iexp]+'/FC2DIAG/expmgfs.nc'
               nc_fid = Dataset(nc_file, "r", format="NETCDF4")
               for fcTDelta in np.arange(0,mu.fcRange+mu.interval,mu.interval):
                   varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                   # data: exp-GFSANA
                   data = np.array( nc_fid.variables[''.join(varNamesList)][:])
                   alldata = np.append(alldata, data)
                   xlabeltime = np.append(xlabeltime,fcTDelta)

               alldata = alldata.reshape(len(xlabeltime),len(data)).T
               ylevels = list(np.arange(0,len(data),1))
               varNamesListUse = 'expmgfs_fc_'+ latBand +'_'+ varName + '_' + metrics 
               if (iexp == 0):
                   arraylist = [alldata]
               else:
                   arraylist= arraylist + [alldata]
             dmin = np.amin(arraylist) #min(arraylist.all()) #np.amin(arraylist)
             dmax = np.amax(arraylist)
             nx = mu.nExp
             ny = 1 #nVars
             nsubplots = nx * ny
             subplot_size = 2.4
             aspect = 0.55
             FULL_SUBPLOT_LABELS =True
             interiorLabels = True
             fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)
             title_tmp = varNamesListUse #''.join(varNamesListAll[i])
             region = ''.join(title_tmp.split("_")[2:][:-2])
             var    = ''.join(title_tmp.split("_")[3:][:-1])
             stats  = ''.join(title_tmp.split("_")[4:])
             statDiagLabel = stats
             if (stats == 'Mean'):
                 centralValue = 0.0
             else:
                 centralValue = None

             #print(region,var,stats)
             indepLabel = 'Vertical level'
             sciTicks = False
             logScale = False
             iplot = 0

             for k in arraylist:
                valuemin = np.amin(k)
                valuemax = np.amax(k)
                title = mu.expNames[iplot]  +' var:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')\
 '+region+' min=' +str(round(valuemin,3))+' max='+str(round(valuemax,3))

                BasicPF.plot2D(fig,
                    xlabeltime,ylevels, k,
                    title, statDiagLabel, indepLabel,
                    BasicPF.defaultIndepConfig,
                    BasicPF.defaultIndepConfig,
                    sciTicks, logScale, centralValue,
                    ny, nx, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels,
                )
                iplot = iplot +1
                filename = varNamesListUse+'_TS_2d'
                pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS, 0.6)

             #plot diff between target and control exp for RMS:
             if (mu.diff2exp == 'True' and varNamesListUse[-3:] == 'RMS'):
                 for iexp in range(1,mu.nExp):
                     #             target_exp - control_exp
                     alldiffdata = arraylist[iexp]-arraylist[0]
                     BasicPF.plotTimeserial2D(alldiffdata,xlabeltime,ylevels,varNamesListUse+mu.expNames[iexp]+'-RMS'+mu.expNames[0])

def main():
    readdata()

if __name__ == '__main__': main()
