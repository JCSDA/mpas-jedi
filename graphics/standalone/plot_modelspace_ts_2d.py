import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import datetime as dt
from datetime import datetime, timedelta
import basic_plot_functions as BasicPF
import plot_utils as pu
import var_utils as vu
import modelsp_utils as mu
def readdata():

    for fcTDelta in np.arange(0,mu.fcRange+mu.interval,mu.interval):
        varNamesListAll = []
        for metrics in mu.allFileStats: #allmetrics:
            for varName in mu.varNames3d:
                for latBand in mu.latBands:
                        varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                        varNamesListAll.append(varNamesList)
        for i in range(0,len(varNamesListAll)): 
          arraylist = []
          for iexp, expName in enumerate(mu.expNames):

            alldata = []
            alldata2 = []
            alldiffdata = []
            alldiffdata_rmsdiv = []
            xlabeltime  = []
            TDATE = mu.SDATE

            while TDATE <= mu.EDATE:

                data1 = []
                diffdata = []
                diffdata_rmsdiv = []
                date = TDATE.strftime('%Y%m%d%H')
                xlabeltime = np.append(xlabeltime,TDATE)

                nc_file = mu.expDirectory+'/'+mu.expLongNames[iexp]+'/FC1DIAG/'+date+'/diagnostic_stats/expmgfs.nc'
                nc_fid = Dataset(nc_file, "r", format="NETCDF4")
                #data: exp1-GFSANA
                data = np.array( nc_fid.variables[''.join(varNamesListAll[i])][:] )
                alldata = np.append(alldata, data)

                TDATE += mu.DATEINC

            ylevels = list(np.arange(0,len(data),1))
            alldata = alldata.reshape(len(xlabeltime),len(data)).T
            if (iexp == 0):
                arraylist = [alldata]
            else:
                arraylist= arraylist + [alldata]

          dmin = np.amin(arraylist) #min(arraylist.all()) #np.amin(arraylist)
          dmax = np.amax(arraylist) #max(arraylist.all()) #np.amax(arraylist)
          nx = mu.nExp
          ny = 1 #nVars
          nsubplots = nx * ny
          subplot_size = 2.4
          aspect = 0.55
          FULL_SUBPLOT_LABELS =True
          interiorLabels = True
          fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)
          title_tmp = ''.join(varNamesListAll[i])
          region = ''.join(title_tmp.split("_")[2:][:-2])
          var    = ''.join(title_tmp.split("_")[3:][:-1])
          statistic = ''.join(title_tmp.split("_")[4:])
          FCDay  = ''.join(title_tmp.split("_")[1:][:-3])
          xLabel = stats
          if (statistic == 'Mean'):
              centralValue = 0.0
          else:
              centralValue = None
          #print(region,var,statistic)
          yLabel = 'Vertical level'
          sciTicks = False
          logScale = False
          iplot = 0

          for k in arraylist:
                valuemin = np.amin(k)
                valuemax = np.amax(k)
                title = mu.expNames[iplot] +' var:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')\
 '+region+' min=' +str(round(valuemin,3))+' max='+str(round(valuemax,3))

                BasicPF.plot2D(fig,
                    xlabeltime, ylevels, k,
                    title, xLabel, yLabel,
                    BasicPF.defaultIndepConfig,
                    BasicPF.defaultIndepConfig,
                    sciTicks, logScale, centralValue,
                    ny, nx, nsubplots, iplot,
                    dmin = dmin, dmax = dmax,
                    interiorLabels = interiorLabels,
                )

                iplot = iplot +1
                filename = ''.join(varNamesListAll[i])
                pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS, 0.6)

          #plot diff between two exp for RMS:
          if (mu.diff2exp == 'True' and ''.join(varNamesListAll[i])[-3:] == 'RMS'):
              for iexp in range(1,mu.nExp):
              #                 target_exp - control_exp
                  alldiffdata = arraylist[iexp]-arraylist[0]
                  is_all_zero = np.all((alldiffdata == 0))
                  if is_all_zero:
                      print(''.join(varNamesListAll[i]), ': no difference between experiments:', mu.expNames[0], ' and ',mu.expNames[iexp] )
                  else:
                      iplot = 0
                      nx = mu.nExp - 1
                      ny = 1
                      nsubplots = nx * ny
                      signDefinite = False
                      valuemin = np.amin(alldiffdata)
                      valuemax = np.amax(alldiffdata)
                      title = 'Diff var:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+') min=' +str(round(valuemin,3))+' max='+str(round(valuemax,3))
                      fig = pu.setup_fig(nx, ny, subplot_size, aspect, FULL_SUBPLOT_LABELS)

                      BasicPF.plot2D(fig,
                          xlabeltime, ylevels, alldiffdata,
                          title, xLabel, yLabel,
                          BasicPF.defaultIndepConfig,
                          BasicPF.defaultIndepConfig,
                          sciTicks, logScale, centralValue=0.0,
                          ny, nx, nsubplots, iplot,
                          dmin = valuemin, dmax = valuemax,
                          interiorLabels = interiorLabels,
                      )
                      iplot = iplot + 1
                      filename = ''.join(varNamesListAll[i])+mu.expNames[iexp]+'-RMS'+mu.expNames[0] #''.join(varNamesListAll[i])+'_TS_2d'
                      pu.finalize_fig(fig, filename, 'png', FULL_SUBPLOT_LABELS, 0.6)
def main():
    readdata()

if __name__ == '__main__': main()   
