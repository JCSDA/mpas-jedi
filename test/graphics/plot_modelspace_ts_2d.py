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
import basic_plot_functions as BasicPF

EXP_DIR1 = os.getenv('FCDIAG_WORK_DIR1','FC1DIAG DIR OR FC2DIAG DIR FOR EXP1')
EXP_DIR2 = os.getenv('FCDIAG_WORK_DIR2','FC1DIAG DIR OR FC2DIAG DIR FOR CONTROL')

initDate = os.getenv('start_init', '2018041500')
endDate  = os.getenv('end_init', '2018051400')
diff2exp = os.getenv('diff2exp', 'False')

fcHours = os.getenv('fcHours', 6)
intervalHours = os.getenv('intervalHours', 6)

SDATE = datetime.strptime(initDate,"%Y%m%d%H")
EDATE = datetime.strptime(endDate,"%Y%m%d%H")
DATEINC = dt.timedelta(days=0.25)

varNames = ['pressure','theta','qv','uReconstructZonal','uReconstructMeridional']
allmetrics =['Mean','RMS','MS','STD'] #,'Min','Max']
latBands = ['NXTro','Tro','SXTro']

fcRange = int(fcHours)/24.
interval = int(intervalHours)/24.

def readdata():

    for fcTDelta in np.arange(0,fcRange+interval,interval):
        varNamesListAll = []
        for metrics in allmetrics:
            for varName in varNames:
                for latBand in latBands:
                    for fcTDelta in np.arange(0, fcRange+interval,interval):
                    
                        varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                        varNamesListAll.append(varNamesList)
        for i in range(0,len(varNamesListAll)): 
            alldata = []
            alldiffdata = []
            alldiffdata_rmsdiv = []
            xlabeltime  = []
            TDATE = SDATE 
            while TDATE <= EDATE:
                data1 = []
                diffdata = []
                diffdata_rmsdiv = []
                date = TDATE.strftime('%Y%m%d%H')
                xlabeltime = np.append(xlabeltime,date[4:][:-2])
                nc_file1 = EXP_DIR1+'/'+date+'/diagnostic_stats/expmgfs.nc'          
                nc_fid1 = Dataset(nc_file1, "r", format="NETCDF4") 
                #data1: exp1-GFSANA 
                data1 = np.array( nc_fid1.variables[''.join(varNamesListAll[i])][:] )
                alldata = np.append(alldata, data1)

                if (diff2exp == 'True' and ''.join(varNamesListAll[i])[-3:] == 'RMS'):
                    nc_file2 = EXP_DIR2+'/'+date+'/diagnostic_stats/expmgfs.nc'
                    nc_fid2 = Dataset(nc_file2, "r", format="NETCDF4")
                    data2 = np.array( nc_fid2.variables[''.join(varNamesListAll[i])][:] )
                    #diffdata: (exp1-GFSANA)-(exp2-GFSANA) 
                    diffdata = data1 - data2
                    alldiffdata = np.append(alldiffdata, diffdata)

                    diffdata_rmsdiv = (data1-data2)/data2 *100
                    alldiffdata_rmsdiv = np.append(alldiffdata_rmsdiv, diffdata_rmsdiv)

                TDATE += DATEINC
            ylevels = list(np.arange(0,len(data1),1))
            alldata = alldata.reshape(len(xlabeltime),len(data1)).T
            BasicPF.plotTimeserial2D(alldata,xlabeltime,ylevels,''.join(varNamesListAll[i]))
       
            if (diff2exp == 'True' and ''.join(varNamesListAll[i])[-3:] == 'RMS'):
                alldiffdata = alldiffdata.reshape(len(xlabeltime),len(diffdata)).T
                BasicPF.plotTimeserial2D(alldiffdata,xlabeltime,ylevels,''.join(varNamesListAll[i])+'diff2exp')
                alldiffdata_rmsdiv = alldiffdata_rmsdiv.reshape(len(xlabeltime),len(data1)).T
                BasicPF.plotTimeserial2D(alldiffdata_rmsdiv,xlabeltime,ylevels,''.join(varNamesListAll[i])+'dividediff2exp')

def main():
    readdata()

if __name__ == '__main__': main()   
