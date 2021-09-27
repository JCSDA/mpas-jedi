import os
import sys
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('pdf')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import datetime as dt
from datetime import datetime, timedelta
import numpy.random as npr
import var_utils as vu
import modelsp_utils as mu

expStats = 'expmgfsci.nc'
allmetrics = ['RMS'] #['MS']

def readdata():

    varNamesListAll = []
    varName = []
    for metrics in allmetrics:
        #only for 3dim variables
        for varName in mu.varNames3d:
            for latBand in mu.latBands:
                for fcTDelta in np.arange(0, mu.fcRange+mu.interval,mu.interval):
                    varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                    varNamesListAll.append(varNamesList)
    for metrics in allmetrics:
        #for 2dim and 3dim variables:
        for varName in mu.varNames:
            for latBand in mu.latBands:
                #for all levels:
                varNamesList2 = ['expmgfs_'+ latBand +'_'+ varName + '_' + metrics]
                varNamesListAll.append(varNamesList2)

    for i in range(0,len(varNamesListAll)): 
        nc_file = mu.EXP_DIR2+'/'+mu.endDate+'/diagnostic_stats/expmgfsci.nc'
        nc_fid = Dataset(nc_file, "r", format="NETCDF4")

        data = np.array( nc_fid.variables[''.join(varNamesListAll[i])][:] )
        x, y = data.shape
        t = list(np.arange(0,x,1))
        varName = ''.join(varNamesListAll[i]) 
        if (x == 55 or x == 56):
           region = ''.join(varName.split("_")[2:][:-2])
           var    = '_'.join(varName.split("_")[3:][:-1])
           stats  = ''.join(varName.split("_")[4:])
           plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
           plt.plot(data[:,0],t)
           plt.fill_betweenx( t, data[:,1],data[:,2],alpha=0.3,  linestyle='-.')
           plt.plot([0,0],[0, x], ls="--", c=".3")  # zero line
           plt.xlabel('RMS('+mu.exp2Name+')-RMS('+mu.exp1Name+')',fontsize=15)
           plt.ylabel('Levels',fontsize=15)
        else:
           region = ''.join(varName.split("_")[1:][:-2])
           var    = '_'.join(varName.split("_")[2:][:-1])
           stats  = ''.join(varName.split("_")[3:])
           plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
           plt.plot(t,data[:,0])
           plt.fill_between( t,data[:,1],data[:,2],alpha=0.3,  linestyle='-.')
           plt.plot([0, x],[0,0], ls="--", c=".3")  # zero line
           plt.ylabel('RMS('+mu.exp2Name+')-RMS('+mu.exp1Name+')',fontsize=15)
           plt.xlabel('Lead Time (day)',fontsize=15)
     
        plt.grid(True)
        plt.savefig(''.join(varNamesListAll[i])+'_ci.png',dpi=200,bbox_inches='tight')
        plt.close() 

def main():
    readdata()

if __name__ == '__main__': main()   
