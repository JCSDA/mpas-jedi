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
import var_utils as vu

EXP_DIR1 = os.getenv('FCDIAG_WORK_DIR1','/gpfs/fs1/scratch/jban/pandac/conv60it/FC1/')
EXP_DIR2 = os.getenv('FCDIAG_WORK_DIR2','/gpfs/fs1/scratch/jban/pandac/conv60it/FC1/')

initDate = os.getenv('start_init', '2018041500')
endDate  = os.getenv('end_init', '2018051400')
diff2exp = os.getenv('diff2exp', 'True')

fcHours = os.getenv('fcHours', 6)
intervalHours = os.getenv('intervalHours', 6)

SDATE = datetime.strptime(initDate,"%Y%m%d%H")
EDATE = datetime.strptime(endDate,"%Y%m%d%H")

DATEINC = dt.timedelta(days=0.25)

expStats = 'expmgfsci.nc'
varNames = ['pressure','theta','uReconstructZonal','uReconstructMeridional','qv']

allmetrics = ['RMS'] #['MS']
latBands = ['NXTro','Tro','SXTro']

fcRange = int(fcHours)/24.
interval = int(intervalHours)/24.

def readdata():

    varNamesListAll = []
    varName = []
    for metrics in allmetrics:
        for varName in varNames:
            for latBand in latBands:
                for fcTDelta in np.arange(0, fcRange+interval,interval):
                    
                    varNamesList = ['expmgfs_day'+str(fcTDelta)+'_'+ latBand +'_'+ varName + '_' + metrics]
                    varNamesListAll.append(varNamesList)
                    print('check',varNamesList)
    for metrics in allmetrics:
        for varName in varNames:
            for latBand in latBands:
                #for all levels:
                varNamesList2 = ['expmgfs_'+ latBand +'_'+ varName + '_' + metrics]
                varNamesListAll.append(varNamesList2)

    for i in range(0,len(varNamesListAll)): 
        nc_file1 = EXP_DIR1+'/'+endDate+'/diagnostic_stats/expmgfsci.nc'          
        nc_fid1 = Dataset(nc_file1, "r", format="NETCDF4") 

        data1 = np.array( nc_fid1.variables[''.join(varNamesListAll[i])][:] )
        x, y = data1.shape
        t = list(np.arange(0,x,1))
        varName = ''.join(varNamesListAll[i]) 
        if x == 55: 
           region = ''.join(varName.split("_")[2:][:-2])
           var    = ''.join(varName.split("_")[3:][:-1])
           stats  = ''.join(varName.split("_")[4:])
           print('jban check',varName,region,var,stats)
           plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
           plt.plot(data1[:,0],t)
           plt.fill_betweenx( t, data1[:,1],data1[:,2],alpha=0.3,  linestyle='-.')
           plt.plot([0,0],[0, x], ls="--", c=".3")  # zero line
           plt.xlabel('RMS(AMSUA)-RMS(CONV)',fontsize=15)
           plt.ylabel('Levels',fontsize=15)
        else:
           #works for 10day FC
           region = ''.join(varName.split("_")[1:][:-2])
           var    = ''.join(varName.split("_")[2:][:-1])
           stats  = ''.join(varName.split("_")[3:])
           print('jban check',region,var,stats)
           plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
           plt.plot(t,data1[:,0])
           plt.fill_between( t,data1[:,1],data1[:,2],alpha=0.3,  linestyle='-.')
           plt.plot([0, x],[0,0], ls="--", c=".3")  # zero line
           plt.ylabel('RMS(AMSUA)-RMS(CONV)',fontsize=15)
           plt.xlabel('Lead Time (day)',fontsize=15)
     
        plt.grid(True)
        plt.savefig(''.join(varNamesListAll[i])+'_ci.png',dpi=200,bbox_inches='tight')
        plt.close() 

def main():
    readdata()

if __name__ == '__main__': main()   
