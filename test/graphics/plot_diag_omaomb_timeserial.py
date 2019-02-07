import os
from netCDF4 import Dataset
import numpy
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import datetime as dt

#file path and file name example for cycling run (on cheyenne): 
#filenc = '/gpfs/fs1/scratch/jban/pandac/DA_3dvar/2018041500/Data/obsout_3dvar_sonde_0000.nc4'
#                EXP_DIR                   EXP_NAME      EXP_TIME             DA_METHOD OBS_TYPE tile

EXP_DIR = os.getenv('EXP_DIR','/gpfs/fs1/scratch/jban/pandac/')

SDATE = dt.datetime(2018,4,15,0,0,0)
EDATE = dt.datetime(2018,4,15,18,0,0)
DATEINC = dt.timedelta(days=0.25)

EXP_NAME=os.getenv('EXP_NAME','DA_3dvar')
#EXP_NAME=os.getenv('EXP_NAME','DA_3denvar')

VAR_NAME=os.getenv('VAR_NAME','air_temperature')
#VAR_NAME=os.getenv('VAR_NAME','eastward_wind')
#VAR_NAME=os.getenv('VAR_NAME','northward_wind')

OBS_TYPE=os.getenv('OBS_TYPE','sonde')
#OBS_TYPE=os.getenv('OBS_TYPE','aircraft')

DA_METHOD=os.getenv('DA_METHOD','3dvar')
#DA_METHOD=os.getenv('DA_METHOD','3denvar')

#total tiles 
tiles = 36

def readdata():

    obsnums = []
    RMSEombs = []
    RMSEomas = []
    MEANombs = []
    MEANomas = []
    xlabeltime  = []

    TDATE = SDATE

    while TDATE <= EDATE:
        yyyymmddhh = TDATE.strftime('%Y%m%d%H')
        print(yyyymmddhh)
        obsncs = []
        ombncs = []
        omancs = []
        prencs = []
        qcncs  = []

        for tile in range(tiles):
#           filenc = '/gpfs/fs1/scratch/jban/pandac/DA_3dvar/2018041500/Data/obsout_3dvar_sonde_0000.nc4' 
            filenc = EXP_DIR + EXP_NAME +'/'+ yyyymmddhh +'/Data/'+ 'obsout_'+DA_METHOD+'_'+OBS_TYPE+'_'+'{0:04}'.format(tile)+'.nc4'
#            print(filenc)
            nc = Dataset(filenc, 'r')
        
            prenc = nc.variables['air_pressure@MetaData']
            obsnc = nc.variables['%s@ObsValue'%VAR_NAME]
            ombnc = nc.variables['%s@ombg'%VAR_NAME]
            omanc = nc.variables['%s@oman'%VAR_NAME]
            qcnc  = nc.variables['%s@EffectiveQC'%VAR_NAME]

            prencs = np.append(prencs,prenc) 
            obsncs = np.append(obsncs, obsnc)
            ombncs = np.append(ombncs,ombnc)
            omancs = np.append(omancs,omanc)
            qcncs  = np.append(qcncs,qcnc)       

#       @EffectiveQC, 0: good;  
#                     1: missing;  
#                     10: rejected by first-guess check
#       keep good data
        obsncs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN
        ombncs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN
        omancs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN

        RMSEomb = np.sqrt(np.nanmean(ombncs**2))
        RMSEoma = np.sqrt(np.nanmean(omancs**2))
        meanomb = np.nanmean(ombncs)
        meanoma = np.nanmean(omancs)
        obsnum = len(obsncs)-np.isnan(obsncs).sum()

        obsnums = np.append(obsnums,obsnum)
        RMSEombs = np.append(RMSEombs,RMSEomb)
        RMSEomas = np.append(RMSEomas,RMSEoma)
        MEANombs = np.append(MEANombs,meanomb)
        MEANomas = np.append(MEANomas,meanoma)

        xlabeltime = np.append(xlabeltime,yyyymmddhh)

        TDATE += DATEINC

    plotTimeSerial(obsnums,RMSEombs,RMSEomas,MEANombs,MEANomas,xlabeltime)

def plotTimeSerial(obsnums,RMSEombs,RMSEomas,MEANombs,MEANomas,xlabeltime):

    fig,ax1 = plt.subplots(3,sharex=True)
    xarray = range(len(xlabeltime))
    major_ticks = np.arange(0, len(xlabeltime), 1)

    ax1[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax1[0].plot(xarray,obsnums,'k-*',markersize=5)
    ax1[0].set_xticks(major_ticks)
    ax1[0].set_ylabel('Obs. Num.',fontsize=15)  # obs num for @EffectiveQC=0

    ax1[1].plot(xarray,MEANombs,'g-*',markersize=5)
    ax1[1].plot(xarray,MEANomas,'r-*',markersize=5)
    ax1[1].legend(('OMB','OMA'), loc='upper right',fontsize=10,frameon=False)
    ax1[1].plot([0, len(xlabeltime)-1], [0, 0], ls="--", c=".3")  # zero line
    ax1[1].set_ylabel('MEAN',fontsize=15)

    ax1[2].plot(xarray,RMSEombs,'g-*',markersize=5)
    ax1[2].plot(xarray,RMSEomas,'r-*',markersize=5)
    ax1[2].legend(('OMB','OMA'), loc='upper right',fontsize=10,frameon=False)
    ax1[2].set_ylabel('RMSE',fontsize=15)
    ax1[2].set_xlabel('TIME',fontsize=15)
    ax1[2].set_xticklabels(xlabeltime)

    plt.savefig('TS_%s_%s_%s.png'%(EXP_NAME,VAR_NAME,OBS_TYPE),dpi=200,bbox_inches='tight')
   
def main():
    readdata()

if __name__ == '__main__': main()
