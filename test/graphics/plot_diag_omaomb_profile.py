from netCDF4 import Dataset
import os
import numpy
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes

#for plot cycling run results (on cheyenne):file path and file name example 
#filenc = '/gpfs/fs1/scratch/jban/pandac/DA_3dvar/2018041500/Data/obsout_3dvar_sonde_0000.nc4'
#             EXP_DIR                   EXP_NAME  EXP_TIME             DA_METHOD OBS_TYPE tile
 
#EXP_DIR = os.getenv('EXP_DIR','/gpfs/fs1/scratch/jban/pandac/')

EXP_NAME=os.getenv('EXP_NAME','DA_ctest')
#EXP_NAME=os.getenv('EXP_NAME','DA_3dvar')
#EXP_NAME=os.getenv('EXP_NAME','DA_3denvar')

#EXP_TIME=os.getenv('EXP_TIME','2018041500')

DA_METHOD=os.getenv('DA_METHOD','3dvar')
#DA_METHOD=os.getenv('DA_METHOD','3dvar_bumpcov')
#DA_METHOD=os.getenv('DA_METHOD','3denvar')

OBS_TYPE=os.getenv('OBS_TYPE','sonde')
#OBS_TYPE=os.getenv('OBS_TYPE','aircraft')

VAR_NAME=os.getenv('VAR_NAME','air_temperature')
#VAR_NAME=os.getenv('VAR_NAME','eastward_wind')
#VAR_NAME=os.getenv('VAR_NAME','northward_wind')

#total tiles 
tiles = 1  # singularity
#tiles = 36  # cheyenne

def readdata():

    obsncs = []
    ombncs = []
    omancs = []
    prencs = []
    qcncs  = []

    for tile in range(tiles):

#       for plot cycling run results (on cheyenne):
#       filenc = '/gpfs/fs1/scratch/jban/pandac/DA_3dvar/2018041500/Data/obsout_3dvar_sonde_0000.nc4' 
#       filenc = EXP_DIR + EXP_NAME +'/'+ EXP_TIME +'/Data/'+ 'obsout_'+DA_METHOD+'_'+OBS_TYPE+'_'+'{0:04}'.format(tile)+'.nc4' 

#       for plot ctest results:
#       filenc = obsout_3dvar_sonde_0000.nc4 
        filenc = 'obsout_'+DA_METHOD+'_'+OBS_TYPE+'_'+'{0:04}'.format(tile)+'.nc4'
        print(filenc)

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

#   @EffectiveQC, 1: missing; 0: good; 10: rejected by first-guess check
#   keep data @EffectiveQC=0
    obsncs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN
    ombncs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN
    omancs[numpy.logical_or(qcncs == 1, qcncs == 10)] = numpy.NaN

#   assign bins and calculate rmse:
    RMSEombs = []
    RMSEomas = []
    obsnums  = []
    bins =   [1050.,950.,850.,750.,650.,550.,450.,350.,250.,150.,50.,0.] 
    binsfory=  [1000.,900.,800.,700.,600.,500.,400.,300.,200.,100.,0]
    bins = numpy.asarray(bins)

    for j in range(0, len(bins)-1):  
        obsncbin = deepcopy(obsncs)
        ombncbin = deepcopy(ombncs)
        omancbin = deepcopy(omancs)

        obsncbin[numpy.logical_or(prencs <bins[j+1], prencs >=bins[j])] = numpy.NaN
        ombncbin[numpy.logical_or(prencs <bins[j+1], prencs >=bins[j])] = numpy.NaN
        omancbin[numpy.logical_or(prencs <bins[j+1], prencs >=bins[j])] = numpy.NaN

        RMSEomb = np.sqrt(np.nanmean(ombncbin**2))
        RMSEoma = np.sqrt(np.nanmean(omancbin**2))
        
        obsnum = len(obsncbin)-np.isnan(obsncbin).sum()
        obsnums = np.append(obsnums,obsnum) 
        RMSEombs = np.append(RMSEombs,RMSEomb) 
        RMSEomas = np.append(RMSEomas,RMSEoma)

    plotrmse(RMSEombs,RMSEomas,binsfory,obsnums)

def plotrmse(var1,var2,binsfory,obsnums): 
    fig, ax1 = plt.subplots()
#   reverse left y-axis
    plt.gca().invert_yaxis()
    plt.grid(True)
    ax1.plot(var1,binsfory,'g-*',markersize=5)
    ax1.plot(var2,binsfory,'r-*',markersize=5)
    ax1.set_xlabel('RMSE',fontsize=15)
    ax1.set_xlim([0,6])
#    ax1.set_xlim([0,2])  
    ax1.set_ylim([1000,0])
    major_ticks = np.arange(0, 1000, 100)
    ax1.set_yticks(major_ticks)
    ax1.set_ylabel('Pressure (hPa)',fontsize=15)

#   set right y-axis and reverse it
    ax2 = ax1.twinx()
    ax2.set_yticks(binsfory)
#   reverse right y-axis
    ax2.set_yticklabels(reversed(obsnums.astype(np.int)))
    ax2.set_ylabel('Observation Number(EffectiveQC=0)',fontsize=15) 

    ax1.legend(('OMB','OMA'), loc='lower left',fontsize=15)

    plt.savefig('RMSE_%s_%s_%s.png'%(EXP_NAME,VAR_NAME,OBS_TYPE),dpi=200,bbox_inches='tight')
    #plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
