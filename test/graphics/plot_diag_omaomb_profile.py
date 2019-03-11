from netCDF4 import Dataset
import os
import numpy
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch
import math

def readdata():

    obsoutfiles = []
    for files in os.listdir('../Data/'):
        #print(files)
        if fnmatch.fnmatch(files, 'obsout_*_0000.nc4'):   # 1tile
            obsoutfiles.append('../Data/'+files)
    print(obsoutfiles)

    for file_name in obsoutfiles:
        print(file_name)
        nc = Dataset(file_name, 'r')
        obstype = file_name.split("_")[-2:][:-1]
        #print(obstype)
        expt_obs = ''.join(file_name[15:][:-9])
        #print(expt_obs)
        if ''.join(obstype) in ('sonde','aircraft','satwind'):
             prenc = nc.variables['air_pressure@MetaData']
             prenc = numpy.asarray(prenc)
        varlist = nc.variables.keys()
        #select variables with the suffix 'ObsValue'
        obslist = [obs for obs in varlist if (obs[-8:] == 'ObsValue')]
        #print(obslist)
        for var in obslist:
            varname = ''.join(var.split("@")[:-1])
            omb=''.join(var.split("@")[:-1])+'@ombg'
            oma=''.join(var.split("@")[:-1])+'@oman'
            qc = ''.join(var.split("@")[:-1])+'@EffectiveQC'
            #print("var=",var,"omb=",omb,"oma=",oma)
            obsnc = nc.variables[var]
            ombnc = nc.variables[omb]
            omanc = nc.variables[oma]
            qcnc  = nc.variables[qc]

            obsnc = numpy.asarray(obsnc)
            ombnc = numpy.asarray(ombnc)
            omanc = numpy.asarray(omanc)
            qcnc  = numpy.asarray(qcnc)
            #@EffectiveQC, 1: missing; 0: good; 10: rejected by first-guess check
            #keep data @EffectiveQC=0
            obsnc[numpy.logical_or(qcnc == 1, qcnc == 10)] = numpy.NaN
            ombnc[numpy.logical_or(qcnc == 1, qcnc == 10)] = numpy.NaN
            omanc[numpy.logical_or(qcnc == 1, qcnc == 10)] = numpy.NaN

            if ''.join(obstype) in ('amsua'):
                print('Radiances: Scatter plots will be added soon.')
            elif ''.join(obstype) in ('sonde','aircraft','satwind'):
                #assign bins and calculate rmse:
                RMSEombs = []
                RMSEomas = []
                obsnums  = []
                bins =   [1050.,950.,850.,750.,650.,550.,450.,350.,250.,150.,50.,0.]
                binsfory=  [1000.,900.,800.,700.,600.,500.,400.,300.,200.,100.,0]
                bins = numpy.asarray(bins)

                for j in range(0, len(bins)-1):
                    obsncbin = deepcopy(obsnc)
                    ombncbin = deepcopy(ombnc)
                    omancbin = deepcopy(omanc)

                    obsncbin[numpy.logical_or(prenc <bins[j+1], prenc >=bins[j])] = numpy.NaN
                    ombncbin[numpy.logical_or(prenc <bins[j+1], prenc >=bins[j])] = numpy.NaN
                    omancbin[numpy.logical_or(prenc <bins[j+1], prenc >=bins[j])] = numpy.NaN

                    RMSEomb = np.sqrt(np.nanmean(ombncbin**2))
                    RMSEoma = np.sqrt(np.nanmean(omancbin**2))
        
                    obsnum = len(obsncbin)-np.isnan(obsncbin).sum()
                    obsnums = np.append(obsnums,obsnum)
                    RMSEombs = np.append(RMSEombs,RMSEomb)
                    RMSEomas = np.append(RMSEomas,RMSEoma)

                plotrmsepro(RMSEombs,RMSEomas,binsfory,obsnums,expt_obs,varname)

def plotrmsepro(var1,var2,binsfory,obsnums,EXP_NAME,VAR_NAME):
    fig, ax1 = plt.subplots()
#   reverse left y-axis
    plt.gca().invert_yaxis()
    plt.grid(True)
    ax1.plot(var1,binsfory,'g-*',markersize=5)
    ax1.plot(var2,binsfory,'r-*',markersize=5)
    ax1.set_xlabel('RMSE',fontsize=15)

    if VAR_NAME in 'specific_humidity':
        ax1.set_xlim([0,np.nanmax(var1)])
    else:
        ax1.set_xlim([0,math.ceil(np.nanmax(var1))])
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

    plt.savefig('RMSE_%s_%s.png'%(EXP_NAME,VAR_NAME),dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
