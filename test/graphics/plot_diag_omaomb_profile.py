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
    profile_group = ['sonde','aircraft','satwind']
    radiance_group = ['amsua_n19']

    all_groups = []
    all_groups.append(profile_group)
    all_groups.append(radiance_group)

    obsoutfiles = []
    for files in os.listdir('../Data/'):
        #print(files)
        if fnmatch.fnmatch(files, 'obsout_*_0000.nc4'):   # 1tile
            obsoutfiles.append('../Data/'+files)
    print(obsoutfiles)

    for file_name in obsoutfiles:
        print(file_name)
        nc = Dataset(file_name, 'r')
        expt_parts = file_name.split("_")[1:][:-1]
        nstr = len(expt_parts)
        obstype = 'none'
        for i in range(0,nstr):
            obstype_ = '_'.join(expt_parts[i:nstr])
            for group in all_groups:
                if ''.join(obstype_) in group:
                    obstype = obstype_
        expt_obs = '_'.join(expt_parts)
        #print(obstype)
        #print(expt_obs)
        if obstype == 'none':
           print('obstype not selected, skipping file')
           continue
        if ''.join(obstype) in profile_group:
             prenc = nc.variables['air_pressure@MetaData']
             prenc = numpy.asarray(prenc)
        varlist = nc.variables.keys()

        #select variables with the suffix 'ombg' (required for OMB)
        bglist = [var for var in varlist if (var[-4:] == 'ombg')]
        #print(bglist)
        for omb in bglist:
            varname = ''.join(omb.split("@")[:-1])
            obs=''.join(omb.split("@")[:-1])+'@ObsValue'
            oma=''.join(omb.split("@")[:-1])+'@oman'
            qc = ''.join(omb.split("@")[:-1])+'@EffectiveQC'
            #print("obs=",obs,"omb=",omb,"oma=",oma)
            obsnc = nc.variables[obs]
            ombnc = nc.variables[omb]
            omanc = nc.variables[oma]
            qcnc  = nc.variables[qc]

            obsnc = numpy.asarray(obsnc)
            ombnc = numpy.asarray(ombnc)
            omanc = numpy.asarray(omanc)
            qcnc  = numpy.asarray(qcnc)
            #@EffectiveQC, 1: missing; 0: good; 10: rejected by first-guess check
            #keep data @EffectiveQC=0
            obsnc[numpy.logical_not(qcnc == 0)] = numpy.NaN
            ombnc[numpy.logical_not(qcnc == 0)] = numpy.NaN
            omanc[numpy.logical_not(qcnc == 0)] = numpy.NaN

            if ''.join(obstype) in radiance_group:
                print('Radiances: Scatter plots will be added soon.')
            elif ''.join(obstype) in profile_group:
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
            else:
                print('obstype = '+obstype+' not supported yet')


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

    fname = 'RMSE_%s_%s.png'%(EXP_NAME,VAR_NAME)
    print('Saving figure to '+fname)
    plt.savefig(fname,dpi=200,bbox_inches='tight')
    plt.close()
   
def main():
    readdata()

if __name__ == '__main__': main()
