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
    profile_group  = ['sonde','aircraft','satwind']
    radiance_group = ['amsua_n19']
    #dummy_group   = ['dummy_obstype1']

    all_groups = []
    all_groups.append(profile_group)
    all_groups.append(radiance_group)
    #all_groups.append(dummy_group)
    #all_groups.append(['dummy_obstype2'])

    fileprefix = 'obsout_'
    filesuffix = '_0000.nc4'
    obsoutfiles = []
    for files in os.listdir('../Data/'):
        #print(files)
        if fnmatch.fnmatch(files, fileprefix+'*'+filesuffix):   # 1tile
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

        if ''.join(obstype) in radiance_group:
            chlist = []
            for bgname in bglist:
                var = ''.join(bgname.split("@")[:-1])
                chlist.append(int(''.join(var.split("_")[-1:])))
            chlist.sort()

        nvars = len(bglist)
        for ivar,omb in enumerate(bglist):
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

            if ''.join(obstype) in profile_group:
                #PROFILE OBS
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
                # Maximum number of variables/channels per figure
                maxsubplts = 12

                if ''.join(obstype) in radiance_group:
                    #RADIANCE OBS
                    maxsubplts = 16
                    ch = ''.join(varname.split("_")[-1:])
                    ifig = chlist.index(int(ch)) + 1
                    shortname = 'BT, ch. '+ch
                    units = 'K'
                else:
                    print('NOTIFICATION: obstype = '+obstype+' not defined, using default scatter plots')
                    #Generic scatter verification for any obstype
                    ## (similar to radiance_group):
                    ifig = ivar
                    shortname = varname
                    units = ''

                if ivar == 0:
                    # scatter_verification yields 2 figure types
                    nfigtypes = 2
                    nx_subplt, ny_subplt, figs, subplt_cnt = \
                       init_subplts(bglist, nfigtypes, maxsubplts)

                subplt_cnt = \
                scatter_verification(ifig, shortname, units, ivar, nvars, \
                                     maxsubplts, subplt_cnt, \
                                     obsnc, ombnc, omanc, \
                                     nx_subplt, ny_subplt, \
                                     nfigtypes, figs, expt_obs)

                if ivar == nvars-1:
                    # Close figs in reverse order to avoid seg fault
                    for fig in reversed(figs):
                        plt.close(fig)


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


def init_subplts(subpltlist, nfigtypes, maxsubplts):
#================================================================
#INPUTS:
# subpltlist - a list object defining all subplots required
# nfigtypes  - the number of figure types for each subpltlist member
# maxsubplts - maximum number of subplot objects in a figure type
#
#OUTPUTS: (all to be used externally)
# figs       - list of matplotlib.pyplot figure objects (figs)
# nx_subplt, - maximum dim. of subplots in all figs members
# ny_subplt
# subplt_cnt - counter of subplots in each figure object
#
#PURPOSE: Initialize multiple figure objects and descriptors
#================================================================

    nnfigs = (len(subpltlist) / maxsubplts + 1)
    nsubplts = len(subpltlist)
    if nsubplts > maxsubplts : nsubplts = maxsubplts

    ny_subplt = []
    nx_subplt = []

    nx = np.ceil(np.sqrt(nsubplts))
    ny = np.ceil(np.true_divide(nsubplts,nx))

    ny_subplt.append(nx)
    nx_subplt.append(ny)

    figs = []
    for ifig in range(0,nnfigs*nfigtypes):
        fig = plt.figure()
        inch_size = 1.9
        fig.set_size_inches(nx_subplt[0]*inch_size,ny_subplt[0]*inch_size)
        figs.append(fig)
    subplt_cnt = np.zeros(nnfigs*nfigtypes)

    return nx_subplt, ny_subplt, figs, subplt_cnt

def scatter_verification(ifig, varname, varunits, ivar, nvars, \
                         maxsubplts, subplt_cnt, \
                         obs, omb, oma, \
                         nx_subplt, ny_subplt, \
                         nfigtypes, figs, EXP_NAME):
#================================================================
#INPUTS: 
# ifig       - subplot number
# varname    - variable name
# varunits   - variable units
# ivar       - variable number
# nvars      - total number of variables
# maxsubplts - maximum number of subplots per figure
# subplt_cnt - counter of subplots in all figs members
# obs        - single list of Observation
# omb        - single list of Observation minus background
# oma        - single list of Observation minus analysis
# nx_subplt  - subplot x-dimension of each figs member
# ny_subplt  - subplot y-dimension of each figs member
# nfigtypes  - number of figure types associated with figs list
# figs       - list of matplotlib.pyplot figure objects
# EXP_NAME   - experiment name
#
#OUTPUT: subplt_cnt - updated counter
#
#PURPOSE: Generate verification subplots for varname amongst
#         subplots containing multiple variables
#================================================================

    nnfigs = len(figs) / nfigtypes
    jfig = (ifig-1)/maxsubplts
    kfig = np.mod(ifig-1,maxsubplts)+1
    subplt_cnt[jfig] += 1

    if jfig == nnfigs-1 :
        numsubplts = np.mod(nvars,maxsubplts)
    else :
        numsubplts = maxsubplts

    offset = 0

    for iifig in range(nfigtypes):
        if iifig == 0:
            xlab = 'y'
            ylab = 'h(x)'
        elif iifig == 1:
            xlab = 'OMB'
            ylab = 'OMA'
        else:
            print('WARNING: scatter_verification has no definitions for nfigtypes == ',nfigtypes)
            continue

        #Comment these 2 lines to put x/y labels on all subplts
        if kfig <= min(nvars,numsubplts) - nx_subplt[0] : xlab = ''
        if np.mod(kfig,nx_subplt[0]) != 1 : ylab = ''

        ax = figs[offset+jfig].add_subplot(nx_subplt[0],ny_subplt[0],kfig)

        if iifig == 0:
            #Add scatter plot for h(x) vs. y
            fname = 'XB_XA_%s_%d.png'%(EXP_NAME,jfig)
            scatter_one2ones(obs, [obs-omb , obs-oma],
                             ['x_b' , 'x_a'], xlab, ylab,
                             varname, varunits, ax)
        if iifig == 1:
            #Add scatter plot for OMA vs. OMB
            fname = 'OMB_OMA_%s_%d.png'%(EXP_NAME,jfig)
            scatter_one2ones(omb , [oma],
                             [], xlab, ylab,
                             varname, varunits, ax)

        if (ivar == nvars-1 or subplt_cnt[jfig] == numsubplts):
            #Save the figure to file
            print('Saving figure to '+fname)
            figs[offset+jfig].subplots_adjust(wspace=0.35,hspace=0.35)
            figs[offset+jfig].savefig(fname,dpi=200,bbox_inches='tight')

        offset += nnfigs

    return subplt_cnt

def scatter_one2ones(XVAL,YVALS,LEG,XLAB,YLAB,VAR_NAME,UNITS,ax):
#================================================================
#INPUTS:
# XVAL     - single list of x-coordinates
# YVALS    - list of lists of y-coordinates
# LEG      - list of legend entries
# XLAB     - xlabel string
# YLAB     - ylabel string
# VAR_NAME - variable name for text label
# UNITS    - variable units
# ax       - matplotlib.pyplot axes object
#
#OUTPUTS: none
#
#PURPOSE: Create a one-to-one scatter plot on ax using XVAL and 
#         YVALS, including:
#         + unique markers for each list contained in YVALS
#         + linear regressions for each list contained in YVALS
#         + a one-to-one line
#================================================================

    colors = ['g','r','b','c','m','y','k']
    markers = ['*','+','.','o']
    for i,y in enumerate(YVALS):
        if len(XVAL) != len(y):
            print('ERROR: Incorrect usage of scatter_one2ones, YVALS must be list of arrays.')
            os._exit()
        ax.plot(XVAL,y,colors[np.mod(i,len(colors))]+markers[i/len(colors)],markersize=2)
        #It would be nice to add a linear regression here as well.

    if XLAB != '':
        ax.set_xlabel(XLAB+' ('+UNITS+')',fontsize=6)
    if YLAB != '':
        ax.set_ylabel(YLAB+' ('+UNITS+')',fontsize=6)
    ax.text(0.97, 0.03, VAR_NAME,
        {'color': 'k', 'fontsize': 8}, 
        ha='right', va='bottom', transform=ax.transAxes)

    if len(LEG) > 0:
       ax.legend(LEG, loc='upper left',fontsize=5)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    xymin=min(xmin,ymin)
    xymax=max(xmax,ymax)

    predictor = np.arange(xymin, xymax, (xymax - xymin) / 10.0)
    for i,y in enumerate(YVALS):
        p = np.polyfit(XVAL,y,1)
        predicted = p[1] + p[0] * predictor
        ax.plot(predictor,predicted,colors[np.mod(i,len(colors))]+'-',lw=0.5)

    #Could try to adjust limits/ticks for better aesthetics
    #round_nmbr = 5
    #xymin=np.floor(min(xmin,ymin) / round_nmbr) * round_nmbr
    #xymax=np.ceil(max(xmax,ymax) / round_nmbr) * round_nmbr

    ax.plot([xymin,xymax],[xymin,xymax],'k--',lw=0.5)
    ticks = ax.get_yticks()
    ticks = ticks[np.logical_and(ticks >= xymin, ticks<=xymax)]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(5) 
        tick.label.set_rotation('vertical')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(5) 

    ax.set_xlim(xymin,xymax)
    ax.set_ylim(xymin,xymax)
    ax.grid()

def main():
    readdata()

if __name__ == '__main__': main()
