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

# columns:    var_name    unit_used  abbr.
vardict = {'air_temperature': ['(K)','T'],
    'eastward_wind': ['(m/s)','U'],
    'northward_wind': ['(m/s)','V'],
    'specific_humidity': ['(kg/kg)','Q'],
    'refractivity': ['(%)','Ref'],    #plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
    'bending_angle': ['(%)','Bnd'],   #plot RMSE of OMB/O and OMA/O; bending_angle unit: rad
    'brightness_temperature': ['(K)','BT'] }

def readdata():

    print_fmt = 'png' #lower fidelity, faster
    #print_fmt = 'pdf' #higher fidelity, slower

    profile_group  = ['sonde','aircraft','satwind','gnssro']
    radiance_group = ['amsua_n19']
    #dummy_group   = ['dummy_obstype1']

    all_groups = []
    all_groups.append(profile_group)
    all_groups.append(radiance_group)
    #all_groups.append(dummy_group)
    #all_groups.append(['dummy_obstype2'])

    fileprefix = 'obsout_'
    filesuffix = '_*.nc4'
    # * in filesuffix is 4-digit processor rank [0000 to XXXX]
    #npedigits = 4

    # Suffixes for nc variables
    omb_suf = 'ombg'
    oma_suf = 'oman'
    obs_suf = 'ObsValue'
    qc_suf  = 'EffectiveQC'

    obsoutfiles = []
    for files in os.listdir('../Data/'):
        #print(files)
        if fnmatch.fnmatch(files, fileprefix+'*'+filesuffix):   # 1tile
            obsoutfiles.append('../Data/'+files)
    #print(obsoutfiles)

    #Group files by experiment_obstype (e.g., 3dvar_aircraft),
    # where each group contains files from all PE's
    file_groups = [[]]
    for j, file_name in enumerate(obsoutfiles):
        # file_group_name excludes everything outside the first/final '_'
        file_group_name =  '_'.join(file_name.split("_")[1:][:-1])
        for i, file_group in enumerate(file_groups):
            if file_group_name in file_group:
                # If file_group with file_group_name exists, add new file to it
                update_group = file_group
                update_group.append(file_name)
                file_groups[i][:] = update_group
                break
            elif i == len(file_groups)-1:
                # If group with file_group_name does not exist, add one
                new_group = [file_group_name]
                new_group.append(file_name)
                file_groups.append(new_group)
                break
    file_groups = file_groups[1:][:]

    # Loop over unique experiment-obstype combinations
    for file_group in file_groups:
        expt_obs = file_group[0]
        print(expt_obs)

        # Determine obstype from expt_obstype string
        expt_parts = expt_obs.split("_")
        nstr = len(expt_parts)
        obstype = 'none'
        for i in range(0,nstr):
            obstype_ = '_'.join(expt_parts[i:nstr])
            for group in all_groups:
                if ''.join(obstype_) in group:
                    obstype = obstype_

        if obstype == 'none':
           print('obstype not selected, skipping data: '+expt_obs)
           continue

        # Select variables with the suffix omb_suf (required for OMB)
        nc = Dataset(file_group[1], 'r')
        varlist = nc.variables.keys()
        bglist = [var for var in varlist if (var[-4:] == omb_suf)]

        # Define a channel list for radiance_group
        if ''.join(obstype) in radiance_group:
            chlist = []
            for bgname in bglist:
                var = ''.join(bgname.split("@")[:-1])
                chlist.append(int(''.join(var.split("_")[-1:])))
            chlist.sort()

        # Loop over variables with omb suffix
        nvars = len(bglist)
        for ivar, omb in enumerate(bglist):
            varname = ''.join(omb.split("@")[:-1])
            print(varname)
            obs=''.join(omb.split("@")[:-1])+'@'+obs_suf
            oma=''.join(omb.split("@")[:-1])+'@'+oma_suf
            qc = ''.join(omb.split("@")[:-1])+'@'+qc_suf
            #print("obs=",obs,"omb=",omb,"oma=",oma)

            obsnc = np.asarray([])
            ombnc = np.asarray([])
            omanc = np.asarray([])
            qcnc  = np.asarray([])
            prenc = np.asarray([])

            # Build up arrays in loop over file_group, 
            # excluding category in file_group[0]
            for file_name in file_group[1:]:
                nc = Dataset(file_name, 'r')
                #file_rank = file_name[-(4+npedigits):-4]

                if ''.join(obstype) in profile_group:
                    if ''.join(obstype) == 'gnssro':
                        prenc = np.append(prenc, nc.variables['altitude@GroupUndefined'])
                    else:
                        prenc =  np.append(prenc, nc.variables['air_pressure@MetaData'])

                obsnc = np.append( obsnc, nc.variables[obs] )
                ombnc = np.append( ombnc, nc.variables[omb] )
                omanc = np.append( omanc, nc.variables[oma] )
                qcnc  = np.append( qcnc,  nc.variables[qc]  )

            #@EffectiveQC, 1: missing; 0: good; 10: rejected by first-guess check
            #keep data @EffectiveQC=0
            obsnc[numpy.logical_not(qcnc == 0)] = numpy.NaN
            ombnc[numpy.logical_not(qcnc == 0)] = numpy.NaN
            omanc[numpy.logical_not(qcnc == 0)] = numpy.NaN
            ombnc[ombnc < -1.0e+15 ] = numpy.NaN

            if ''.join(obstype) == 'gnssro':
                ombnc = (ombnc/obsnc)*100
                omanc = (omanc/obsnc)*100

            if ''.join(obstype) in profile_group:
                #PROFILE OBS
                #assign bins and calculate rmse:
                RMSEombs = []
                RMSEomas = []
                obsnums  = []

                if ''.join(obstype) == 'gnssro':
                    bins = list(np.arange(50000.0, -1000, -2000.))
                    binsfory= list(np.arange(49500.0,-500.0,-2000.))
                else:
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

                plotrmsepro(RMSEombs,RMSEomas,binsfory,obsnums,expt_obs,varname,print_fmt)
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
                                     nfigtypes, figs, expt_obs,print_fmt)

                if ivar == nvars-1:
                    # Close figs in reverse order to avoid seg fault
                    for fig in reversed(figs):
                        plt.close(fig)


def plotrmsepro(var1,var2,binsfory,obsnums,EXP_NAME,VAR_NAME,fmt):
    fig, ax1 = plt.subplots()
#   reverse left y-axis
    plt.gca().invert_yaxis()
    plt.grid(True)
    ax1.plot(var1,binsfory,'g-*',markersize=5)
    ax1.plot(var2,binsfory,'r-*',markersize=5)

    if VAR_NAME in 'specific_humidity':
        ax1.set_xlim([0,np.nanmax(var1)])
    else:
        ax1.set_xlim([0,math.ceil(np.nanmax(var1))])

    ax2 = ax1.twinx()
    ax2.set_yticks(binsfory)
    ax2.set_ylabel('Observation Number(EffectiveQC=0)',fontsize=15)

    if EXP_NAME[-6:] == 'gnssro':
        ax1.set_ylim([0,49500])
        ax1.set_ylabel('Altitude (m)',fontsize=15)
        ax1.set_xlabel(VAR_NAME+'  RMSE of OMB/O & OMA/O '+ vardict[VAR_NAME][0],fontsize=15)
        ax1.legend(('OMB/O','OMA/O'), loc='upper right',fontsize=15)
        ax2.set_yticklabels(obsnums.astype(np.int))
    else:
        ax1.set_ylim([1000,0])
        major_ticks = np.arange(0, 1000, 100)
        ax1.set_yticks(major_ticks)
        ax1.set_ylabel('Pressure (hPa)',fontsize=15)
        ax1.set_xlabel(VAR_NAME+'  RMSE  '+ vardict[VAR_NAME][0],fontsize=15)
        ax1.legend(('OMB','OMA'), loc='lower left',fontsize=15)
        ax2.set_yticklabels(reversed(obsnums.astype(np.int)))

    fname = 'RMSE_%s_%s.'%(EXP_NAME,VAR_NAME)+fmt
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
    nx_subplt.append(nx)
    ny_subplt.append(ny)

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
                         nfigtypes, figs, EXP_NAME, fmt):
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

        # Uncomment these 2 lines to put x/y labels only on peripheral subplts
        #if kfig <= min(nvars,numsubplts) - nx_subplt[0] : xlab = ''
        #if np.mod(kfig,nx_subplt[0]) != 1 : ylab = ''

        ax = figs[offset+jfig].add_subplot(ny_subplt[0],nx_subplt[0],kfig)

        if iifig == 0:
            #Add scatter plot for h(x) vs. y
            fname = 'XB_XA_%s'%(EXP_NAME)
            scatter_one2ones(obs, [obs-omb , obs-oma],
                             ['x_b' , 'x_a'], True, xlab, ylab,
                             varname, varunits, ax)
        if iifig == 1:
            #Add scatter plot for OMA vs. OMB
            fname = 'OMB_OMA_%s'%(EXP_NAME)
            scatter_one2ones(omb , [oma],
                             [], False, xlab, ylab,
                             varname, varunits, ax)

        if nnfigs > 1: fname=fname+'_%d-of-%d'%(jfig,nnfig)
        fname=fname+'.'+fmt

        if (ivar == nvars-1 or subplt_cnt[jfig] == numsubplts):
            #Save the figure to file
            print('Saving figure to '+fname)
            figs[offset+jfig].subplots_adjust(wspace=0.35,hspace=0.35)
            figs[offset+jfig].savefig(fname,dpi=200,bbox_inches='tight')

        offset += nnfigs

    return subplt_cnt

def scatter_one2ones(XVAL,YVALS,LEG,show_stats,XLAB,YLAB,VAR_NAME,UNITS,ax):
#================================================================
#INPUTS:
# XVAL       - single list of x-coordinates
# YVALS      - list of lists of y-coordinates
# LEG        - list of legend entries
# show_stats - boolean, show slope and RMSE for each line
# XLAB       - xlabel string
# YLAB       - ylabel string
# VAR_NAME   - variable name for text label
# UNITS      - variable units
# ax         - matplotlib.pyplot axes object
#
#OUTPUTS: none, modifies ax object to include one-to-one plots
#
#PURPOSE: Create a one-to-one scatter plot on ax using XVAL and 
#         YVALS, including:
#         + unique markers for each list contained in YVALS
#         + linear regressions for each list contained in YVALS
#         + a one-to-one line
#================================================================

    colors = [ \
              [0.0000, 0.4470, 0.7410], \
              [0.8500, 0.3250, 0.0980], \
              [0.9290, 0.6940, 0.1250], \
              [0.4940, 0.1840, 0.5560], \
              [0.4660, 0.6740, 0.1880], \
              [0.3010, 0.7450, 0.9330], \
              [0.6350, 0.0780, 0.1840], \
             ]
    markers = ['*','+','o','.']
    msizes  = [0.5,0.5,0.5, 3 ]
    for i,YVAL in enumerate(YVALS):
        if len(XVAL) != len(YVAL):
            print('ERROR: Incorrect usage of scatter_one2ones, YVALS must be list of arrays.')
            os._exit()
        col = colors[np.mod(i,len(colors))]
        mind = np.mod(i,len(markers))
        ax.plot( XVAL, YVAL, markers[mind], color = col, \
                 markersize = msizes[mind], alpha=0.5 )

    if XLAB != '':
        label = XLAB
        if UNITS != '': label = label+' ('+UNITS+')'
        ax.set_xlabel(label,fontsize=6)
    if YLAB != '':
        label = YLAB
        if UNITS != '': label = label+' ('+UNITS+')'
        ax.set_ylabel(label,fontsize=6)
    if len(LEG) > 0:
       ax.legend(LEG, loc='upper left',fontsize=5)

    fsize1 = 6.0
    fsize2 = 4.0
    ax.text(0.03, 0.75, VAR_NAME,
        {'color': 'k', 'fontsize': fsize1}, 
        ha='left', va='top', transform=ax.transAxes)


    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    xymin=min(xmin,ymin)
    xymax=max(xmax,ymax)

    # Add linear regression and statistics
    predictor = np.arange(xymin, xymax, (xymax - xymin) / 10.0)
    tx = 0.98
    ty = 0.02
    nline = ''
    for j,YVAL in enumerate(reversed(YVALS)):
        i = len(YVALS) - j - 1
        p = np.polyfit(XVAL,YVAL,1)
        predicted = p[0] * predictor + p[1]
        col0 = colors[np.mod(i,len(colors))]
#        bright = 0.5
#        col = bright * np.asarray([1.,1.,1.]) + (1. - bright) * np.asarray(col0)
        dimmer = 0.35
        col = (1. - dimmer) * np.asarray(col0)

        ax.plot( predictor, predicted, '-', color = col, lw=1.2 )
        stat = 'slope: %0.2f  '%(p[0])

        if show_stats:
            RMSE = np.sqrt( np.sum( np.square(YVAL - XVAL) ) / len(XVAL) )
            BIAS = np.sum( YVAL - XVAL ) / len(XVAL)
            stat = stat+'\nRMSE: %0.2f \nBIAS: %0.2f'%(RMSE,BIAS)
        ax.text(tx, ty, stat+nline, \
            {'color': col, 'fontsize': fsize2}, 
            ha='right', va='bottom', backgroundcolor=[1,1,1,0.2], \
            clip_on=True, transform=ax.transAxes)
        nline = nline + ''.join('\n' * stat.count('\n')) + '\n'

    ax.text(tx, ty, 'N = %d'%(len(XVAL))+nline,
        {'color': 'k', 'fontsize': fsize2}, 
        ha='right', va='bottom', transform=ax.transAxes)


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
