import os
import numpy as np
from copy import deepcopy
import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch
import math
import basic_plot_functions
import h5py as h5
import config as conf
import var_utils as vu

'''
Directory Structure for ctest:
test/
 ├── Data/os/
 │        ├── obsout_3denvar_bump_sondes_0000.nc4
 │        ├── obsout_3denvar_bump_sondes_0001.nc4
 │        ├── ...
 ├── graphics/
 │   ├── plot_diag_omaomb.py
 │   ├── basic_plot_functions.py
 │   ├── ...

Directory Structure for cycling:
test/
 ├── Data/
 │   ├── obsout_3denvar_bump_sondes_0000.h5
 │   ├── obsout_3denvar_bump_sondes_0001.h5
 │   ├── ...
 ├── graphics/
 │   ├── plot_diag_omaomb.py
 │   ├── basic_plot_functions.py
 │   ├── ...
'''

# columns: var_name            unit_used   abbr.
vardict = { \
    'air_temperature':        [ '(K)',     'T'   ] \
  , 'virtual_temperature':    [ '(K)',     'Tv'  ] \
  , 'eastward_wind':          [ '(m/s)',   'U'   ] \
  , 'northward_wind':         [ '(m/s)',   'V'   ] \
  , 'specific_humidity':      [ '(kg/kg)', 'Q'   ] \
  , 'refractivity':           [ '(%)',     'Ref' ] \
  , 'bending_angle':          [ '(%)',     'Bnd' ] \
  , 'brightness_temperature': [ '(K)',     'BT'  ] \
  , 'surface_pressure':       [ '(Pa)',    'Ps'  ] \
    }
# need add Ps
    #Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
    #Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bending_angle unit: rad

def readdata():

    print_fmt = 'png' #lower fidelity, faster
    #print_fmt = 'pdf' #higher fidelity, slower

    plot_allinOneDistri = True   # plot profile_group (includes all levels) and sfc_group.
    plot_eachLevelDistri = True  # plot every level separately for profile_group.

    profile_group  = [
        'sondes',
        'aircraft',
        'satwind',
        'satwnd',
        'gnssroref',
        'gnssrobndropp1d',
    ]
    sfc_group  = [
        'sfc',
    ]
    radiance_group = [
        'abi_g16',
        'ahi_himawari8',
        'amsua_aqua',
        'amsua_metop-a',
        'amsua_n15',
        'amsua_n18',
        'amsua_n19',
        'amsua_n19--hydro',
        'amsua_n19--nohydro',
    ]

    all_groups = []
    all_groups.append(profile_group)
    all_groups.append(sfc_group)
    all_groups.append(radiance_group)

    # Diagnostic oma/omb files located in diagdir
    diagdir    = '../Data/'
    diagprefix = 'obsout_'
    diagsuffix = '_*.nc4'  #for ctests 
    #diagsuffix = '_*.h5'   #for cycling
    # * in diagsuffix is 4-digit processor rank [0000 to XXXX]
    #npedigits = 4

    # Suffixes to required nc variables
    #observations
    obs_var = 'ObsValue'

    #departures
    depbg_var = 'depbg'
    depan_var = 'depan'

    #quality control
    NOUTER=str(os.getenv('NOUTER',2)) #set equal to number of outer iterations
    qcbg_var  = 'EffectiveQC0' #EffectiveQCi, where i is the iteration for depbg_var
    qcan_var  = 'EffectiveQC'+NOUTER #EffectiveQCi, where i is the iteration for depan_var

    err_obs  = 'ObsError'
    errstart_var = 'EffectiveError0'
    errend_var = 'EffectiveError'+NOUTER
    obsoutfiles = []
    for files in os.listdir(diagdir):
        #print(files)
        if fnmatch.fnmatch(files, diagprefix+'*'+diagsuffix):   # 1tile
            obsoutfiles.append(diagdir+files)
    #print(obsoutfiles)

    # Group files by experiment-obstype combination
    #  (e.g., 3dvar_aircraft), where each group 
    #  contains files from all PE's
    exob_groups = [[]]
    for j, file_name in enumerate(obsoutfiles):
        # exob_group_name excludes everything outside the first/final '_'
        exob_group_name =  '_'.join(file_name.split("_")[1:][:-1])
        for i, exob_group in enumerate(exob_groups):
            if exob_group_name in exob_group:
                # If exob_group with exob_group_name exists, add new file to it
                update_group = exob_group
                update_group.append(file_name)
                exob_groups[i][:] = update_group
                break
            elif i == len(exob_groups)-1:
                # If group with exob_group_name does not exist, add one
                new_group = [exob_group_name]
                new_group.append(file_name)
                exob_groups.append(new_group)
                break
    exob_groups = exob_groups[1:][:]

    # Loop over unique experiment-obstype groups
    for exob_group in exob_groups:
        expt_obs = exob_group[0]
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

        # Select variables with the suffix depbg_var (required for OMB)
        nc = h5.File(exob_group[1], 'r')
        varlist = []
        for node in nc:
          if type(nc[node]) is h5._hl.group.Group:
            for var in nc[node]:
              varlist += [node+'/'+var]

        bglist = [var for var in varlist if (var[:][:5] == depbg_var)]
        # Define a channel list for radiance_group
        if ''.join(obstype) in radiance_group:
            chlist = []
            for bgname in bglist:
                var = ''.join(bgname.split("/")[1:])
                chlist = nc['nchans']
            ObsSpaceInfo = conf.DiagSpaceConfig.get(obstype,conf.nullDiagSpaceInfo)
            channels = ObsSpaceInfo.get('channels',[vu.miss_i])

        # Loop over variables with omb suffix
        nvars = len(bglist)
        for ivar, depbg in enumerate(bglist):
            varname = ''.join(depbg.split("/")[1:])
            obs=obs_var+'/'+''.join(depbg.split("/")[1:])
            depan=depan_var+'/'+''.join(depbg.split("/")[1:])
            qcb = qcbg_var+'/'+''.join(depbg.split("/")[1:])
            qca = qcan_var+'/'+''.join(depbg.split("/")[1:])
            obserror = err_obs+'/'+''.join(depbg.split("/")[1:])
            errstart = errstart_var+'/'+''.join(depbg.split("/")[1:])
            errend   = errend_var+'/'+''.join(depbg.split("/")[1:])
            print("obs=",obs,"depbg=",depbg,"depan=",depan)

            obsnc = np.asarray([])
            bkgnc = np.asarray([])
            ananc = np.asarray([])
            ombnc = np.asarray([])
            omanc = np.asarray([])
            qcbnc = np.asarray([])
            qcanc = np.asarray([])
            prenc = np.asarray([])
            latnc = np.asarray([])
            lonnc = np.asarray([])
            recordNum_baknc = np.asarray([])
            recordNum_ananc = np.asarray([])
            stationId_baknc = np.asarray([])
            stationId_ananc = np.asarray([])
            obserrornc = np.asarray([])
            errstartnc = np.asarray([])
            errendnc = np.asarray([])
            # Build up arrays in loop over exob_group, 
            # excluding category in exob_group[0]
            for file_name in exob_group[1:]:
                nc = h5.File(file_name, 'r')
                if ''.join(obstype)[:6] == 'gnssro':
                    prenc = np.append(prenc, nc['MetaData/altitude'])
                    station_id = nc['MetaData/occulting_sat_id']
                    recordNum = nc['MetaData/record_number']
                    recordNum_baknc = np.append( recordNum_baknc, recordNum)
                    recordNum_ananc = np.append( recordNum_ananc, recordNum)
                elif  ''.join(obstype) not in radiance_group:
                    tmp = nc['MetaData/air_pressure']
                    if np.max(tmp) > 10000.0:
                        tmp = np.divide(tmp,100.0)
                    prenc = np.append( prenc, tmp )
                    station_id = nc['MetaData/station_id']
                    station_id = [ ''.join(i[0:5]) for i in nc['MetaData/station_id'][:]]
                    stationId_baknc = np.append( stationId_baknc, station_id )
                    stationId_ananc = np.append( stationId_ananc, station_id )
                obsnc = np.append( obsnc, nc[obs] )
                ombnc = np.append( ombnc, np.negative( nc[depbg] ) ) # omb = (-) depbg
                omanc = np.append( omanc, np.negative( nc[depan] ) ) # oma = (-) depan
                qcbnc  = np.append( qcbnc,  nc[qcb]  )
                qcanc  = np.append( qcanc,  nc[qca]  )
                obserrornc = np.append(obserrornc, nc[obserror])
                errstartnc = np.append(errstartnc, nc[errstart])
                errendnc = np.append(errendnc, nc[errend])
                latnc = np.append( latnc, nc['MetaData/latitude'] )
                lonnc = np.append( lonnc, nc['MetaData/longitude'] )
                bkgnc = obsnc - ombnc
                ananc = obsnc - omanc

            #@EffectiveQC, 10: missing; 0: good; 19: rejected by first-guess check
            #keep data @EffectiveQC=0
            obsnc[qcbnc != 0] = np.NaN
            ombnc[np.less(ombnc,-1.0e+15)] = np.NaN
            ombnc[qcbnc != 0] = np.NaN
            omanc[np.less(omanc,-1.0e+15)] = np.NaN
            omanc[qcanc != 0] = np.NaN
            bkgnc[np.less(bkgnc,-1.0e+15)] = np.NaN
            bkgnc[qcbnc != 0] = np.NaN
            ananc[np.less(ananc,-1.0e+15)] = np.NaN
            ananc[qcanc != 0] = np.NaN
            obserrornc[qcbnc != 0] = np.NaN
            errstartnc[qcbnc != 0] = np.NaN
            errendnc[qcanc != 0] = np.NaN
            if ''.join(obstype) not in radiance_group:
                stationId_baknc[qcbnc != 0] = np.NaN
                stationId_ananc[qcanc != 0] = np.NaN
                if ''.join(obstype)[:6] == 'gnssro':
                    recordNum_baknc[qcbnc != 0] = np.NaN
                    recordNum_ananc[qcanc != 0] = np.NaN
                if np.isnan(ombnc).all() == True:
                    print('all values are NaN in', varname, obstype)
                    continue

            if ''.join(obstype)[:6] == 'gnssro':
                ombnc = (ombnc/obsnc)*100
                omanc = (omanc/obsnc)*100

            if ''.join(obstype) not in radiance_group:
                if (plot_allinOneDistri):
                  if ''.join(obstype)[:6] == 'gnssro':
                    n_stationId_bak = len(np.unique(stationId_baknc[~np.isnan(stationId_baknc)]))
                    n_stationId_ana = len(np.unique(stationId_ananc[~np.isnan(stationId_ananc)]))
                    n_record_bak   = len(np.unique(recordNum_baknc[~np.isnan(recordNum_baknc)]))
                    n_record_ana   = len(np.unique(recordNum_ananc[~np.isnan(recordNum_ananc)]))
                  else:
                    if 'nan' in stationId_baknc:
                        n_stationId_bak = len(set(stationId_baknc)) -1
                    else:
                        n_stationId_bak = len(set(stationId_baknc))
                    if 'nan' in stationId_ananc:
                        n_stationId_ana = len(set(stationId_ananc)) -1
                    else:
                        n_stationId_ana = len(set(stationId_ananc))

                # plot oma, omb from all vertical levels
                  if ''.join(obstype)[:6] == 'gnssro':
                    basic_plot_functions.plotDistri(latnc,lonnc,ombnc,obstype,varname,vardict[varname][0],expt_obs,n_record_bak,"omb_allLevels")
                    basic_plot_functions.plotDistri(latnc,lonnc,omanc,obstype,varname,vardict[varname][0],expt_obs,n_record_ana,"oma_allLevels")
                  else:
                    basic_plot_functions.plotDistri(latnc,lonnc,ombnc,obstype,varname,vardict[varname][0],expt_obs,n_stationId_bak,"omb_allLevels")
                    basic_plot_functions.plotDistri(latnc,lonnc,omanc,obstype,varname,vardict[varname][0],expt_obs,n_stationId_ana,"oma_allLevels")

            if ''.join(obstype) in profile_group:
                RMSEombs = []
                RMSEomas = []
                obsnums  = []
                ombnums  = []
                omanums  = []
                latncs   = []
                lonncs   = []
                Meanobserrors = []
                Meanerrstarts = []
                obserrornums = []
                errstartnums = []
                if ''.join(obstype)[:6] == 'gnssro':
                    bins = list(np.arange(50000.0, -1000, -2000.))
                    binsfory= list(np.arange(49500.0,-500.0,-2000.))
                else:
                    bins =   [1050.,950.,850.,750.,650.,550.,450.,350.,250.,150.,50.,0.]
                    binsfory=  [1000.,900.,800.,700.,600.,500.,400.,300.,200.,100.,0]
                bins = np.asarray(bins)

                ombnum = len(ombnc)-np.isnan(ombnc).sum()
                ombnums = np.append(ombnums,ombnum)

                omanum = len(omanc)-np.isnan(omanc).sum()
                omanums = np.append(omanums,omanum)

                for j in range(0, len(bins)-1):
                    #print('jban check segmentation fault',j)
                    obsncbin = deepcopy(obsnc)
                    ombncbin = deepcopy(ombnc)
                    omancbin = deepcopy(omanc)
                    latncbin = deepcopy(latnc)
                    lonncbin = deepcopy(lonnc)
                    obserrorbin = deepcopy(obserrornc)
                    errstartbin = deepcopy(errstartnc)
                    errendbin = deepcopy(errendnc)
                    stationbin_bak = deepcopy(stationId_baknc)
                    stationbin_ana = deepcopy(stationId_ananc)
                    recordbin_bak = deepcopy(recordNum_baknc)
                    recordbin_ana = deepcopy(recordNum_ananc)

                    obsncbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    ombncbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    omancbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    latncbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    lonncbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    stationbin_bak[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    stationbin_ana[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    obserrorbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    errstartbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    errendbin[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                    RMSEomb = np.sqrt(np.nanmean(ombncbin**2))
                    RMSEoma = np.sqrt(np.nanmean(omancbin**2))
                    Meanobserror = np.nanmean(obserrorbin)
                    Meanerrstart = np.nanmean(errstartbin)

                    obsnum = len(obsncbin)-np.isnan(obsncbin).sum()
                    obsnums = np.append(obsnums,obsnum)

                    ombnum = len(ombncbin)-np.isnan(ombncbin).sum()
                    ombnums = np.append(ombnums,ombnum)

                    omanum = len(omancbin)-np.isnan(omancbin).sum()
                    omanums = np.append(omanums,omanum)

                    obserrornum = len(obserrorbin)-np.isnan(obserrorbin).sum()
                    obserrornums = np.append(obserrornums,obserrornum)

                    errstartnum = len(errstartbin)-np.isnan(errstartbin).sum()
                    errstartnums = np.append(errstartnums,errstartnum)

                    RMSEombs = np.append(RMSEombs,RMSEomb)
                    RMSEomas = np.append(RMSEomas,RMSEoma)
                    Meanobserrors = np.append(Meanobserrors,Meanobserror)
                    Meanerrstarts = np.append(Meanerrstarts,Meanerrstart)
                    print(obstype)
                    # plot oma, omb from every bin range
                    if (plot_eachLevelDistri):
                      if ''.join(obstype)[:6] == 'gnssro':
                        recordbin_bak[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                        recordbin_ana[np.logical_or(prenc <bins[j+1], prenc >=bins[j])] = np.NaN
                        n_record_bak = len(np.unique(recordbin_bak[~np.isnan(recordbin_bak)]))
                        n_record_ana = len(np.unique(recordbin_ana[~np.isnan(recordbin_ana)]))
                        n_stationId_bak = len(np.unique(stationbin_bak[~np.isnan(stationbin_bak)]))
                        n_stationId_ana = len(np.unique(stationbin_ana[~np.isnan(stationbin_ana)]))
                        print('check every level n_record_bak, n_record_ana=',n_record_bak, n_record_ana)
                      else:
                        if 'nan' in stationbin_bak:
                            n_stationId_bak = len(set(stationbin_bak)) -1
                        else:
                            n_stationId_bak = len(set(stationbin_bak))
                        if 'nan' in stationbin_ana:
                            n_stationId_ana = len(set(stationbin_ana)) -1
                        else:
                            n_stationId_ana = len(set(stationbin_ana))
                      if ''.join(obstype)[:6] == 'gnssro':
                        basic_plot_functions.plotDistri(latncbin,lonncbin,ombncbin,obstype,varname,vardict[varname][0],expt_obs,n_record_bak,"omb_vertbin"+str(j))
                        basic_plot_functions.plotDistri(latncbin,lonncbin,omancbin,obstype,varname,vardict[varname][0],expt_obs,n_record_ana,"oma_vertbin"+str(j))
                      else:
                        basic_plot_functions.plotDistri(latncbin,lonncbin,ombncbin,obstype,varname,vardict[varname][0],expt_obs,n_stationId_bak,"omb_vertbin"+str(j))
                        basic_plot_functions.plotDistri(latncbin,lonncbin,omancbin,obstype,varname,vardict[varname][0],expt_obs,n_stationId_ana,"oma_vertbin"+str(j))

                plotrmsepro(RMSEombs,RMSEomas,binsfory,ombnums,omanums,expt_obs,varname,print_fmt)
                ploterrpro(Meanobserrors,Meanerrstarts,binsfory,obserrornums,errstartnums,expt_obs,varname,print_fmt,"Mean")

            elif ''.join(obstype) in radiance_group:
                obsnc = obsnc.reshape(len(latnc),len(chlist))
                bkgnc = bkgnc.reshape(len(latnc),len(chlist))
                ananc = ananc.reshape(len(latnc),len(chlist))
                ombnc = ombnc.reshape(len(latnc),len(chlist))
                omanc = omanc.reshape(len(latnc),len(chlist))
                # Generate scatter plots
                # Maximum number of variables/channels per figure
                maxsubplts = 16
                nvars = len(channels)
                for channel in channels:
                  ifig = channels.index(channel) + 1
                  varval = vardict.get(varname,['',varname])
                  units = vu.varDictObs[var][0]
                  shortname = '_'+str(channel)
                  shortname = varval[1] +  shortname
                  ivar = channels.index(channel)
                  if ivar == 0:
                    # scatter_verification yields 2 figure types
                    nfigtypes = 2
                    nx_subplt, ny_subplt, figs, subplt_cnt = \
                       init_subplts(channels, nfigtypes, maxsubplts)

                  subplt_cnt = \
                  scatter_verification(ifig, shortname, units, ivar, nvars, \
                                     maxsubplts, subplt_cnt, \
                                     obsnc[:,channel-1], ombnc[:,channel-1], omanc[:,channel-1], \
                                     nx_subplt, ny_subplt, \
                                     nfigtypes, figs, expt_obs,print_fmt)

                # Horizontal distribution of radiance OBS, BCKG, ANA, OMB, OMA
                varval = vardict.get(varname,['',varname])
                units = vu.varDictObs[var][0]
                dotsize = 5.0
                if ''.join(obstype) in radiance_group:
                    color = "BT"
                else:
                    color = "gist_ncar"
                for channel in channels:
                    shortname = '_'+str(channel)
                    shortname = varval[1] + shortname
                    basic_plot_functions.plotDistri(latnc,lonnc,obsnc[:,channel-1], \
                                                obstype,shortname,units,expt_obs,0,"obs", \
                                                None,None,dotsize,color)
                    basic_plot_functions.plotDistri(latnc,lonnc,bkgnc[:,channel-1], \
                                                obstype,shortname,units,expt_obs,0,"bkg", \
                                                None,None,dotsize,color)
                    basic_plot_functions.plotDistri(latnc,lonnc,ananc[:,channel-1], \
                                                obstype,shortname,units,expt_obs,0,"ana", \
                                                None,None,dotsize,color)
                    dmin = -30
                    dmax = 30
                    color = "hsv"
                    basic_plot_functions.plotDistri(latnc,lonnc,ombnc[:,channel-1], \
                                                obstype,shortname,units,expt_obs,0,"omb", \
                                                dmin,dmax,dotsize,color)
                    basic_plot_functions.plotDistri(latnc,lonnc,omanc[:,channel-1], \
                                                obstype,shortname,units,expt_obs,0,"oma", \
                                                dmin,dmax,dotsize,color)

def plotrmsepro(var1,var2,binsfory,ombnums,omanums,EXP_NAME,VAR_NAME,fmt):
    fig, ax1 = plt.subplots()
#   reverse left y-axis
    plt.gca().invert_yaxis()
    plt.grid(True)
    ax1.plot(var1,binsfory,'g-*',markersize=5)
    ax1.plot(var2,binsfory,'r-*',markersize=5)
    print(var1)
    print(np.nanmax(var1))
    if VAR_NAME in 'specific_humidity':
        ax1.set_xlim([0,np.nanmax(var1)])
    else:
        ax1.set_xlim([0,math.ceil(np.nanmax(var1))])

    ax2 = ax1.twinx()
    ax2.spines['right'].set_position(('axes', 1.0))
    ax2.set_yticks(binsfory)
    ax2.set_ylabel('Observation Number(EffectiveQCbg=0)',fontsize=15)

    ax3 = ax1.twinx()
    ax3.set_yticks(binsfory)
    ax3.spines['right'].set_position(('axes', 1.2))
    ax3.set_ylabel('Observation Number(EffectiveQCan=0)',fontsize=15)

    if ''.join(EXP_NAME.split("_")[-1:])[:6] == 'gnssro':
        ax1.set_ylim([0,49500])
        ax1.set_ylabel('Altitude (m)',fontsize=15)
        ax1.set_xlabel(VAR_NAME+'  RMSE of OMB/O & OMA/O '+ vardict[VAR_NAME][0],fontsize=15)
        ax1.legend(('(OMB/O)*100%','(OMA/O)*100%'), loc='upper right',fontsize=15)
        ax2.set_yticklabels(ombnums.astype(np.int))
        ax3.set_yticklabels(omanums.astype(np.int))
    else:
        ax1.set_ylim([1000,0])
        major_ticks = np.arange(0, 1000, 100)
        ax1.set_yticks(major_ticks)
        ax1.set_ylabel('Pressure (hPa)',fontsize=15)
        ax1.set_xlabel(VAR_NAME+'  RMSE  '+ vardict[VAR_NAME][0],fontsize=15)
        ax1.legend(('OMB','OMA'), loc='lower left',fontsize=15)
        ax2.set_yticklabels(reversed(ombnums.astype(np.int)))
        ax3.set_yticklabels(reversed(omanums.astype(np.int)))
    fname = 'RMSE_%s_%s.'%(EXP_NAME,VAR_NAME)+fmt
    print('Saving figure to '+fname)
    plt.savefig(fname,dpi=200,bbox_inches='tight')
    plt.close()

def ploterrpro(var1,var2,binsfory,ombnums,omanums,EXP_NAME,VAR_NAME,fmt,metrics):
    fig, ax1 = plt.subplots()
#   reverse left y-axis
    plt.gca().invert_yaxis()
    plt.grid(True)
    ax1.plot(var1,binsfory,'g-*',markersize=5)
    ax1.plot(var2,binsfory,'r-*',markersize=5)
    print(np.nanmax(var1))
    if VAR_NAME in 'specific_humidity':
        ax1.set_xlim([0,np.nanmax(var1)])
    else:
        ax1.set_xlim([0,math.ceil(np.nanmax(var2))])

    ax2 = ax1.twinx()
    ax2.spines['right'].set_position(('axes', 1.0))
    ax2.set_yticks(binsfory)
    #ax2.set_ylabel('Obs Num for obserror(PreQC>-10 and <3)',fontsize=8)
    ax2.set_ylabel('Obs Num for obserror (EffectiveQCbg=0)',fontsize=8)

    ax3 = ax1.twinx()
    ax3.set_yticks(binsfory)
    ax3.spines['right'].set_position(('axes', 1.2))
    ax3.set_ylabel('Obs Num for EffectiveError0(EffectiveQCbg=0)',fontsize=8)

    if ''.join(EXP_NAME.split("_")[-1:])[:6] == 'gnssro':
        ax1.set_ylim([0,49500])
        ax1.set_ylabel('Altitude (m)',fontsize=15)
        #ax1.set_xlabel(VAR_NAME+'  RMSE of OMB/O & OMA/O '+ vardict[VAR_NAME][0],fontsize=15)
        #ax1.legend(('(OMB/O)*100%','(OMA/O)*100%'), loc='upper right',fontsize=15)
        ax2.set_yticklabels(ombnums.astype(np.int))
        ax3.set_yticklabels(omanums.astype(np.int))
    else:
        ax1.set_ylim([1000,0])
        major_ticks = np.arange(0, 1000, 100)
        ax1.set_yticks(major_ticks)
        ax1.set_ylabel('Pressure (hPa)',fontsize=15)
        #ax1.set_xlabel(VAR_NAME+' '+ metrics +' '+ vardict[VAR_NAME][0],fontsize=15)
        #ax1.legend(('obserr','EffectiveError0'), loc='lower left',fontsize=15)
        ax2.set_yticklabels(reversed(ombnums.astype(np.int)))
        ax3.set_yticklabels(reversed(omanums.astype(np.int)))
    ax1.set_xlabel(VAR_NAME+' '+ metrics +' '+ vardict[VAR_NAME][0],fontsize=15)
    ax1.legend(('obserr','EffectiveError0'), loc='lower left',fontsize=15)

    fname = metrics+'_%s_%s.'%(EXP_NAME,VAR_NAME)+fmt
    #fname = 'RMSE_%s_%s.'%(EXP_NAME,VAR_NAME)+fmt
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
    nnfigs = (int(len(subpltlist) / maxsubplts) + 1)
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
    jfig = int((ifig-1)/maxsubplts)
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
            xlab = 'h(xb) - y'
            ylab = 'h(xa) - y'
        else:
            print('WARNING: scatter_verification has no definitions for nfigtypes == ',nfigtypes)
            continue

        # Uncomment these 2 lines to put x/y labels only on peripheral subplts
        #if kfig <= min(nvars,numsubplts) - nx_subplt[0] : xlab = ''
        #if np.mod(kfig,nx_subplt[0]) != 1 : ylab = ''
        ax = figs[int(offset+jfig)].add_subplot(ny_subplt[0],nx_subplt[0],kfig)

        if iifig == 0:
            #Add scatter plot for h(x) vs. y
            fname = 'XB_XA_%s'%(EXP_NAME)
            stat = scatter_one2ones( \
                obs, [obs-omb , obs-oma], \
                ['x_b' , 'x_a'], True, xlab, ylab, \
                varname, varunits, ax)
        if iifig == 1:
            #Add scatter plot for OMA vs. OMB
            fname = 'OMB_OMA_%s'%(EXP_NAME)
            stat = scatter_one2ones( \
                -omb , [-oma], \
                [], False, xlab, ylab, \
                varname, varunits, ax)

        if stat != 0:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(0.5, 0.5, '[NO DATA]',
                {'color': 'k', 'fontsize': 12}, 
                ha='center', va='center', transform=ax.transAxes)

        if nnfigs > 1: fname=fname+'_%d-of-%d'%(jfig,nnfigs)
        fname=fname+'.'+fmt

        if (ivar == nvars-1 or subplt_cnt[jfig] == numsubplts):
            #Save the figure to file
            print('Saving figure to '+fname)
            figs[int(offset+jfig)].subplots_adjust(wspace=0.35,hspace=0.35)
            figs[int(offset+jfig)].savefig(fname,dpi=200,bbox_inches='tight')

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
    fsize_leg = 6.0
    fsize_lab = 4.5

    ax.text(0.03, 0.97 - len(LEG) * 0.125, VAR_NAME,
        {'color': 'k', 'fontsize': fsize_leg}, 
        ha='left', va='top', transform=ax.transAxes)

    if len(XVAL) == 0: 
        print('WARNING in scatter_one2ones: len(XVAL)==0; skipping this dataset')
        return 1
    NVALS = np.asarray([])
    for i,YVAL in enumerate(YVALS):
        if len(XVAL) != len(YVAL):
            print('ERROR: Incorrect usage of scatter_one2ones, YVALS must be list of arrays.')
            os._exit()
        not_nan = np.isfinite(XVAL) & np.isfinite(YVAL)
        NVALS = np.append(NVALS,np.sum(not_nan))

    if np.all(NVALS == 0):
        print('WARNING in scatter_one2ones: all(XVAL/YVAL) are non-finite; skipping this dataset')
        return 1

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
        col = colors[np.mod(i,len(colors))]
        mind = np.mod(i,len(markers))
        ax.plot( XVAL, YVAL, markers[mind], color = col, \
                 markersize = msizes[mind], alpha=0.5 )

    if XLAB != '':
        label = XLAB
        if UNITS != '': label = label+' '+UNITS
        ax.set_xlabel(label,fontsize=6)
    if YLAB != '':
        label = YLAB
        if UNITS != '': label = label+' '+UNITS
        ax.set_ylabel(label,fontsize=6)
    if len(LEG) > 0:
       ax.legend(LEG, loc='upper left',fontsize=5)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    xymin=min(xmin,ymin)
    xymax=max(xmax,ymax)

    #Could adjust limits/ticks for better aesthetics
    #round_nmbr = 5
    #xymin=np.floor(min(xmin,ymin) / round_nmbr) * round_nmbr
    #xymax=np.ceil(max(xmax,ymax) / round_nmbr) * round_nmbr

    # Add linear regression and statistics
    predictor = np.arange(xymin, xymax, (xymax - xymin) / 10.0)
    tx = 0.98
    ty = 0.02
    nline = ''
    for j,YVAL in enumerate(reversed(YVALS)):
        i = len(YVALS) - j - 1
        if NVALS[j]==0: continue
        del not_nan
        not_nan = np.isfinite(XVAL) & np.isfinite(YVAL)

        # Add linear fit to YVAL vs. XVAL
        p = np.polyfit(XVAL[not_nan],YVAL[not_nan],1)
        predicted = p[0] * predictor + p[1]
        col0 = colors[np.mod(i,len(colors))]
#        bright = 0.5
#        col = bright * np.asarray([1.,1.,1.]) + (1. - bright) * np.asarray(col0)
        dimmer = 0.35
        col = (1. - dimmer) * np.asarray(col0)
        ax.plot( predictor, predicted, '-', color = col, lw=1.2 )

        # Add statistics for YVAL vs. XVAL
        stat = 'N = %d\nslope: %0.2f  '%(NVALS[j],p[0])
        if show_stats:
            RMSE = np.sqrt( np.sum( np.square(YVAL[not_nan] - XVAL[not_nan]) ) / NVALS[j] )
            BIAS = np.sum( YVAL[not_nan] - XVAL[not_nan] ) / NVALS[j]
            stat = stat+'\nRMSE: %0.2f \nBIAS: %0.2f'%(RMSE,BIAS)
        ax.text(tx, ty, stat+nline, \
            {'color': col, 'fontsize': fsize_lab}, 
            ha='right', va='bottom', backgroundcolor=[1,1,1,0.2], \
            clip_on=True, transform=ax.transAxes)
        nline = nline + ''.join('\n' * stat.count('\n')) + '\n'

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
    return 0

def main():
    readdata()

if __name__ == '__main__': main()
