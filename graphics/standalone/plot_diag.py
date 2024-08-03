import os
import sys
sys.path.insert(1, '../')
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
import netCDF4 as nc4
import config as conf
import var_utils as vu
import JediDB

'''
Directory Structure for ctest:
test/
 ├── Data/os/
 │        ├── obsout_3denvar_bump_sondes_0000.nc4
 │        ├── obsout_3denvar_bump_sondes_0001.nc4
 │        ├── ...
 ├── graphics/
 │   ├── basic_plot_functions.py
 │   ├── ...
 │   ├── standalone
 │   │   ├── plot_diag.py
 │   │   ├── ...

Directory Structure for cycling:
test/
 ├── dbOut/
 │   ├── obsout_3denvar_bump_sondes_0000.h5
 │   ├── obsout_3denvar_bump_sondes_0001.h5
 │   ├── ...
 ├── graphics/
 │   ├── plot_diag.py
 │   ├── basic_plot_functions.py
 │   ├── ...
'''

# need add Ps
    #Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
    #Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bending_angle unit: rad

def readdata():

  imageFmt = 'png' #lower fidelity, faster
  #imageFmt = 'pdf' #higher fidelity, slower

  # Diagnostic omb,oma and hofx files located in diagdir
  # assume obsout_hofx_*.h5 or obsout_da_*.h5 are located '../../' of current standalone folder
  diagdir    = '../../'
  diagprefix = 'obsout_'
  #diagsuffix = '_*.nc4'  #for ctests
  diagsuffix = '_*.h5'   #for cycling
  # * in diagsuffix is 4-digit processor rank [0000 to XXXX]
  #npedigits = 4

  # controls for generating 2D maps of locations
  makeDistributionPlots = True
  plot_allinOneDistri = True   # plot profileObsTypes (includes all levels) and sfcObsTypes.
  plot_eachLevelDistri = False # plot every level separately for profileObsTypes.

  # NOUTER: number of outer iterations
  # can be set as env variable 'NOUTER'
  # irrelevant for hofx applications
  NOUTER=str(os.getenv('NOUTER', 2))

  # observation types to analyze
  profileObsTypes = [
    'sondes',
    'aircraft',
    'satwind',
    'satwnd',
    'gnssroref',
    'gnssrobndropp1d',
  ]
  sfcObsTypes = [
    'sfc',
  ]
  radianceObsTypes = [
    'abi_g16',
    'ahi_himawari8',
    'amsua_aqua',
    'amsua_metop-a',
    'amsua_metop-b',
    'amsua_n15',
    'amsua_n18',
    'amsua_n19',
    'amsua_n19--hydro',
    'amsua_n19--nohydro',
    'amsua-cld_aqua',
    'amsua-cld_metop-a',
    'amsua-cld_metop-b',
    'amsua-cld_n15',
    'amsua-cld_n18',
    'amsua-cld_n19',
    'mhs_metop-a',
    'mhs_metop-b',
    'mhs_n19',
    'mhs_n18',
  ]

  analyzedObsTypeGroups = []
  analyzedObsTypeGroups.append(profileObsTypes)
  analyzedObsTypeGroups.append(sfcObsTypes)
  analyzedObsTypeGroups.append(radianceObsTypes)

  # application-dependent ObsGroups for variational and hofx applications
  variationalApp = 'variational'
  hofxApp = 'hofx'
  ApplicationObsGroups = {
    variationalApp: {
      'obsGroup': 'ObsValue',
      'ombgGroup': 'ombg',
      'omanGroup': 'oman',
      'qcbgGroup': 'EffectiveQC0',
      'qcanGroup': 'EffectiveQC'+NOUTER,
      'obserrGroup': 'ObsError',
      'errstartGroup': 'EffectiveError0',
    },
    hofxApp: {
      'obsGroup': 'ObsValue',
      'hofxGroup': 'hofx',
      'qcbgGroup': 'EffectiveQC',
      'obserrGroup': 'ObsError',
      'errstartGroup': 'EffectiveError',
    },
  }

  # fixed MetaData variable names
  metaGroup = 'MetaData'
  latitude = metaGroup+'/latitude'
  longitude = metaGroup+'/longitude'
  pressure = metaGroup+'/pressure'
  altitude = metaGroup+'/height'
  occultID = metaGroup+'/occulting_sat_is'
  satelliteId = metaGroup+'/satelliteIdentifier'
  stationID = metaGroup+'/stationIdentification'
  recordNumber = metaGroup+'/sequenceNumber'

  # settings for binning coordinates
  binningCoordinates = {
    altitude: {'varName': 'Altitude (m)'},
    pressure: {'varName': 'Pressure (hPa)'},
  }
  binningCoordinates[altitude]['edges'] = np.arange(50000.0, -1000., -2000.)
  binningCoordinates[pressure]['edges'] = \
    np.array([1050., 950., 850., 750., 650., 550., 450., 350., 250., 150., 50., 0.])

  # collect all file names that fit file name format
  obsoutfiles = []
  for file in os.listdir(diagdir):
    if fnmatch.fnmatch(file, diagprefix+'*'+diagsuffix):   # 1tile
      obsoutfiles.append(diagdir+file)

  # Group files by experiment-obstype combination
  #  (e.g., 3dvar_aircraft), where each group
  #  contains files from all PE's
  exob_groups = [[]]
  for j, file in enumerate(obsoutfiles):
    if JediDB.IODAFileIsRanked(file):
      # exob_group_name excludes everything outside the first/final '_'
      exob_group_name =  '_'.join(file.split("_")[1:][:-1])
    else:
      # Remove suffix .nc or .h5 or .nc4
      remove_suffix_name = '.'.join(file.split(".")[:-1])
      # exob_group_name excludes everything outside the first '_'
      exob_group_name = '_'.join(remove_suffix_name.split("_")[1:])
    for i, exob_group in enumerate(exob_groups):
      if exob_group_name in exob_group:
        # If exob_group with exob_group_name exists, add new file to it
        update_group = exob_group
        update_group.append(file)
        exob_groups[i][:] = update_group
        break
      elif i == len(exob_groups)-1:
        # If group with exob_group_name does not exist, add one
        new_group = [exob_group_name]
        new_group.append(file)
        exob_groups.append(new_group)
        break
  exob_groups = exob_groups[1:][:]

  # Loop over unique experiment-obstype groups
  for exob_group in exob_groups:
    expt_obs = exob_group[0]
    print("Processing ", expt_obs)

    # Determine obstype from expt_obstype string
    expt_parts = expt_obs.split("_")
    nstr = len(expt_parts)
    obstype = 'none'
    for i in range(0, nstr):
      obstype_ = '_'.join(expt_parts[i:nstr])
      for obsTypeGroup in analyzedObsTypeGroups:
        if obstype_ in obsTypeGroup:
          obstype = obstype_

    if obstype == 'none':
      print('obstype not selected, skipping data: '+expt_obs)
      continue

    # Sort files based on PE
    obsFiles = np.array(deepcopy(exob_group[1:]))
    if JediDB.IODAFileIsRanked(obsFiles[0]):
      PEs = []
      for file in obsFiles:
        print('check file=',file)
        PEs.append(int(file.split('_')[-1].split('.')[0]))
      obsFiles = obsFiles[np.argsort(np.array(PEs))]

    # Determine total nlocs
    nlocs = 0
    for file in obsFiles:
      ncDB = nc4.Dataset(file, 'r')
      ncDB.set_auto_mask(False)
      nlocs += ncDB.dimensions['Location'].size
      ncDB.close()

    # get some basic info from obsFiles[0]
    ncDB = nc4.Dataset(obsFiles[0], 'r')
    ncDB.set_auto_mask(False)

    # define a channel list, if available
    if 'Channel' in ncDB.dimensions:
      nchans = ncDB.variables['Channel'][:]
      nch = len(nchans)

    # applicationType: type of jedi application (hofx, variational)
    # assume hofx at first, then check for presence of 'ombgGroup'
    applicationType = hofxApp
    ncVarList = []
    for group in ncDB.groups:
      if group == ApplicationObsGroups[variationalApp]['ombgGroup']:
        applicationType = variationalApp
      for var in ncDB.groups[group].variables:
        ncVarList+= [group+'/'+var]

    ncDB.close()

    obsGroup = ApplicationObsGroups[applicationType].get('obsGroup', None)
    ombgGroup = ApplicationObsGroups[applicationType].get('ombgGroup', None)
    omanGroup = ApplicationObsGroups[applicationType].get('omanGroup', None)
    hofxGroup = ApplicationObsGroups[applicationType].get('hofxGroup', None)
    qcbgGroup = ApplicationObsGroups[applicationType].get('qcbgGroup', None)
    qcanGroup = ApplicationObsGroups[applicationType].get('qcanGroup', None)
    obserrGroup = ApplicationObsGroups[applicationType].get('obserrGroup', None)
    errstartGroup = ApplicationObsGroups[applicationType].get('errstartGroup', None)

    # select all variables that were simulated
    # only simulated variables are relevant to OMB & OMA diagnostics
    simulatedGroup = None
    for g in [ombgGroup, hofxGroup, omanGroup, errstartGroup]:
      if g is not None:
        simulatedGroup = g
        simulatedVariables = [
          ''.join(var.split("/")[1:])
          for var in ncVarList
          if (g+'/' == var[:len(g)+1])
        ]
        break

    assert simulatedGroup is not None, 'simulatedGroup was not selected'

    # establish location variables to read from files
    coordVars = [latitude, longitude]

    # observation-type-specific coordVars
    binCoord = None
    station = None
    record = None
    isGNSSRO = 'gnssro' in obstype
    if isGNSSRO:
      binCoord = altitude
      if metaGroup+'/occulting_sat_is' not in ncVarList:
        station = satelliteId
      else:
        station = occultID
      record = recordNumber
      coordVars += [binCoord, station, record]
    elif obstype in profileObsTypes:
      binCoord = pressure
      station = stationID
      coordVars += [binCoord, station]
    elif obstype in sfcObsTypes:
      station = stationID
      coordVars += [station]

    # loop over simulated variables
    for ivar, varName in enumerate(simulatedVariables):
      print("Working on ", varName)

      # assume obs, obserr, qcbg, and errstart ObsGroups are always present
      obs = obsGroup+'/'+varName
      obserror = obserrGroup+'/'+varName
      qcb = qcbgGroup+'/'+varName
      errstart = errstartGroup+'/'+varName

      obsVars = [obs, obserror, qcb, errstart]

      # read ObsGroups for which presence is application-dependent
      ombg = None
      if ombgGroup is not None:
        ombg = ombgGroup+'/'+varName
        obsVars += [ombg]

      oman = None
      if omanGroup is not None:
        oman = omanGroup+'/'+varName
        obsVars += [oman]

      qca = None
      if qcanGroup is not None:
        qca = qcanGroup+'/'+varName
        obsVars += [qca]

      hofx = None
      if hofxGroup is not None:
        hofx = hofxGroup+'/'+varName
        obsVars += [hofx]

      # generate list of required variables to read
      readVars = coordVars + obsVars

      # Build up in-memory database in loops over obsFiles and readVars
      ss = 0
      db = {}
      for file in obsFiles:
        ncDB = nc4.Dataset(file, 'r')
        ncDB.set_auto_mask(False)
        nn = ncDB.dimensions['Location'].size
        ee = ss + nn

        for varGrp in readVars:
          if varGrp is None: continue
          var, grp = vu.splitObsVarGrp(varGrp)
          # recordNumber is not available in ctest data files
          # In that case, assume record and station counts are identical
          if (varGrp == record and
            (var not in ncDB.groups[grp].variables)):
            if record in db: del db[record]
            if record in readVars: readVars[readVars.index(record)] = None
            record = None
            continue

          assert var in ncDB.groups[grp].variables, (
            varGrp,' not in ncDB.groups[grp].variables: ',ncDB.groups[grp].variables)

          varHasNChansDim = 'Channel' in ncDB.groups[grp].variables[var].dimensions

          # initialize each varGrp
          if varGrp not in db:
            dtype = ncDB.groups[grp].variables[var].dtype
            if 'str' in str(dtype) or 'byte' in str(dtype): dtype = np.object_
            if varHasNChansDim:
              db[varGrp] = np.empty((nlocs, nch), dtype=dtype)
            else:
              db[varGrp] = np.empty(nlocs, dtype=dtype)

          # read
          try:
            if varHasNChansDim:
              db[varGrp][ss:ee,:] = ncDB.groups[grp].variables[var][:,:]
            else:
              db[varGrp][ss:ee] = ncDB.groups[grp].variables[var][:]
          except UnicodeDecodeError as p:
            if varGrp == stationID:
              # stationID is mangled in ioda-v1 files, sometimes converted to
              # ioda-v2 files; handle UnicodeDecodeError exception
              if stationID in db: del db[stationID]
              if stationID in readVars: readVars[readVars.index(stationID)] = None
              if station is not None and stationID == station: station = None
              pass
            else:
              print('Incorrect unicode format for '+varGrp+' in '+file)
              raise UnicodeDecodeError(p)

        ss = ee
        ncDB.close()

      # create equivalent qcan for hofx application
      if qca is None and qcb is not None:
        qca = 'qcan/'+varName
        db[qca] = deepcopy(db[qcb])

      # only keep data @EffectiveQC==0
      failbg = db[qcb] != 0
      failan = db[qca] != 0

      # convert pressure to hPa
      if pressure in db:
        if np.max(db[pressure]) > 10000.0:
          db[pressure] = np.divide(db[pressure], 100.0)

      # background qc checks
      for var in [obs, ombg, hofx, obserror, errstart]:
        if var is None: continue
        db[var][np.greater(np.abs(db[var]), 1.0e+15)] = np.NaN
        db[var][failbg] = np.NaN

      # analysis QC checks
      for var in [oman]:
        if var is None: continue
        db[var][np.greater(np.abs(db[var]), 1.0e+15)] = np.NaN
        db[var][failan] = np.NaN

      # create equivalent bg/an departures for hofx application
      if hofx is not None:
        if ombg is None:
          ombg = 'ombg/'+varName
          db[ombg] = db[obs] - db[hofx]

        if oman is None:
          oman = 'oman/'+varName
          db[oman] = db[obs] - db[hofx]

      # convert departures to omb/oma and bkg/ana
      omb = 'omb/'+varName
      db[omb] = db[ombg]
      bkg = 'bkg/'+varName
      db[bkg] = db[obs] - db[omb]
      passbg = np.logical_and(~failbg, np.isfinite(db[omb]))
      for var in [obs, ombg, obserror, errstart, omb, bkg]:
        db[var][~passbg] = np.NaN

      oma = 'oma/'+varName
      db[oma] = db[oman]
      ana = 'ana/'+varName
      db[ana] = db[obs] - db[oma]
      passan = np.logical_and(~failan, np.isfinite(db[oma]))
      for var in [oman, oma, ana]:
        db[var][~passan] = np.NaN

      if not np.isfinite(db[omb]).any():
        print('all values are NaN of inf: ', varName, obstype)
        continue

      # diagnose relative omb/oma for GNSSRO, which varies
      # across large orders of magnitude
      if isGNSSRO:
        db[omb] = (db[omb]/db[obs])*100.
        db[oma] = (db[oma]/db[obs])*100.

      # plot oma, omb from all vertical levels
      if makeDistributionPlots and plot_allinOneDistri and obstype not in radianceObsTypes:

        nProfile_bak = 0
        nProfile_ana = 0

        if station is not None:
          if 'int' in str(db[station].dtype):
            nProfile_bak = len(np.unique(db[station][passbg]))
            nProfile_ana = len(np.unique(db[station][passan]))
          else:
            nProfile_bak = len(set(db[station][passbg]))
            if 'nan' in db[station][passbg]:
                nProfile_bak -= 1

            nProfile_ana = len(set(db[station][passan]))
            if 'nan' in db[station][passan]:
                nProfile_ana -= 1

        # record overrides station when available
        if record is not None:
          nProfile_bak = len(np.unique(db[record][passbg]))
          nProfile_ana = len(np.unique(db[record][passan]))

        basic_plot_functions.plotDistri(db[latitude], db[longitude], db[omb], obstype, varName, vu.varDictObs[varName][0], expt_obs, nProfile_bak, "omb_allLevels")
        basic_plot_functions.plotDistri(db[latitude], db[longitude], db[oma], obstype, varName, vu.varDictObs[varName][0], expt_obs, nProfile_ana, "oma_allLevels")

      if binCoord is not None and binCoord in db:
        binVar = binningCoordinates[binCoord]['varName']
        binEdges = binningCoordinates[binCoord]['edges']
        binLeft = binEdges[:-1]
        binRight = binEdges[1:]
        binCenters = 0.5*(binLeft+binRight)

        nBins = len(binEdges)-1

        binnedVars = [obs, omb, oma, obserror, errstart]
        Count = {}
        RMS = {}
        dbBinned = {}
        for var in binnedVars:
          Count[var] = np.empty(nBins, dtype=int)
          RMS[var] = np.empty(nBins, dtype=np.float64)

        for iBin, (edge1, edge2) in enumerate(list(zip(binLeft, binRight))):
          minEdge = np.min(np.array([edge1, edge2]))
          maxEdge = np.max(np.array([edge1, edge2]))

          outsideLevel = np.logical_or(db[binCoord] <  minEdge,
                                       db[binCoord] >= maxEdge)

          for var in binnedVars:
            dbBinned[var] = deepcopy(db[var])
            dbBinned[var][outsideLevel] = np.NaN

            Count[var][iBin] = np.isfinite(dbBinned[var]).sum()
            if Count[var][iBin] > 0:
              RMS[var][iBin] = np.sqrt(np.nanmean(np.square(dbBinned[var])))
            else:
              RMS[var][iBin] = np.NaN

          # plot oma, omb from every bin range
          if makeDistributionPlots and plot_eachLevelDistri:

            passbgBin = np.logical_and(~outsideLevel, passbg)
            passanBin = np.logical_and(~outsideLevel, passan)

            nProfile_bak = 0
            nProfile_ana = 0

            if station is not None:
              if 'int' in str(db[station].dtype):
                nProfile_bak = len(np.unique(db[station][passbgBin]))
                nProfile_ana = len(np.unique(db[station][passanBin]))
              else:
                nProfile_bak = len(set(db[station][passbgBin]))
                if 'nan' in db[station][passbgBin]:
                  nProfile_bak -= 1

                nProfile_ana = len(set(db[station][passanBin]))
                if 'nan' in db[station][passanBin]:
                  nProfile_ana -= 1

            # record overrides station when available
            if record is not None:
              nProfile_bak = len(np.unique(db[record][passbgBin]))
              nProfile_ana = len(np.unique(db[record][passanBin]))
              #print('check every level nProfile_bak, nProfile_ana=', nProfile_bak, nProfile_ana)

            basic_plot_functions.plotDistri(
              db[latitude], db[longitude], dbBinned[omb],
              obstype, varName, vu.varDictObs[varName][0], expt_obs,
              nProfile_bak, "omb_vertbin"+str(iBin))
            basic_plot_functions.plotDistri(
              db[latitude], db[longitude], dbBinned[oma],
              obstype, varName, vu.varDictObs[varName][0], expt_obs,
              nProfile_ana, "oma_vertbin"+str(iBin))

        plotprofile(RMS[omb], 'OMB',
                    RMS[oma], 'OMA',
                    binCenters, binVar,
                    Count[omb], Count[oma],
                    expt_obs, varName, imageFmt, "RMS")
        ploterrpro(RMS[obserror], 'file ObsError',
                   RMS[errstart], 'UFO EffectiveError',
                   binCenters, binVar,
                   Count[obserror], Count[errstart],
                   expt_obs, varName, imageFmt, "RMS")

      elif obstype in radianceObsTypes:
        # Generate scatter plots
        # Maximum number of variables/channels per figure
        maxsubplts = 16
        var = varName
        varval = vu.varDictObs.get(varName,['', varName])
        units = vu.varDictObs[var][0]

        # plot all channels that are in both ncDB and config
        ObsSpaceInfo = conf.DiagSpaceConfig.get(obstype, conf.nullDiagSpaceInfo)
        confChans = ObsSpaceInfo.get('channels',[vu.miss_i])
        plotChannels = set(confChans).union(set(nchans))
        nvars = len(plotChannels)
        for ivar, channel in enumerate(plotChannels):
          ich = list(nchans).index(channel)
          shortname = varval[1] + str(channel)
          if ivar == 0:
            # scatter_verification yields 2 figure types
            nfigtypes = 2
            nx_subplt, ny_subplt, figs, subplt_cnt = \
               init_subplts(plotChannels, nfigtypes, maxsubplts)

          subplt_cnt = \
          scatter_verification(
            ivar+1, shortname, units, ivar, nvars,
            maxsubplts, subplt_cnt,
            db[obs][:,ich], db[omb][:,ich], db[oma][:,ich],
            nx_subplt, ny_subplt,
            nfigtypes, figs, expt_obs, imageFmt)

        # Horizontal distribution of radiance OBS, BCKG, ANA, OMB, OMA
        if makeDistributionPlots:
          dotsize = 3.0
          for channel in plotChannels:
            ich = list(nchans).index(channel)
            shortname = varval[1] + str(channel)

            color = "BT"
            basic_plot_functions.plotDistri(db[latitude], db[longitude], db[obs][:,ich],
                                        obstype, shortname, units, expt_obs, 0, "obs",
                                        None, None, dotsize, color)
            basic_plot_functions.plotDistri(db[latitude], db[longitude], db[bkg][:,ich],
                                        obstype, shortname, units, expt_obs, 0, "bkg",
                                        None, None, dotsize, color)
            basic_plot_functions.plotDistri(db[latitude], db[longitude], db[ana][:,ich],
                                        obstype, shortname, units, expt_obs, 0, "ana",
                                        None, None, dotsize, color)
            dmin = -30
            dmax = 30
            color = "hsv"
            basic_plot_functions.plotDistri(db[latitude], db[longitude], db[omb][:,ich],
                                        obstype, shortname, units, expt_obs, 0, "omb",
                                        dmin, dmax, dotsize, color)
            basic_plot_functions.plotDistri(db[latitude], db[longitude], db[oma][:,ich],
                                        obstype, shortname, units, expt_obs, 0, "oma",
                                        dmin, dmax, dotsize, color)

def plotprofile(xVals1, xLabel1,
                xVals2, xLabel2,
                yVals, yLabel,
                counts1, counts2,
                EXP_NAME, varName, fmt, metric):

  fig, ax1 = plt.subplots()
  plt.grid(True)
  ax1.plot(xVals1, yVals,'b-o', markersize=5)
  ax1.plot(xVals2, yVals,'r--*', markersize=5)
  if varName in 'specificHumidity':
    ax1.set_xlim([0, np.nanmax(xVals1)])
  else:
    ax1.set_xlim([0, math.ceil(np.nanmax(xVals1))])

  ax1.set_ylim([min(yVals), max(yVals)])
  ax1.set_ylabel(yLabel, fontsize=15)
  varUnits = '('+vu.varDictObs[varName][0]+')'
  if 'gnssro' in EXP_NAME:
    ax1.set_xlabel(varName+' '+metric+'[(y-h(x))/y] '+varUnits, fontsize=15)
  else:
    ax1.set_xlabel(varName+' '+metric+'[y-h(x)] '+varUnits, fontsize=15)

  ax2 = ax1.twinx()
  ax2.spines['right'].set_position(('axes', 1.0))
  ax2.set_yticks(yVals)
  ax2.set_ylabel('BG Count', fontsize=12)
  ax2.set_yticklabels(counts1.astype(int))

  ax3 = ax1.twinx()
  ax3.set_yticks(yVals)
  ax3.spines['right'].set_position(('axes', 1.2))
  ax3.set_ylabel('AN Count', fontsize=12)
  ax3.set_yticklabels(counts2.astype(int))

  if 'Pressure' in yLabel:
    ax1.set_yticks(yVals)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()

  ax1.legend((xLabel1, xLabel2), loc='best', fontsize=15)

  fname = metric+'ofOMM_%s_%s.'%(EXP_NAME, varName)+fmt
  print('Saving figure to '+fname)
  plt.savefig(fname, dpi=200, bbox_inches='tight')
  plt.close()

def ploterrpro(xVals1, xLabel1,
               xVals2, xLabel2,
               yVals, yLabel,
               counts1, counts2,
               EXP_NAME, varName, fmt, metric):

  fig, ax1 = plt.subplots()
  plt.grid(True)
  ax1.plot(xVals1, yVals,'b-o', markersize=5)
  ax1.plot(xVals2, yVals,'r--*', markersize=5)
  if varName in 'specificHumidity':
    ax1.set_xlim([0, np.nanmax(xVals1)])
  else:
    ax1.set_xlim([0, math.ceil(np.nanmax(xVals2))])
  ax1.set_ylim([min(yVals), max(yVals)])
  ax1.set_ylabel(yLabel, fontsize=15)
  if 'gnssro' in EXP_NAME:
    varUnits = '(unitless)'
  else:
    varUnits = '('+vu.varDictObs[varName][0]+')'
  ax1.set_xlabel(varName+' '+metric+'(σ_o) '+varUnits, fontsize=15)

  ax2 = ax1.twinx()
  ax2.spines['right'].set_position(('axes', 1.0))
  ax2.set_yticks(yVals)
  ax2.set_ylabel('BG Count', fontsize=12)
  ax2.set_yticklabels(counts1.astype(int))

  ax3 = ax1.twinx()
  ax3.set_yticks(yVals)
  ax3.spines['right'].set_position(('axes', 1.2))
  ax3.set_ylabel('AN Count', fontsize=12)
  ax3.set_yticklabels(counts2.astype(int))

  if 'Pressure' in yLabel:
    ax1.set_yticks(yVals)
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()

  ax1.legend((xLabel1, xLabel2), loc='best', fontsize=15)

  fname = metric+'ofObsError_%s_%s.'%(EXP_NAME, varName)+fmt
  print('Saving figure to '+fname)
  plt.savefig(fname, dpi=200, bbox_inches='tight')
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

  nx = int(np.ceil(np.sqrt(nsubplts)))
  ny = int(np.ceil(np.true_divide(nsubplts, nx)))
  nx_subplt.append(nx)
  ny_subplt.append(ny)
  figs = []
  for ifig in range(0, nnfigs*nfigtypes):
    fig = plt.figure()
    inch_size = 1.9
    fig.set_size_inches(nx_subplt[0]*inch_size, ny_subplt[0]*inch_size)
    figs.append(fig)
  subplt_cnt = np.zeros(nnfigs*nfigtypes)

  return nx_subplt, ny_subplt, figs, subplt_cnt

def scatter_verification(ifig, varName, varUnits, ivar, nvars,
                         maxsubplts, subplt_cnt,
                         obs, omb, oma,
                         nx_subplt, ny_subplt,
                         nfigtypes, figs, EXP_NAME, fmt):
#================================================================
#INPUTS:
# ifig       - subplot number
# varName    - variable name
# varUnits   - variable units
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
#PURPOSE: Generate verification subplots for varName amongst
#         subplots containing multiple variables
#================================================================

  nnfigs = len(figs) / nfigtypes
  jfig = int((ifig-1)/maxsubplts)
  kfig = np.mod(ifig-1, maxsubplts)+1
  subplt_cnt[jfig] += 1
  if jfig == nnfigs-1 :
    numsubplts = np.mod(nvars, maxsubplts)
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
      print('WARNING: scatter_verification has no definitions for nfigtypes == ', nfigtypes)
      continue

    # Uncomment these 2 lines to put x/y labels only on peripheral subplts
    #if kfig <= min(nvars, numsubplts) - nx_subplt[0] : xlab = ''
    #if np.mod(kfig, nx_subplt[0]) != 1 : ylab = ''
    ax = figs[int(offset+jfig)].add_subplot(ny_subplt[0], nx_subplt[0], kfig)

    if iifig == 0:
      #Add scatter plot for h(x) vs. y
      fname = 'XB_XA_%s'%(EXP_NAME)
      stat = scatter_one2ones(
          obs, [obs-omb , obs-oma],
          ['x_b' , 'x_a'], True, xlab, ylab,
          varName, varUnits, ax)
    if iifig == 1:
      #Add scatter plot for OMA vs. OMB
      fname = 'OMB_OMA_%s'%(EXP_NAME)
      stat = scatter_one2ones(
          -omb , [-oma],
          [], False, xlab, ylab,
          varName, varUnits, ax)

    if stat != 0:
      ax.set_xticks([])
      ax.set_yticks([])
      ax.text(0.5, 0.5, '[NO DATA]',
              {'color': 'k', 'fontsize': 12},
              ha='center', va='center', transform=ax.transAxes)

    if nnfigs > 1: fname=fname+'_%d-of-%d'%(jfig, nnfigs)
    fname=fname+'.'+fmt

    if (ivar == nvars-1 or subplt_cnt[jfig] == numsubplts):
      #Save the figure to file
      print('Saving figure to '+fname)
      figs[int(offset+jfig)].subplots_adjust(wspace=0.35, hspace=0.35)
      figs[int(offset+jfig)].savefig(fname, dpi=200, bbox_inches='tight')

    offset += nnfigs

  return subplt_cnt

def scatter_one2ones(XVAL, YVALS, LEG, show_stats, XLAB, YLAB, VAR_NAME, UNITS, ax):
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
  for i, YVAL in enumerate(YVALS):
    if len(XVAL) != len(YVAL):
      print('ERROR: Incorrect usage of scatter_one2ones, YVALS must be list of arrays.')
      os._exit()
    not_nan = np.isfinite(XVAL) & np.isfinite(YVAL)
    NVALS = np.append(NVALS, np.sum(not_nan))

  if np.all(NVALS == 0):
    print('WARNING in scatter_one2ones: all(XVAL/YVAL) are non-finite; skipping this dataset')
    return 1

  colors = [
    [0.0000, 0.4470, 0.7410],
    [0.8500, 0.3250, 0.0980],
    [0.9290, 0.6940, 0.1250],
    [0.4940, 0.1840, 0.5560],
    [0.4660, 0.6740, 0.1880],
    [0.3010, 0.7450, 0.9330],
    [0.6350, 0.0780, 0.1840],
  ]
  linespecs = ['-','--','-.',':']
  markers = ['*','+','o','.']
  msizes  = [0.5, 0.5, 0.5, 3 ]

  for i, YVAL in enumerate(YVALS):
    col = colors[np.mod(i, len(colors))]
    mind = np.mod(i, len(markers))
    ax.plot(XVAL, YVAL, markers[mind], color = col,
            markersize = msizes[mind], alpha=0.5)

  if XLAB != '':
    label = XLAB
    if UNITS != '': label = label+' ('+UNITS+')'
    ax.set_xlabel(label, fontsize=6)
  if YLAB != '':
    label = YLAB
    if UNITS != '': label = label+' ('+UNITS+')'
    ax.set_ylabel(label, fontsize=6)
  if len(LEG) > 0:
    ax.legend(LEG, loc='upper left', fontsize=5)

  ymin, ymax = ax.get_ylim()
  xmin, xmax = ax.get_xlim()
  xymin=min(xmin, ymin)
  xymax=max(xmax, ymax)

  #Could adjust limits/ticks for better aesthetics
  #round_nmbr = 5
  #xymin=np.floor(min(xmin, ymin) / round_nmbr) * round_nmbr
  #xymax=np.ceil(max(xmax, ymax) / round_nmbr) * round_nmbr

  # Add linear regression and statistics
  predictor = np.arange(xymin, xymax, (xymax - xymin) / 15.0)
  tx = 0.98
  ty = 0.02
  nline = ''
  for j, YVAL in enumerate(reversed(YVALS)):
    i = len(YVALS) - j - 1
    if NVALS[j]==0: continue
    del not_nan
    not_nan = np.isfinite(XVAL) & np.isfinite(YVAL)

    # Add linear fit to YVAL vs. XVAL
    p = np.polyfit(XVAL[not_nan], YVAL[not_nan], 1)
    predicted = p[0] * predictor + p[1]
    col0 = colors[np.mod(i, len(colors))]
    linespec = linespecs[np.mod(j, len(linespecs))]
#    bright = 0.5
#    col = bright * np.asarray([1., 1., 1.]) + (1. - bright) * np.asarray(col0)
    dimmer = 0.35
    col = (1. - dimmer) * np.asarray(col0)
    ax.plot(predictor, predicted, linespec, color = col, lw=1.2)

    # Add statistics for YVAL vs. XVAL
    stat = 'N = %d\nslope: %0.2f  '%(NVALS[j], p[0])
    if show_stats:
      RMSE = np.sqrt( np.sum( np.square(YVAL[not_nan] - XVAL[not_nan]) ) / NVALS[j] )
      BIAS = np.sum( YVAL[not_nan] - XVAL[not_nan] ) / NVALS[j]
      stat = stat+'\nRMSE: %0.2f \nBIAS: %0.2f'%(RMSE, BIAS)
    ax.text(tx, ty, stat+nline,
      {'color': col, 'fontsize': fsize_lab},
      ha='right', va='bottom', backgroundcolor=[1, 1, 1, 0.2],
      clip_on=True, transform=ax.transAxes)
    nline = nline + ''.join('\n' * stat.count('\n')) + '\n'

  ax.plot([xymin, xymax],[xymin, xymax],'k:', lw=0.5)
  ticks = ax.get_yticks()
  ticks = ticks[np.logical_and(ticks >= xymin, ticks<=xymax)]
  ax.set_xticks(ticks)
  ax.set_yticks(ticks)
  for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(5)
    tick.label.set_rotation('vertical')
  for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(5)

  ax.set_xlim(xymin, xymax)
  ax.set_ylim(xymin, xymax)
  ax.grid()

  return 0

def main():
  readdata()

if __name__ == '__main__': main()
