import os
import sys
sys.path.insert(1, '../')
import numpy as np
from copy import deepcopy
import fnmatch
import config as conf
import basic_plot_functions
import var_utils as vu
import h5py as h5

'''
Directory structure and file names for ctest:
test/
 ├── Data/
 │   ├──ufo/
 │       ├─testinput_tier_1/
 │           ├── sondes_obs_2018041500_m.nc4 (or sondes_obs_2018041500.h5)
 │           ├── satwind_obs_2018041500_m.nc4 (or satwind_obs_2018041500.h5)
 │           ├── ...
 ├── graphics/
 │   ├── basic_plot_functions.py
 │   │──...
 │   ├── standalone
 │   │   ├── plot_obs_loc.py
 │   │   ├── ...

Directory structure and file names for cycling:
test/
 ├── dbIn/
 │   ├── sondes_obs_2018041500.nc4 (or sondes_obs_2018041500.h5)
 │   ├── satwind_obs_2018041500.nc4 (or satwind_obs_2018041500.h5)
 │   ├── ...
 ├── graphics/
 │   ├── basic_plot_functions.py
 │   │──...
 │   ├── standalone
 │   │   ├── plot_obs_loc.py
 │   │   ├── ...


How to run it:
python plot_obs_loc.py cycling 2018041500
python plot_obs_loc.py ctest   2018041500
'''
test = str(sys.argv[1])
Date = str(sys.argv[2])

def readdata():
    '''
    Observation file name used in ctest:  |  Observation file name used in cycling:
        aircraft_obs_2018041500_m.nc4     |      aircraft_obs_2018041500.h5
        amsua_n19_obs_2018041500_m.nc4    |      amsua_n19_obs_2018041500.h5
        satwind_obs_2018041500_m.nc4      |      satwind_obs_2018041500.h5
        sondes_obs_2018041500_m.nc4       |      sondes_obs_2018041500.h5
        gnssro_obs_2018041500_s.nc4       |      gnssro_obs_2018041500.h5
    '''
    obsfiles = []
    if (test == 'ctest'):
        suffix = '_m'
        string = '_obs_'+Date+suffix+'.nc4'
        nChar1 = -21  # remove "_obs_2018041500_m.nc4" to get obstype string
        nChar2 = -4   # remove '.nc4' to get output string
        dataPath = '../../Data/ufo/testinput_tier_1/'
    elif (test == 'cycling'):
        suffix = ''
        string = '_obs_'+Date+suffix+'.h5'
        nChar1 = -18  # remove "_obs_2018041500.h5" to get obstype string
        nChar2 = -3   # remove '.h5' to get output string
        dataPath = '../../dbIn/'
    PreQCMaxvalueConv = 3
    PreQCMinvalueConv = -90
    PreQCMaxvalueAmsua = 0

    #search all obs files and get obs types from  Data dir:
    allobstypes = []
    for files in os.listdir(dataPath):
        if fnmatch.fnmatch(files, '*'+string):
            allobstypes.append(files[:nChar1])

    #get obs types with 'process': True from dictionary conf.DiagSpaceConfig
    ObsSpaceDict = {}
    obsfiles_prefix = []
    for (key,baseval) in conf.DiagSpaceConfig.items():
        if baseval['process'] and baseval['DiagSpaceGrp'] != conf.model_s:
            ObsSpaceDict = deepcopy(baseval)
            obsfiles_prefix.append(key)

    #get the file name we want to plot:
    match = set(allobstypes) & set(obsfiles_prefix)
    print('match=', match)
    obsfiles = [x + string for x in match]

    if (test == 'ctest'):
        obsfiles = [f.replace('gnssro_obs_2018041500_m.nc4', 'gnssro_obs_2018041500_s.nc4') for f in obsfiles]
    print(obsfiles)

    for file_name in obsfiles:
        nc = h5.File(dataPath+file_name, 'r')
        print('Plotting:', file_name)
        varlist = []
        for node in nc:
          if type(nc[node]) is h5._hl.group.Group:
            for var in nc[node]:
              varlist += [node+'/'+var]
        stationidnc = []
        if 'MetaData/stationIdentification' in varlist:
            stationidnc = nc['MetaData/stationIdentification']
        #for gnss:
        elif 'MetaData/occulting_sat_is' in varlist:
            stationidnc = nc['MetaData/occulting_sat_is']
            if (test == 'cycling'):
                recordNum = nc['MetaData/sequenceNumber']
        elif 'MetaData/satelliteIdentifier' in varlist:
            stationidnc = nc['MetaData/satelliteIdentifier']
        #for radiances:
        else:
            nstation = 0

        obstype = file_name[:nChar1]
        latnc = nc['MetaData/latitude']
        lonnc = nc['MetaData/longitude']

        ObsSpaceInfo = conf.DiagSpaceConfig.get(obstype,conf.nullDiagSpaceInfo)
        channels = ObsSpaceInfo.get('channels',[vu.miss_i])
        #select variables with the suffix 'ObsValue'
        obslist = [obs for obs in varlist if (obs[:8] == 'ObsValue')]
        #print('obslist=',obslist)

        obs_type = file_name[:nChar1]
        out_name = file_name[:nChar2]

        for var in obslist:
            var_name = var[9:]  # e.g. remove 'ObsValue/' from 'ObsValue/air_temperature'
            PreQC = 'PreQC/'+var_name

            if var_name == 'refractivity':
                var_unit = 'N-unit'
            elif var_name == 'bending_angle':
                var_unit = 'Radians'
            else:
                var_unit = vu.varDictObs[var_name][0]

            levbin='all'

            if channels[0] == vu.miss_i:
                obsnc = nc[var]
                obsnc = np.asarray(obsnc)
                stationidnc_array = []
                recordnc_array = []
                if (obstype == 'gnssro' or obstype == 'gnssroref'):
                    obsnc[np.less(obsnc, -999)] = np.NaN
                    stationidnc_array=np.asarray(stationidnc).astype(str)
                    if (test == 'cycling'):
                        recordnc_array=np.asarray(recordNum).astype(str)
                        recordnc_array[np.isnan(obsnc)]= np.NaN
                        nrecord = len(set(recordnc_array)) -1
                else:
                    PreQCnc = nc[PreQC]
                    obsnc[np.greater(PreQCnc, PreQCMaxvalueConv)] = np.NaN
                    obsnc[np.less(PreQCnc,PreQCMinvalueConv)] = np.NaN
                    stationidnc_array=np.asarray(stationidnc)
                stationidnc_array[np.isnan(obsnc)]= np.NaN
                nstation = len(set(stationidnc_array)) -1 # -1: 'nan' is also included, so remove it
                if (obstype == 'satwind' or obstype == 'satwnd'):
                    nstation = 0
                if ((obstype == 'gnssro' or obstype == 'gnssroref') and test == 'cycling'):
                    basic_plot_functions.plotDistri(latnc,lonnc,obsnc,obs_type,var_name,var_unit,out_name,nrecord,levbin)
                else:
                    basic_plot_functions.plotDistri(latnc,lonnc,obsnc,obs_type,var_name,var_unit,out_name,nstation,levbin)
            else:
                for channel in channels:
                    obsnc = nc[var][:,channel-1]
                    PreQCnc = nc[PreQC][:,channel-1]
                    obsnc = np.asarray(obsnc)
                    obsnc[np.greater(PreQCnc, PreQCMaxvalueAmsua)] = np.NaN
                    var_name = var_name +'_ch'+ str(channel)
                    nstation = 0
                    basic_plot_functions.plotDistri(latnc,lonnc,obsnc,obs_type,var_name,var_unit,out_name,nstation,levbin)
                    var_name = var[9:]
def main():
    readdata()

if __name__ == '__main__': main()
