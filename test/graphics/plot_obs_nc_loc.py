import os, json
import sys
import numpy
import numpy as np
from copy import deepcopy
from netCDF4 import Dataset
import matplotlib
matplotlib.use('AGG')
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch
import basic_plot_functions
import plot_utils as pu

'''
Directory Structure:
test/
 ├── Data/
 │   ├── sondes_obs_2018041500_m.nc4
 │   ├── satwind_obs_2018041500_m.nc4
 │   ├── ...
 ├── graphics/
 │   ├── plot_obs_nc_loc.py
 │   ├── ...
How to run it:
python plot_obs_nc_loc.py cycling 2018041500
python plot_obs_nc_loc.py ctest   2018041500
'''
test = str(sys.argv[1])
Date = str(sys.argv[2])

def readdata():
    '''
    Observation file name used in ctest:  |  Observation file name used in cycling:
        aircraft_obs_2018041500_m.nc4     |      aircraft_obs_2018041500.nc4
        amsua_n19_obs_2018041500_m.nc4    |      amsua_n19_obs_2018041500.nc4
        satwind_obs_2018041500_m.nc4      |      satwind_obs_2018041500.nc4
        sondes_obs_2018041500_m.nc4       |      sondes_obs_2018041500.nc4
        gnssro_obs_2018041500_s.nc4       |      gnssro_obs_2018041500.nc4
    '''
    obsfiles = []
    if (test == 'ctest'):
        suffix = '_m'
    elif (test == 'cycling'):
        suffix = ''
    string = '_obs_'+Date+suffix+'.nc4'

    #search all obs files and get obs types from  Data dir:
    allobstypes = []
    for files in os.listdir('../Data/'):
        if fnmatch.fnmatch(files, '*'+string):
            if (test == 'ctest'):
                allobstypes.append(files[:-21])
            elif (test == 'cycling'):
                allobstypes.append(files[:-19])

    #get obs types with 'process': True from dictionary DiagSpaceDict
    ObsSpaceDict = {}
    obsfiles_prefix = []
    for (key,baseval) in pu.DiagSpaceDict.items():
        if baseval['process'] and baseval['DiagSpaceGrp'] != pu.model_s:
            ObsSpaceDict = deepcopy(baseval)
            obsfiles_prefix.append(key)

    #get the file name we want to plot:
    match = set(allobstypes) & set(obsfiles_prefix)
    print('match=', match)
    obsfiles = [x + string for x in match]

    if (test == 'ctest'):
        obsfiles = [obsfiles.replace('gnssro_obs_2018041500_m.nc4', 'gnssro_obs_2018041500_s.nc4') for obsfiles in obsfiles]
    print(obsfiles)

    for index, file_name in enumerate(obsfiles):
        nc = Dataset("../Data/"+file_name, 'r')
        print('Plotting:', file_name)

        varlist = nc.variables.keys()
        if 'station_id@MetaData' in varlist:
            stationidnc =[ b''.join(i[0:8]) for i in nc.variables['station_id@MetaData'][:]]
        else:
            nstation = 0

        obstype = file_name[:-19]
        #print(obstype)
        latnc = nc.variables['latitude@MetaData']
        lonnc = nc.variables['longitude@MetaData']

        lonnc = numpy.asarray(lonnc)
        for i in range(len(lonnc)):
            if lonnc[i] > 180:
                lonnc[i] = lonnc[i]-360

        ObsSpaceInfo = pu.DiagSpaceDict.get(obstype,pu.nullDiagSpaceInfo)
        channels = ObsSpaceInfo.get('channels',[pu.miss_i])
        #select variables with the suffix 'ObsValue'
        if len(channels) == 0:
            continue
        if channels[0] == pu.miss_i:
            obslist = [obs for obs in varlist if (obs[-8:] == 'ObsValue')]
        else:
            obslist = ['brightness_temperature_{0}@ObsValue'.format(i) for i in channels]
        for var in obslist:
            print(var)
            obsnc = nc.variables[var]
            obsnc = numpy.asarray(obsnc)
            obsnc[np.greater(obsnc,1.0e+8)] = np.NaN

            if 'station_id@MetaData' in varlist:
                stationidnc_array=np.asarray(stationidnc)
                stationidnc_array[np.isnan(obsnc)]= np.NaN
                nstation = len(set(stationidnc_array)) -1   # -1: 'nan' is also included, so remove it

            if (test == 'ctest'):
                obs_type = file_name[:-21]
            elif (test == 'cycling'):
                obs_type = file_name[:-19]
            var_name = var[:-9]
            out_name = file_name[:-4]
            if '_'.join(var_name.split("_")[:-1]) == 'brightness_temperature':
                var_unit = pu.varDict['brightness_temperature'][0]
            else:
                var_unit = pu.varDict[var_name][0]

            levbin='all'
            basic_plot_functions.plotDistri(latnc,lonnc,obsnc,obs_type,var_name,var_unit,out_name,nstation,levbin)
def main():
    readdata()

if __name__ == '__main__': main()
