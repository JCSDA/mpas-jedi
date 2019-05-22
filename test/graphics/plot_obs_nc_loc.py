import os, json
import numpy
import numpy as np
from netCDF4 import Dataset
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import fnmatch

# columns: obs_type,       plot or not,    select channels for radiance
plotdict = { \
    'sondes':                 ['T',        '0'   ] \
  , 'aircraft':               ['T',        '0'   ] \
  , 'satwind':                ['T',        '0'   ] \
  , 'gnssro':                 ['T',        '0'   ] \
  , 'airs_aqua':              ['T',        [1,6,7]]\
  , 'amsua_n19':              ['T',        [4,5,6,7,9,10,11,12,13,14]] \
  , 'cris-fsr_npp':           ['T',        [24,26,28,32,37,39]]        \
  , 'hirs4_metop-a':          ['T',        range(1,16)] \
  , 'iasi_metop-a':           ['T',        [16,29,32,35,38,41,44]]     \
  , 'mhs_n19':                ['T',        range(1,6)]  \
  , 'seviri_m08':             ['T',        [5]]         \
  , 'sndrd1_g15':             ['T',        range(1,16)] \
  , 'sndrd2_g15':             ['F',        range(1,16)] \
  , 'sndrd3_g15':             ['F',        range(1,16)] \
  , 'sndrd4_g15':             ['F',        range(1,16)] \
    }

# columns: var_name            unit_used   abbr.
vardict = { \
    'air_temperature':        [ '(K)',     'T'   ] \
  , 'virtual_temperature':    [ '(K)',     'T'   ] \
  , 'eastward_wind':          [ '(m/s)',   'U'   ] \
  , 'northward_wind':         [ '(m/s)',   'V'   ] \
  , 'specific_humidity':      [ '(kg/kg)', 'Q'   ] \
  , 'refractivity':           [ '(N-unit)','Ref' ] \
  , 'bending_angle':          [ '(rad)',   'Bnd' ] \
  , 'brightness_temperature': [ '(K)',     'BT'  ] \
  , 'aerosol_optical_depth_4':[ '   ',     'AOD' ] \
    }

def readdata():
    '''
    Observation file name:
             aircraft_obs_2018041500_m.nc4
             amsua_n19_obs_2018041500_m.nc4
             aod_obs_2018041500_m.nc4
             satwind_obs_2018041500_m.nc4
             sondes_obs_2018041500_m.nc4
             gnssro_obs_2018041500_s.nc4
    '''
    obsfiles = []
    string = '_obs_2018041500_m.nc4'

#   newdict: filter out 'F' lists.
    newplotdict = dict(filter(lambda x: x[1][0] == 'T', plotdict.items()))
    obsfiles = [x + string for x in newplotdict.keys()]
    obsfiles = [obsfiles.replace('gnssro_obs_2018041500_m.nc4', 'gnssro_obs_2018041500_s.nc4') for obsfiles in obsfiles]
    print(obsfiles)

    for index, file_name in enumerate(obsfiles):
        nc = Dataset("../Data/"+file_name, 'r')
        print 'Plotting:', file_name

        varlist = nc.variables.keys()
        if 'station_id@MetaData' in varlist:
            stationidnc =[ ''.join(i)  for i in nc.variables['station_id@MetaData'] ]
            nstation = len(set(stationidnc))
        else:
            nstation = 0

        obstype = str(file_name.split("_")[:1])
        #print(obstype)
        latnc = nc.variables['latitude@MetaData']
        lonnc = nc.variables['longitude@MetaData']

        lonnc = numpy.asarray(lonnc)
        for i in range(len(lonnc)):
            if lonnc[i] > 180:
                lonnc[i] = lonnc[i]-360

        #select variables with the suffix 'ObsValue'
        obslist = [obs for obs in varlist if (obs[-8:] == 'ObsValue')]

        if type(plotdict[newplotdict.keys()[index]][1]) is not str:
            obslist=(['brightness_temperature_{0}@ObsValue'.format(i) for i in plotdict[newplotdict.keys()[index]][1]])
        #print 'check obslist=', obslist
        for var in obslist:
            print(var)
            obsnc = nc.variables[var]
            obsnc = numpy.asarray(obsnc)
            obsnc[np.greater(obsnc,1.0e+8)] = np.NaN

            obs_type = file_name[:-21]
            var_name = var[:-9]
            out_name = file_name[:-4]
            if '_'.join(var_name.split("_")[:-1]) == 'brightness_temperature':
                var_unit = vardict['brightness_temperature'][0]
            else:
                var_unit = vardict[var_name][0]

            plot(latnc,lonnc,obsnc,obs_type,var_name,var_unit,out_name,nstation)

def plot(lats,lons,values,OBS_TYPE,VAR_NAME,var_unit,out_name,nstation):
#set map=======================================================================
    fig,ax=plt.subplots(figsize=(8,8))
    m=Basemap(projection='cyl',     #map projection
        llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90, #the lat lon of leftlow and upper right corner
#       llcrnrlon=-20,llcrnrlat=10,urcrnrlon=60,urcrnrlat=80,   #zoom in
        resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
    m.drawcoastlines()  #draw coastline
#   m.drawmapboundary(fill_color='lightblue') #the whole map is filled with the specified color
#   m.fillcontinents(color='wheat',lake_color='lightblue') #the continents are filled

#set the lat lon grid and the ticks of x y axies ==============================
    parallels=np.arange(-90,90,30)      #lat
    meridians=np.arange(-180,180,60)    #lon
    m.drawparallels(parallels,
        labels=[True,True,True,True])   #Turn the ticks on or off
    m.drawmeridians(meridians,
        labels=[True,True,True,True])

#set title  ===================================================================
    if nstation == 0:
        plt.text(0.5, 1.25, '%s   variable: %s %s nlocs:%s' \
            %(OBS_TYPE,VAR_NAME,var_unit,len(values)),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)
    else:
        plt.text(0.5, 1.25, '%s   variable: %s %s nlocs:%s nstation:%s' \
            %(OBS_TYPE,VAR_NAME,var_unit,len(values),nstation),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)

#draw points onto map =========================================================
    cm=plt.cm.get_cmap('rainbow')
    sc=m.scatter(lons,lats,c=values,s=6,cmap=cm,
        zorder=10) #10 ; zorder determines the order of the layer. If not set,
                        #the point on continent will be blocked

#create axes for colorbar======================================================
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top",
        size="3%", #width of the colorbar
        pad=0.3)   #space between colorbar and graph
    plt.colorbar(sc,cax=cax,ax=ax,orientation='horizontal')
    cax.xaxis.set_ticks_position('top')  #set the orientation of the ticks of the colorbar

    plt.savefig('distri_%s_%s.png'%(VAR_NAME,out_name),dpi=200,bbox_inches='tight')
    plt.close()

def main():
    readdata()

if __name__ == '__main__': main()
