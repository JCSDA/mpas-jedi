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

def readdata():

    # Observation file name format:
    #                *_obs_*_m.nc4
    #         aircraft_obs_2018041500_m.nc4
    #         amsua_obs_n19_2018041500_m.nc4
    #         aod_obs_2018041500_m.nc4
    #         satwind_obs_2018041500_m.nc4
    #         sondes_obs_2018041500_m.nc4
    obsfiles = []
    for files in os.listdir('../Data/'):
        if fnmatch.fnmatch(files, '*_obs_*_m.nc4'):
            obsfiles.append('../Data/'+files)
    for file_name in obsfiles:
        nc = Dataset(file_name, 'r')
        print 'file_name', file_name
        latnc = nc.variables['latitude@MetaData']
        lonnc = nc.variables['longitude@MetaData']

        lonnc = numpy.asarray(lonnc)
        for i in range(len(lonnc)):
            if lonnc[i] > 180:
                lonnc[i] = lonnc[i]-360

        varlist = nc.variables.keys()
        #select variables with the suffix 'ObsValue'
        obslist = [obs for obs in varlist if (obs[-8:] == 'ObsValue')]
        #print(obsvalue)
        for var in obslist:
            print(var)
            obsnc = nc.variables[var]
            obsnc = numpy.asarray(obsnc)
            obsnc [obsnc== 9.96920997e+36] = numpy.NaN
 
            if (file_name[8:][:5] == 'amsua'):
                plot(latnc,lonnc,obsnc,file_name[8:][:-25],var[:-10])
            else:
                plot(latnc,lonnc,obsnc,file_name[8:][:-21],var[:-9])

def plot(lats,lons,values,OBS_TYPE,VAR_NAME):
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
    plt.text(0.5, 1.25, '%s   variable: %s'%(OBS_TYPE,VAR_NAME) ,
        horizontalalignment='center',
        fontsize=12,
        transform = ax.transAxes)

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

    plt.savefig('distri_%s_%s.png'%(OBS_TYPE,VAR_NAME),dpi=200,bbox_inches='tight')
    plt.close()

def main():
    readdata()

if __name__ == '__main__': main()
