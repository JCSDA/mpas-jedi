import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes

#This script includes basic plotting functions.

def plotDistri(lats,lons,values,OBS_TYPE,VAR_NAME,var_unit,out_name,nstation,levbin):
#================================================================
#INPUTS:
# lats     - latitude
# lons     - longitude
# values   - values will be plotted
# OBS_TYPE - observation type
# VAR_NAME - variable name
# var_unit - variable units
# out_name - will be included in output file name. It can be experiment name.
# nstation - station numbers for sondes.
# levbin   - plot all levels together (levbin=all); or plot every level.
#================================================================

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
    sc=m.scatter(lons[:],lats[:],c=values[:],s=6,cmap=cm,
        zorder=10) #10 ; zorder determines the order of the layer. If not set,
                        #the point on continent will be blocked

#create axes for colorbar======================================================
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top",
        size="3%", #width of the colorbar
        pad=0.3)   #space between colorbar and graph
    plt.colorbar(sc,cax=cax,ax=ax,orientation='horizontal')
    cax.xaxis.set_ticks_position('top')  #set the orientation of the ticks of the colorbar

    plt.savefig('distri_%s_%s_%s.png'%(VAR_NAME,out_name,levbin),dpi=200,bbox_inches='tight')
    plt.close()

#TODO:def plot2D():
