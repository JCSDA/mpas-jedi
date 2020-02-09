import numpy as np
import matplotlib
matplotlib.use('AGG')
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import plot_utils as pu

#This script includes basic plotting functions.

def plotDistri(lats,lons,values,ObsType,VarName,var_unit,out_name,nstation,levbin):
#================================================================
#INPUTS:
# lats     - latitude
# lons     - longitude
# values   - values will be plotted
# ObsType - observation type
# VarName - variable name
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
            %(ObsType,VarName,var_unit,len(values)),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)
    else:
        plt.text(0.5, 1.25, '%s   variable: %s %s nlocs:%s nstation:%s' \
            %(ObsType,VarName,var_unit,len(values),nstation),    \
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

    plt.savefig('distri_%s_%s_%s.png'%(VarName,out_name,levbin),dpi=200,bbox_inches='tight')
    plt.close()

def plotTimeserial2D(Stats,xlabeltime,ylevels,VarName):
#================================================================
#INPUTS:
# Stats      - statistics
# xlabeltime - time labels for x-axis 
# ylevels    - vertical levels for y-axis
# VarName    - variable name
#================================================================
    zgrid = np.loadtxt("/glade/work/jban/pandac/fix_input/graphics/zgrid_v55.txt")

    fig, ax1 = plt.subplots()

    xarray = range(len(xlabeltime))
    valuemin = np.amin(Stats)
    valuemax = np.amax(Stats)

    if (valuemin > 0 or valuemax < 0):
        color = 'rainbow'
        plt.contourf(xarray,ylevels,Stats,40,vmin=valuemin, vmax=valuemax,cmap=color)
    else:
        cmap = 'coolwarm'
        norm = matplotlib.colors.DivergingNorm(vmin=valuemin, vcenter=0, vmax=valuemax)
        plt.contourf(xarray,ylevels,Stats,40,vmin=valuemin, vmax=valuemax,norm=norm,cmap=cmap)

    xarray = range(len(xlabeltime))
    major_ticks = np.arange(0, 56, 5)
    ax1.set_yticks(major_ticks)
    ax1.set_ylim([0,54])
    ax1.set_ylabel('Vertical level',fontsize=15)

    ax2 = ax1.twinx()
    ax2.set_yticks(major_ticks-1)
    ax2.set_yticklabels((zgrid[::5]).astype(int))

    ax2.set_ylabel('Height (m)',fontsize=13)

    FCDay = ''.join(VarName.split("_")[1:][:-3])
    if (FCDay == 'day0.0'):
        ax1.set_xlabel('Analysis Time',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    elif (FCDay == 'day0.25'):
        ax1.set_xlabel( '6h Forecast',fontsize=15)
        ax1.set_xticks(xarray[::4])
        ax1.set_xticklabels(xlabeltime[::4],rotation=90)
    else:
        ax1.set_xlabel( 'Lead Time',fontsize=15)

    plt.colorbar(extend='both',orientation="horizontal",pad=0.2)
    ax1.grid(True)
    region = ''.join(VarName.split("_")[2:][:-2])
    var    = ''.join(VarName.split("_")[3:][:-1])
    stats  = ''.join(VarName.split("_")[4:])
    plt.title(stats+'  variable:'+pu.varDictModel[var][1]+'('+ pu.varDictModel[var][0]+')  '+region, fontsize = 12)
    plt.savefig(VarName+'_TS_2d.png',dpi=200,bbox_inches='tight')
    plt.close()
