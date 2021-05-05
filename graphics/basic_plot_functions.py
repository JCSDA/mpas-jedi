#!/usr/bin/env python3

import datetime as dt
import logging
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib
matplotlib.use('AGG')
import matplotlib.axes as maxes
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import plot_utils as pu
import var_utils as vu
import os

_logger = logging.getLogger(__name__)

cmGray     = plt.cm.get_cmap("gist_gray")
cmRainbow  = plt.cm.get_cmap("gist_rainbow")
cmSpectral = plt.cm.get_cmap("nipy_spectral")
cmHeat     = plt.cm.get_cmap("gist_heat")
cmOcean    = plt.cm.get_cmap("ocean")
cmNCAR     = plt.cm.get_cmap("gist_ncar")

WhiteBlack1 = cmGray(np.linspace(1.0,0.0,17)) # white to black (-90 to -74 C)
BlackRed    = cmHeat(np.linspace(0.0,0.5,10)) #black to red (-74 to -65 C)
ROYG        = cmSpectral(np.linspace(0.9,0.43,27)) # red, orange, yellow, green, blue (-65 to -39 C)
#GreenBlue   = cmNCAR(np.linspace(0.05,0.1,8)) # green to blue (-39 to -32 C)
#BlueCyan    = cmRainbow(np.linspace(0.8,0.6,13)) # blue to cyan (-32 to -20 C)
GreenBlueCyan = cmNCAR(np.linspace(0.05,0.2,20)) # green to blue (-39 to -20 C)
#WhiteBlack2 = cmGray(np.linspace(0.9,0.0,51)) # white to black (-20 to 30 C)
MVW = cmNCAR(np.linspace(0.8,0.98,21)) # magenta to violet to white (-20 to 0 C)
WhiteBlack2 = cmGray(np.linspace(0.9,0.0,31)) # white to black (0 to 30 C)

#btcolors = np.concatenate((WhiteBlack1, BlackRed, ROYG, GreenBlue, BlueCyan, WhiteBlack2))
#btcolors = np.concatenate((WhiteBlack1, BlackRed, ROYG, GreenBlueCyan, WhiteBlack2))
btcolors = np.concatenate((WhiteBlack1, BlackRed, ROYG, GreenBlueCyan, MVW, WhiteBlack2))

btCMap = colors.ListedColormap(btcolors)

#This script includes basic plotting functions.

distriZooms = {}

#Full Earth
distriZooms['default'] = {
    'cLon': None,
    'minLon': -180,
    'maxLon': 180,
    'minLat': -90,
    'maxLat': 90,
}
distriZooms['abi'] = {
    'cLon': -75.2,
    'minLon': None,
    'maxLon': None,
    'minLat': None,
    'maxLat': None,
}
distriZooms['ahi'] = {
    'cLon': 140.7,
    'minLon': None,
    'maxLon': None,
    'minLat': None,
    'maxLat': None,
}

def plotDistri(lats,lons,values, \
               ObsType,VarName,var_unit,out_name,nstation,levbin, \
               dmin=None,dmax=None,dotsize=6,color="rainbow"):
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
# dmin, dmax  - min/max values of colorbars, optional
# dotsize  - dot size, optional
# color    - color scheme, optional
#================================================================
# For some plots that need to change longitude from [-180,180] to [0,360]
#    tmp = np.logical_not(lons > 0)
#    lons[tmp] = lons[tmp] + 360

#set map=======================================================================
    cLon = distriZooms['default']['cLon']
    minLon = distriZooms['default']['minLon']
    maxLon = distriZooms['default']['maxLon']
    minLat = distriZooms['default']['minLat']
    maxLat = distriZooms['default']['maxLat']

    for key, val in distriZooms.items():
        if key in ObsType:
            cLon = val['cLon']
            minLon = val['minLon']
            maxLon = val['maxLon']
            minLat = val['minLat']
            maxLat = val['maxLat']

    parallels=np.arange(-90,90,30)      #lat

    if cLon is not None:
        fig,ax=plt.subplots(figsize=(5,5))
        m=Basemap(projection='nsper',     #map projection
            lon_0 = cLon,
            lat_0 = 0.0,
            resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
        #set the lat lon grid and the ticks of x y axies ==============================
        m.drawparallels(parallels, linewidth=0.2, dashes=[1,3])
        meridians=np.arange(-180,180,30)    #lon
        m.drawmeridians(meridians, linewidth=0.2, dashes=[1,3])
    else:
        fig,ax=plt.subplots(figsize=(8,8))
        m=Basemap(projection='cyl',     #map projection
            llcrnrlon=minLon, llcrnrlat=minLat, #the lat lon of leftlower corner
            urcrnrlon=maxLon, urcrnrlat=maxLat, #the lat lon of rightupper corner
            resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
        #set the lat lon grid and the ticks of x y axies ==============================
        m.drawparallels(parallels, linewidth=0.2, dashes=[1,3],
            labels=[True,True,True,True])   #Turn the ticks on or off
        meridians=np.arange(-180,180,60)    #lon
        m.drawmeridians(meridians, linewidth=0.2, dashes=[1,3],
            labels=[True,True,True,True])

    m.drawcoastlines(linewidth=0.2,zorder=10)  #draw coastline
#   m.drawmapboundary(fill_color='lightblue') #the whole map is filled with the specified color
#   m.fillcontinents(color='wheat',lake_color='lightblue') #the continents are filled


#set title  ===================================================================
    if nstation == 0:
        plt.text(0.5, 1.25, '%s   %s %s nlocs:%s' \
            %(ObsType,VarName,var_unit,len(values[~np.isnan(values)])),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)
    else:
        plt.text(0.5, 1.25, '%s   %s %s nlocs:%s nstation:%s' \
            %(ObsType,VarName,var_unit,len(values[~np.isnan(values)]),nstation),    \
            horizontalalignment='center', \
            fontsize=12, transform = ax.transAxes)

#draw points onto map =========================================================
    if color == "BT":
        if ("abi" in ObsType or "ahi" in ObsType):
            cm = btCMap
            if dmin is None: dmin = 183
            if dmax is None: dmax = 303
        else:
            cm = plt.cm.get_cmap("gist_ncar")
            if dmin is None: dmin = 190
            if dmax is None: dmax = 270
    else:
        cm = plt.cm.get_cmap(color)

    finite = np.isfinite(values)
    if ((("abi" in ObsType or "ahi" in ObsType)
         and finite.sum() > 4e4)
        or "model" in ObsType):
        # option 1: smoothed contours (note: color bar is not quite right)
        # sc=m.contourf(lons[finite], lats[finite], values[finite],
        #               cm.N, cmap = cm, vmin = dmin, vmax = dmax,
        #               latlon = True, tri = True, extend='both')

        # option 2: pixel contours
        # first sort by longitude to avoid bug for cyclic projections in basemap
        lonsPlot = lons[finite]
        lonsPlot[lonsPlot > 180.0] -= 360.0 # fixes latitude swap bug for cyclic projections
        latsPlot = lats[finite]
        valuesPlot = values[finite]
        lonSort = np.argsort(lonsPlot)

        sc = m.pcolor(lonsPlot[lonSort], latsPlot[lonSort], valuesPlot[lonSort],
                      cmap = cm, vmin = dmin, vmax = dmax,
                      latlon = True, tri = True)

    else:
        sc = m.scatter(lons[finite], lats[finite], c=values[finite],
                       latlon = True,
                       s = dotsize, cmap = cm, vmin = dmin, vmax = dmax,
                       marker = '.', linewidth = 0,
                       zorder=11) #10 ; zorder determines the order of the layer. If not set,
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


def scatterMapFields(
  lons, lats, fields,
  filename,
  minLon = -180., maxLon = 180.,
  minLat = -90., maxLat = 90.,
  cLon = None,
  projection = 'default',
  dmin = None, dmax = None,
  markers = {},
  sizes = {},
  cmap = 'gist_ncar',
  cbarType = None,
  c = {},
  logVLim = 1.e-12,
  ):

  # set up projection
  if cLon is not None:
    fig,ax=plt.subplots(figsize=(5,5))
    if projection == 'default':
      proj = 'nsper'
    else:
      proj = projection
    m=Basemap(projection=proj,
              lon_0 = cLon,
              lat_0 = 0.0,
              resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
    #set the lat lon grid and the ticks of x y axies ==============================
    parallels=np.arange(-90, 90, 30)
    m.drawparallels(parallels, linewidth=0.2, dashes=[1,3], zorder = 3)
    meridians=np.arange(-180, 180, 30)
    m.drawmeridians(meridians, linewidth=0.2, dashes=[1,3], zorder = 3)
  else:
    fig,ax=plt.subplots(figsize=(8,8))
    if projection == 'default':
      proj = 'cyl'
    else:
      proj = projection
    m=Basemap(projection=proj,
              llcrnrlon=minLon, llcrnrlat=minLat, #the lat lon of leftlower corner
              urcrnrlon=maxLon, urcrnrlat=maxLat, #the lat lon of rightupper corner
              resolution='l') #c: crude; l:low; i:intermediate, h:high, f:full; sometimes can only >=h
    #set the lat lon grid and the ticks of x y axies ==============================
    parallels = np.arange(minLat, maxLat, (maxLat - minLat) / 6.)
    m.drawparallels(parallels, linewidth=0.2, dashes=[1,3],
                    labels=[True,True,True,True], zorder = 3)
    meridians = np.arange(minLon, maxLon, (maxLon - minLon) / 6.)
    m.drawmeridians(meridians, linewidth=0.2, dashes=[1,3],
                    labels=[True,True,True,True], zorder = 3)

  m.drawcoastlines(linewidth=0.2, zorder=2)  #draw coastline

    # title
#    plt.text(0.5, 1.25, '%s   %s %s nlocs:%s' \
#             %(ObsType,VarName,var_unit,len(values)),    \
#             horizontalalignment='center', \
#             fontsize=12, transform = ax.transAxes)

  assert (cbarType is None or cbarType in ['Log', 'SymLog']), \
    'scatterMapFields: invalid cbarType: '+cbarType

  for name, field in fields.items():
    finite = np.isfinite(field)
    if dmin is None:
       vmin = field[finite].min()
    else:
       vmin = dmin
    if dmax is None:
       vmax = field[finite].max()
    else:
       vmax = dmax

    if cbarType is None:
      sc = m.scatter(lons[name][finite], lats[name][finite], c=c.get(name, field[finite]),
                     latlon = True,
                     s = sizes.get(name, 1),
                     cmap = cmap, vmin = vmin, vmax = vmax,
                     marker = markers.get(name, '.'), linewidth = 0,
                     zorder=1) #10 ; zorder determines the order of the layer. If not set,
                                    #the point on continent will be blocked
    elif cbarType == 'Log':
      if vmin <= logVLim: vmin = logVLim
      field_ = field[finite]
      field_[field_ < vmin] = vmin
      sc = m.scatter(lons[name][finite], lats[name][finite], c=c.get(name, field_),
                     latlon = True,
                     s = sizes.get(name, 1),
                     cmap = cmap,
                     marker = markers.get(name, '.'), linewidth = 0,
                     norm=colors.LogNorm(vmin=vmin, vmax=vmax),
                     zorder=1) #10 ; zorder determines the order of the layer. If not set,
                                    #the point on continent will be blocked
    elif cbarType == 'SymLog':
      sc = m.scatter(lons[name][finite], lats[name][finite], c=c.get(name, field[finite]),
                     latlon = True,
                     s = sizes.get(name, 1),
                     cmap = cmap,
                     marker = markers.get(name, '.'), linewidth = 0,
                     norm=colors.SymLogNorm(vmin=vmin, vmax=vmax,
                                            linthresh=1.e-4*vmax, linscale=1.0, base=10),
                     zorder=1) #10 ; zorder determines the order of the layer. If not set,
                                    #the point on continent will be blocked

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("top",
    size="3%", #width of the colorbar
    pad=0.3)   #space between colorbar and graph
  plt.colorbar(sc, cax=cax, ax=ax, orientation='horizontal')
  cax.xaxis.set_ticks_position('top')  #set the orientation of the ticks of the colorbar

  plt.savefig(filename, dpi=200, bbox_inches='tight')
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
    plt.title(stats+'  variable:'+vu.varDictModel[var][1]+'('+ vu.varDictModel[var][0]+')  '+region, fontsize = 12)
    plt.savefig(VarName+'_TS_2d.png',dpi=200,bbox_inches='tight')
    plt.close()

maxLegendEntries = 12

###############################################################################
lenWarnSer = 0
nanWarnSer = 0
def plotSeries(fig, \
               linesVals, xVals, \
               linesLabel, \
               title="", dataLabel="y", \
               sciticks=False, logscale= False, signdef=False, \
               indepLabel="x", invert_ind_axis=False, \
               ny=1, nx=1, nplots=1, iplot=0, \
               linesValsMinCI=None, linesValsMaxCI=None, \
               dmin=np.NaN, dmax=np.NaN, \
               lineAttribOffset=0, \
               legend_inside=True,
               interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# xVals            - independent variable on x-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel        - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# logscale         - y-axis is scaled logarithmically, optional, overrides sciticks
# signdef          - whether linesVals is positive/negative definite, optional
# indepLabel       - label for xVals, optional
# invert_ind_axis  - whether to invert x-axis orientation, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = np.asarray([])
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if np.all(np.isnan(lineVals)):
            global nanWarnSer
            if nanWarnSer==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+indepLabel+"; "+linesLabel[iline])
            nanWarnSer=nanWarnSer+1
            continue
        if len(lineVals)!=len(xVals):
            global lenWarnSer
            if lenWarnSer==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+indepLabel+"; "+linesLabel[iline])
            lenWarnSer=lenWarnSer+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColor(len(linesVals),iline+lineAttribOffset)

        ax.plot(xVals, lineVals, \
                color=pColor, \
                label=linesLabel[iline], \
                ls=pu.plotLineStyle(len(linesVals),iline+lineAttribOffset), \
                linewidth=0.5)
        nLines += 1
        plotVals = np.append(plotVals, lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if np.isfinite(x) else -1.0 for x in significant])

            lineArr = np.array(lineVals)
            xArr = np.array(xVals)
            negsiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] < 0.0)],dtype=int)
            if len(negsiginds) > 0:
                ax.plot(xArr[negsiginds], lineArr[negsiginds], \
                        color=pColor, \
                        ls='', \
                        marker='v', \
                        markersize=1.5)

            possiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] > 0.0)],dtype=int)
            if len(possiginds) > 0:
                ax.plot(xArr[possiginds], lineArr[possiginds], \
                        color=pColor, \
                        ls='', \
                        marker='^', \
                        markersize=1.5)

            ax.plot(xVals, linesValsMinCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(xVals, linesValsMaxCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add horizontal zero line for unbounded quantities
    if not signdef:
        ax.plot([xVals[0], xVals[-1]], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    isLogScale = logscale
    if logscale:
        nonzero = np.logical_and(np.greater(np.abs(plotVals), 0.), np.isfinite(plotVals))
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if signdef:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0:
                    ax.set_yscale('log')
                else:
                    isLogScale = False
            else:
                ax.set_yscale('symlog')
        else:
            isLogScale = False

        if isLogScale and np.isfinite(maxdval) and maxdval > 0.:
            ax.set_ylim(None, maxdval)
            if np.abs(vmin) > 0.:
                ax.set_ylim(vmin, None)

    if not isLogScale:
        if sciticks:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if (np.isfinite(mindval) and
            np.isfinite(maxdval)):
            ax.set_ylim(mindval,maxdval)
            if maxdval-mindval < 1.0 or \
               maxdval-mindval > 100.0:
                ax.tick_params(axis='y',rotation=-35)
        ax.yaxis.get_offset_text().set_fontsize(3)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(indepLabel,fontsize=4)
    if interiorLabels or iy == ny-1:
        ax.set_ylabel(dataLabel,fontsize=4)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                           framealpha=0.4,ncol=1)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=3,frameon=False, \
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if invert_ind_axis:
        ax.invert_xaxis()

    ax.grid()

    return

###############################################################################
lenWarnProf = 0
nanWarnProf = 0
def plotProfile(fig, \
                linesVals, yVals, \
                linesLabel, \
                title="", dataLabel="x", \
                sciticks=False, logscale=False, signdef=False, \
                indepLabel="y", invert_ind_axis=False, \
                ny=1, nx=1, nplots=1, iplot=0, \
                linesValsMinCI=None, linesValsMaxCI=None, \
                dmin=np.NaN, dmax=np.NaN, \
                lineAttribOffset=0, \
                legend_inside=True,
                interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# yVals            - independent variable on y-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel        - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# logscale         - x-axis is scaled logarithmically, optional, overrides sciticks
# signdef          - whether linesVals is positive/negative definite, optional
# indepLabel       - label for yVals, optional
# invert_ind_axis  - whether to invert y-axis orientation, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = np.asarray([])
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if np.all(np.isnan(lineVals)):
            global nanWarnProf
            if nanWarnProf==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            nanWarnProf=nanWarnProf+1
            continue
        if len(lineVals)!=len(yVals):
            global lenWarnProf
            if lenWarnProf==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            lenWarnProf=lenWarnProf+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColor(len(linesVals),iline+lineAttribOffset)

        ax.plot(lineVals, yVals, \
                color=pColor, \
                label=linesLabel[iline], \
                ls=pu.plotLineStyle(len(linesVals),iline+lineAttribOffset), \
                linewidth=0.5)
        nLines += 1
        plotVals = np.append(plotVals,lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if np.isfinite(x) else -1.0 for x in significant])

            lineArr = np.array(lineVals)
            yArr = np.array(yVals)
            negsiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] < 0.0)],dtype=int)
            if len(negsiginds) > 0:
                ax.plot(lineArr[negsiginds], yArr[negsiginds], \
                        color=pColor, \
                        ls='', \
                        marker='<', \
                        markersize=1.5)

            possiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] > 0.0)],dtype=int)
            if len(possiginds) > 0:
                ax.plot(lineArr[possiginds], yArr[possiginds], \
                        color=pColor, \
                        ls='', \
                        marker='>', \
                        markersize=1.5)

            ax.plot(linesValsMinCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(linesValsMaxCI[iline], yVals, \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add vertical zero line for unbounded quantities
    if not signdef:
        ax.plot([0., 0.], [yVals[0], yVals[-1]], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    isLogScale = logscale
    if logscale:
        nonzero = np.logical_and(np.greater(np.abs(plotVals), 0.), np.isfinite(plotVals))
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if signdef:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0:
                    ax.set_xscale('log')
                else:
                    isLogScale = False
            else:
                ax.set_xscale('symlog')
        else:
            isLogScale = False

        if isLogScale and np.isfinite(maxdval) and maxdval > 0.:
            ax.set_xlim(None, maxdval)
            if np.abs(mindval) > 0.:
                ax.set_xlim(mindval, None)

    if not isLogScale:
        if sciticks:
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        if (np.isfinite(mindval) and
            np.isfinite(maxdval)):
            ax.set_xlim(mindval,maxdval)
            if maxdval-mindval < 1.0 or \
               maxdval-mindval > 100.0:
                ax.tick_params(axis='x',rotation=-35)
        ax.xaxis.get_offset_text().set_fontsize(3)


    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(dataLabel,fontsize=4)
    if interiorLabels or iy == ny-1:
        ax.set_ylabel(indepLabel,fontsize=4)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                           framealpha=0.4,ncol=1)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=3,frameon=False, \
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if invert_ind_axis:
        ax.invert_yaxis()

    ax.grid()

    return


###############################################################################
lenWarnTS=0
nanWarnTS=0
def plotTimeSeries(fig, \
                   xsDates, linesVals, \
                   linesLabel, \
                   title="", dataLabel="", \
                   sciticks=False, logscale = False, signdef=False, \
                   ny=1, nx=1, nplots=1, iplot=0, \
                   linesValsMinCI=None, linesValsMaxCI=None, \
                   dmin=np.NaN, dmax=np.NaN, \
                   lineAttribOffset=0, \
                   legend_inside=True,
                   interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# xsDates          - date x-values (list/array or list of lists/arrays
#                                   of float seconds, dt.timedelta, dt.datetime)
# linesVals        - dependent variable (list of arrays)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel        - label for linesVals, optional
# sciticks         - whether linesVals needs scientific formatting for ticks, optional
# logscale         - y-axis is scaled logarithmically, optional, overrides sciticks
# signdef          - whether linesVals is positive/negative definite, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum error bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum error bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = np.asarray([])
    nLines = 0
    jline = 0
    for iline, lineVals in enumerate(linesVals):
        if np.all(np.isnan(lineVals)):
            global nanWarnTS
            if nanWarnTS==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            nanWarnTS=nanWarnTS+1
            continue

        #float xVals
        if isinstance(xsDates[0],(list,np.ndarray)):
            xVals = pu.TDeltas2Seconds(xsDates[min([iline,len(xsDates)-1])])
        else:
            xVals = pu.TDeltas2Seconds(xsDates)

        if len(lineVals)!=len(xVals):
            global lenWarnTS
            if lenWarnTS==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            lenWarnTS=lenWarnTS+1
            continue

        if jline == 0:
            minX = xVals[0]
            maxX = xVals[-1]
        else:
            minX = min([xVals[0], minX])
            maxX = max([xVals[-1], maxX])
        jline += 1

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColor(len(linesVals),iline+lineAttribOffset)

        ax.plot(xVals, lineVals, \
                label=linesLabel[iline], \
                color=pColor, \
                ls=pu.plotLineStyle(len(linesVals),iline+lineAttribOffset), \
                linewidth=0.5)
        nLines += 1
        plotVals = np.append(plotVals, lineVals)

        # Add shaded CI regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            # test statistical significance versus zero
            if signdef:
                significant = np.empty(len(lineVals))
                significant[:] = np.NaN
            else:
                significant = np.multiply(linesValsMinCI[iline], linesValsMaxCI[iline])
            significant = np.array([x if np.isfinite(x) else -1.0 for x in significant])

            lineArr = np.array(lineVals)
            xArr = np.array(xVals)
            negsiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] < 0.0)],dtype=int)
            if len(negsiginds) > 0:
                ax.plot(xArr[negsiginds], lineArr[negsiginds], \
                        color=pColor, \
                        ls='', \
                        marker='v', \
                        markersize=1.5)

            possiginds = np.array([i for i,x in enumerate(significant)
                                   if (x > 0.0 and lineArr[i] > 0.0)],dtype=int)
            if len(possiginds) > 0:
                ax.plot(xArr[possiginds], lineArr[possiginds], \
                        color=pColor, \
                        ls='', \
                        marker='^', \
                        markersize=1.5)

            ax.plot(xVals, linesValsMinCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.plot(xVals, linesValsMaxCI[iline], \
                    color=pColor, \
                    alpha=0.4, \
                    ls='-', \
                    linewidth=0.5)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.0, alpha = 0.1)
            ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline], \
                            where=significant > 0.0, \
                            color=pColor, \
                            edgecolor=pColor, \
                            linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # standardize y-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,plotVals,signdef)

    # add horizontal zero line for unbounded quantities
    if not signdef:
        ax.plot([minX, maxX], [0., 0.], ls="--", c=".3", \
            linewidth=0.7,markersize=0)

    #axes settings
    if isinstance(xsDates[0],(list,np.ndarray)):
        pu.format_x_for_dates(ax, xsDates[0])
    else:
        pu.format_x_for_dates(ax, xsDates)

    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)
    isLogScale = logscale
    if logscale:
        nonzero = np.logical_and(np.greater(np.abs(plotVals), 0.), np.isfinite(plotVals))
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if signdef:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0:
                    ax.set_yscale('log')
                else:
                    isLogScale = False
            else:
                ax.set_yscale('symlog')
        else:
            isLogScale = False

        if isLogScale and np.isfinite(maxdval) and maxdval > 0.:
            ax.set_ylim(None, maxdval)
            if np.abs(vmin) > 0.:
                ax.set_ylim(vmin, None)

    if not isLogScale:
        if sciticks:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        if (np.isfinite(mindval) and
            np.isfinite(maxdval)):
            ax.set_ylim(mindval,maxdval)
            if maxdval-mindval < 1.0 or \
               maxdval-mindval > 100.0:
                ax.tick_params(axis='y',rotation=-35)
        ax.yaxis.get_offset_text().set_fontsize(3)

    ax.grid()

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_ylabel(dataLabel,fontsize=4)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            nlcol = np.int(np.ceil(np.sqrt(nLines)))
            lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                           framealpha=0.4,ncol=nlcol)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=3,frameon=False, \
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)


    return


###############################################################################
def plotTimeSeries2D(fig, \
                     xDates, yVals, contourVals, \
                     title="", clabel="", \
                     sciticks=False, logscale=False, signdef=False, \
                     dataLabel="y", invert_ind_axis=False, \
                     ny=1, nx=1, nplots=1, iplot=0, \
                     dmin=np.NaN, dmax=np.NaN,
                     interiorLabels=True):

# ARGUMENTS
# fig           - matplotlib figure object
# xDates        - date x-values (array of float seconds, dt.timedelta, dt.datetime)
# yVals         - second independent variable
# contourVals   - dependent variable (2d array)

# title         - subplot title, optional
# clabel        - label for dependent variable, optional
# sciticks      - whether contourVals needs scientific formatting for ticks, optional
# logscale      - whether contours are spaced logarithmically, optional, overrides sciticks
# signdef       - whether contourVals is positive/negative definite, optional
# dataLabel     - label for yVals, optional
# invert_ind_axis - whether to invert y-axis orientation, optional

# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# dmin, dmax    - min/max values of contourVals, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    xVals = pu.TDeltas2Seconds(xDates)

    # standardize c-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin,dmax,contourVals,signdef)
    if signdef:
        cmapName = 'BuPu'
        nlevs = 18

        # scientific contours
        cint = contourVals.astype(int)
        isInt = np.all((contourVals - cint) == 0)
        if isInt:
          minscid = np.nanmax(np.array([1., dmin]))
        else:
          minscid = maxdval*1.e-5
        lognorm = colors.LogNorm(vmin=minscid, vmax=maxdval)
    else:
        cmapName = 'seismic'
        nlevs = 28

        # scientific contours
        lognorm = colors.SymLogNorm(vmin=mindval, vmax=maxdval,
                    linthresh=1.e-3*maxdval, linscale=1.3, base=10)

    # plot contour
    # option 1: smoothed contours
    #cp = ax.contourf(xVals, yVals, contourVals, nlevs, cmap=cmapName, extend='both', \
    #     vmin=mindval, vmax=maxdval)

    # option 2: pixel contours
    cmap = plt.get_cmap(cmapName)
    cmap.set_bad(color = 'k', alpha = 1.0)
    if logscale:
        norm = lognorm
    else:
        levels = MaxNLocator(nbins=nlevs).tick_values(mindval,maxdval)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    xVals_pcolor, yVals_pcolor = transformXY_for_pcolor(xVals,yVals)
    cp = ax.pcolormesh(xVals_pcolor, yVals_pcolor, contourVals, cmap=cmap, norm=norm)

    #title
    ax.set_title(title,fontsize=5)

    #axes settings
    pu.format_x_for_dates(ax, xDates)
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_ylabel(dataLabel,fontsize=4)
    if interiorLabels or ix == nx-1:
        #colorbar
        m = plt.cm.ScalarMappable(cmap=cmap)
        m.set_array(contourVals)
        m.set_norm(norm)
        if (np.isfinite(mindval) and
            np.isfinite(maxdval) and
            not logscale):
            m.set_clim(mindval,maxdval)
        cb = plt.colorbar(m, ax=ax)
        #scientific formatting
        if sciticks and not logscale:
            cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        cb.ax.yaxis.get_offset_text().set_fontsize(3)

        cb.ax.tick_params(labelsize=3)
        cb.set_label(clabel,fontsize=5)

    if invert_ind_axis:
        ax.invert_yaxis()

    # optionally add a grid
    #ax.grid()

    return


###############################################################################
def transformXY_for_pcolor(xs,ys):
    # adjust centered x and y values to edges to work with pcolormesh 
    # note: works best for regularly spaced data
    xs_diff = xs[1] - xs[0]
    # extend xs by 2
    # fill in first endpoint
    xs_extend = [xs[0]-xs_diff]
    # fill in internal values
    for x in xs: xs_extend.append(x)
    # fill in last endpoint
    xs_extend.append(xs_extend[-1]+(xs[-1]-xs[-2]))
    # calculate the midpoints
    xs_pcolormesh_midpoints = []
    for ii, x in enumerate(xs_extend[:-1]):
        xs_pcolormesh_midpoints.append(x+0.5*(xs_extend[ii+1] - xs_extend[ii]))

    ys_diff = ys[1] - ys[0]
    # extend ys by 2
    # fill in first endpoint
    ys_extend = [ys[0]-ys_diff]
    # fill in internal values
    for y in ys: ys_extend.append(y)
    # fill in last endpoint
    ys_extend.append(ys_extend[-1]+(ys[-1]-ys[-2]))
    # calculate the midpoints
    ys_pcolormesh_midpoints = []
    for ii, y in enumerate(ys_extend[:-1]):
        ys_pcolormesh_midpoints.append(y+0.5*(ys_extend[ii+1] - ys_extend[ii]))

    return xs_pcolormesh_midpoints, ys_pcolormesh_midpoints


###############################################################################
lenWarnPDF = 0
nanWarnPDF = 0
def plotPDF(fig,
             countsVals, xVals,
             countsLabel,
             title="",
             indepLabel="x",
             ny=1, nx=1, nplots=1, iplot=0,
             lineAttribOffset=1,
             legend_inside=True,
             interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# countsVals       - list of arrays, each containing counts across xVals
# xVals            - independent variable on x-axis (array)
# countsLabel      - legend label for countsVals (list)

# title            - subplot title, optional
# indepLabel       - label for xVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# lineAttribOffset - offset for selecting line attributes, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    #title
    ax.set_title(title,fontsize=5)

    #add counts
    plotVals = []
    nPDFs = 0
    for ihist, countVals in enumerate(countsVals):
        if np.all(np.isnan(countVals)):
            global nanWarnPDF
            if nanWarnPDF==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+indepLabel+"; "+countsLabel[ihist])
            nanWarnPDF=nanWarnPDF+1
            continue
        if len(countVals)!=len(xVals):
            global lenWarnPDF
            if lenWarnPDF==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+indepLabel+"; "+countsLabel[ihist])
            lenWarnPDF=lenWarnPDF+1
            continue

        # Plot line for each countVals that has non-missing data

        # assume constant dx between bins
        dx = xVals[1] - xVals[0]

        ax.plot(xVals, np.divide(countVals,np.sum(countVals)*dx),
                color=pu.plotColor(len(countsVals),ihist+lineAttribOffset),
                label=countsLabel[ihist],
                ls=pu.plotLineStyle(len(countsVals),ihist+lineAttribOffset),
                linewidth=0.5)
        nPDFs = nPDFs + 1
        plotVals.append(countVals)


    if nPDFs == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add a standard normal pdf
    from scipy.stats import norm
    ax.plot(xVals, norm.pdf(xVals),
                color='k',
                ls='-',
                linewidth=0.35,
                label='N(0,1)'
    )

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)
    plt.yscale('log')
    ax.set_ylim(bottom=1.e-6)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(indepLabel,fontsize=4)
        ax.set_ylabel('PDF',fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid()

    return


###############################################################################
lenWarnRamp = 0
nanWarnRamp = 0
def plotfitRampComposite(fig,
               xVals,
               countVals,
               meanVals,
               rmsVals,
               stdVals,
               title="", dataLabel="y", \
               indepLabel="x",
               ny=1, nx=1, nplots=1, iplot=0,
               lineAttribOffset=1,
               legend_inside=True,
               interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# countVals        - Count of quantity (array)
# meanVals         - Mean of quantity (array)
# rmsVals          - RMS of quantity (array)
# stdVals          - STD of quantity (array)

# xVals           - independent variable on x-axis (array)

# title            - subplot title, optional
# dataLabel        - label for y-axis, optional
# indepLabel       - label for xVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# lineAttribOffset - offset for selecting line attributes, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)

    #title
    ax.set_title(title,fontsize=5)

    #add lines
    plotVals = []
    nLines = 0
    linesLabel = ['RMS','STD','Mean']
    for iline, lineVals in enumerate([rmsVals,stdVals,meanVals]):
        if np.all(np.isnan(lineVals)):
            global nanWarnRamp
            if nanWarnRamp==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+indepLabel+"; "+linesLabel[iline])
            nanWarnRamp=nanWarnRamp+1
            continue
        if len(lineVals)!=len(xVals):
            global lenWarnRamp
            if lenWarnRamp==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+indepLabel+"; "+linesLabel[iline])
            lenWarnRamp=lenWarnRamp+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pu.plotColor(4,iline+lineAttribOffset)

        ax.plot(xVals, lineVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pu.plotLineStyle(4,iline+lineAttribOffset),
                linewidth=0.6)
        nLines += 1
        plotVals.append(lineVals)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # Add fit for stdVals here using info from countVals
    ind0 = np.argmax(countVals)

    indexMaxX4Std = 0
    for ii, std in enumerate(stdVals):
        if np.isfinite(std): indexMaxX4Std = ii
    indexMaxX = indexMaxX4Std
    maxCount = 0
    for ii, count in enumerate(countVals):
        if count > maxCount: maxCount = count
        if count < 0.002*maxCount:
            indexMaxX = ii
            break
    if indexMaxX > indexMaxX4Std:
        ind1 = np.argmax(stdVals[0:indexMaxX4Std])
    else:
        ind1 = np.argmax(stdVals[0:indexMaxX])

    weights = [0.2]*(ind1-ind0+1)
    weights[0] = 1.0
    p = np.polyfit(xVals[ind0:ind1+1],stdVals[ind0:ind1+1],1,
                   w=weights)

    X0 = xVals[ind0]
    ERR0 = X0 * p[0] + p[1]

    # X1 = xVals[ind1]
    # ERR1 = X1 * p[0] + p[1]
    ERR1 = stdVals[ind1]
    X1 = (ERR1 - p[1]) / p[0]


    ERRfitDict = {
        'bu':{
            'X':  [round(X0,2),  round(X1,2)],
            'ERR': [round(ERR0,2), round(ERR1,2)],
        },
        'YAML':{
            'X0': [round(X0,2)],
            'X1': [round(X1,2)],
            'ERR0': [round(ERR0,2)],
            'ERR1': [round(ERR1,2)],
        },
    }

    fitX  = np.asarray([0.0]  + ERRfitDict['bu']['X']  + [xVals[indexMaxX4Std]])
    fitERR = np.asarray([ERR0] + ERRfitDict['bu']['ERR'] + [ERR1])

    plotVals.append(fitERR)

    pColor = pu.plotColor(4,1+lineAttribOffset)

    ax.plot(fitX, fitERR,
            color=pColor,
            label='Fit-STD',
            ls='-.',
            linewidth=1.2,
            marker='+',
            ms=1.5
    )

    #axes settings
    ax.xaxis.set_tick_params(labelsize=3)
    ax.yaxis.set_tick_params(labelsize=3)

    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(plotVals=plotVals)
    if (np.isfinite(mindval) and
        np.isfinite(maxdval)):
        ax.set_ylim(mindval,maxdval)

    #handle interior subplot ticks/labels
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or ix == 0:
        ax.set_xlabel(indepLabel,fontsize=4)
    if interiorLabels or iy == ny-1:
        ax.set_ylabel(dataLabel,fontsize=4)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=3,frameon=True,\
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=3,frameon=False, \
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid()

    # Add count on RHS y-axis
    ax2 = ax.twinx()
    color = 'black'
    if interiorLabels or ix == nx:
        ax2.set_ylabel('Count',fontsize=4,color=color)
    ax2.plot(xVals[:indexMaxX4Std], countVals[:indexMaxX4Std],
             color=color,
             label='Count',
             ls=':',
             linewidth=0.5)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.yaxis.set_tick_params(labelsize=3)
    plt.yscale('log')
    ax2.set_ylim(bottom=100.)

    return ERRfitDict
