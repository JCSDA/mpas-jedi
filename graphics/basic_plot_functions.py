#!/usr/bin/env python3

import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator, LongitudeLocator)
from collections.abc import Iterable
from copy import deepcopy, copy
import datetime as dt
import logging
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib as mpl
mpl.use('AGG')
import matplotlib.axes as maxes
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import plot_styles as pstyle
import plot_utils as pu
import var_utils as vu
import os


mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['hatch.linewidth'] = 0.02

_logger = logging.getLogger(__name__)

## color maps

# IR brightness temperature color map
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

btCMap = mcolors.ListedColormap(btcolors)

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

def plotDistri(lats,lons,values,
               ObsType,VarName,var_unit,out_name,nstation,levbin,
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

#Skipping the dataset if it's empty ===========================================
    if (np.isfinite(values).sum() == 0):
       print("WARNING in plotDistri: empty for", ObsType, VarName, levbin, "; skipping this dataset")
       return

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

    if cLon is not None:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(projection=ccrs.Orthographic(cLon))
    else:
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(projection=ccrs.PlateCarree())

    ax.set_global()

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

#        p = plt.pcolor(lonsPlot[lonSort], latsPlot[lonSort], valuesPlot[lonSort],
#                       transform = ccrs.PlateCarree(),
#                       cmap = cm, vmin = dmin, vmax = dmax,
#                       latlon = True, tri = True)

        p = plt.tripcolor(lonsPlot[lonSort], latsPlot[lonSort], valuesPlot[lonSort],
                       transform = ccrs.PlateCarree(),
                       cmap = cm, vmin = dmin, vmax = dmax)

    else:
        p=ax.scatter(lons[finite], lats[finite], c=values[finite],
                     transform = ccrs.PlateCarree(),
                     cmap= cm, s = dotsize)
        ax.gridlines(draw_labels=True, xlocs=np.arange(-180,180,60),linestyle='--')

    ax.coastlines()

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom",size="5%", pad=0.3,axes_class=plt.Axes)

    #fig.add_axes(cax)
    plt.colorbar(p,cax=cax,orientation='horizontal') #,cax=cax,ax=ax,orientation='horizontal')

#set title  ===================================================================
    if nstation == 0 or ObsType == 'satwind' or ObsType == 'satwnd':
        plt.text(0.5, 1.15, '%s   %s %s nlocs:%s'
            %(ObsType,VarName,var_unit,len(values[~np.isnan(values)])),
            horizontalalignment='center',
            fontsize=12, transform = ax.transAxes)
    else:
        if ObsType[:6] == 'gnssro':
            plt.text(0.5, 1.15, '%s   %s %s nlocs:%s nprofile:%s'
                %(ObsType,VarName,var_unit,len(values[~np.isnan(values)]),nstation),
                horizontalalignment='center',
                fontsize=12, transform = ax.transAxes)
        elif ObsType == 'aircraft':
            plt.text(0.5, 1.15, '%s   %s %s nlocs:%s nflight:%s'
                %(ObsType,VarName,var_unit,len(values[~np.isnan(values)]),nstation),
                horizontalalignment='center',
                fontsize=12, transform = ax.transAxes)
        else:
            plt.text(0.5, 1.15, '%s   %s %s nlocs:%s nstation:%s'
                %(ObsType,VarName,var_unit,len(values[~np.isnan(values)]),nstation),
                horizontalalignment='center',
                fontsize=12, transform = ax.transAxes)

    plt.savefig('distri_%s_%s_%s.png'%(VarName,out_name,levbin),dpi=200,bbox_inches='tight')
    plt.close()


def scatterMapFields(
  lonVals, latVals, fields,
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

  # setup map
  cLons = np.asarray([])
  lonVals_180 = {}

  for name in lonVals.keys():
    cLon = None

    # 0 < longitude <= 360
    lonVals_360 = deepcopy(lonVals[name])
    while np.max(lonVals_360) >= 360.0:
        lonVals_360[lonVals_360 >= 360.0] -= 360.0
    while np.min(lonVals_360) < 0.0:
        lonVals_360[lonVals_360 < 0.0] += 360.0

    # -180 < longitude <= 180
    lonVals_180[name] = deepcopy(lonVals_360)
    lonVals_180[name][lonVals_180[name] > 180.0] -= 360.0

    for lon in [lonVals_360, lonVals_180[name]]:
      if np.max(lon) - np.min(lon) <= 180.0:
        cLon = 0.5*(np.max(lon) + np.min(lon))

    cLons = np.append(cLons, cLon)

  anycLonNone = np.any([c is None for c in cLons])

  if anycLonNone:
    # plot entire Earth
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(projection=ccrs.Mollweide(0.0))

  else:
    # plot single projected side of Earth
    cLon = cLons[0]
    if cLon > 180.0: cLon-=360.0
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(projection=ccrs.Orthographic(cLon))

  assert (cbarType is None or cbarType in ['Log', 'SymLog']), \
    'scatterMapFields: invalid cbarType: '+cbarType

  for name, field in fields.items():
    f = c=c.get(name, field)
    finite = np.isfinite(f)
    lons = lonVals_180[name][finite]
    lats = latVals[name][finite]
    f = f[finite]

    ## transform to pcolormesh and cartopy conventions
    # longitude monotonically increasing
    lonSort = np.argsort(lons)
    lons = lons[lonSort]
    lats = lats[lonSort]
    f = f[lonSort]

    if dmin is None:
       vmin = f.min()
    else:
       vmin = dmin
    if dmax is None:
       vmax = f.max()
    else:
       vmax = dmax

    if cbarType is None:
      norm = None
    elif cbarType == 'Log':
      if vmin <= logVLim: vmin = logVLim
      f[f < vmin] = vmin
      norm=mcolors.LogNorm(vmin=vmin, vmax=vmax)
    elif cbarType == 'SymLog':
      norm=mcolors.SymLogNorm(vmin=vmin, vmax=vmax,
                             linthresh=1.e-4*vmax, linscale=1.0, base=10)

    sc = ax.scatter(lons, lats, c=f,
                   s = sizes.get(name, 1),
                   cmap = cmap,
                   norm = norm,
                   marker = markers.get(name, '.'), linewidth = 0,
                   transform=ccrs.PlateCarree(),
    )

  # show full projection extent
  ax.set_global()

  # add coastlines
  ax.coastlines()

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("bottom",size="5%", pad=0.3,axes_class=plt.Axes)
  cb = plt.colorbar(sc, cax=cax, orientation='horizontal')

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
    # yonggangyu introduce epsilon and xi for plotting absolutely zero field,
    # solving vmin, vcenter, vmax ascending order issue
    epsilon  = 1.e-8
    if (valuemin > 0 or valuemax < 0):
        color = 'rainbow'
        plt.contourf(xarray,ylevels,Stats,40,vmin=valuemin, vmax=valuemax,cmap=color)
        xi=-1
    else:
        cmap = 'coolwarm'
        if ( -valuemin < epsilon and valuemax < epsilon ):
            xi=1
            valuemin = -epsilon
            valuemax =  epsilon
        elif ( -valuemin < epsilon and valuemax > epsilon ):
            xi=2
            valuemin = -epsilon
        elif ( -valuemin > epsilon and valuemax < epsilon ):
            xi=3
            valuemax =  epsilon
        else:
            xi=4
        #print('xi= '+str(xi)+' valuemin= ',str(valuemin)+' valuemax= ',str(valuemax))
        norm = mcolors.DivergingNorm(vmin=valuemin, vcenter=0, vmax=valuemax)
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


## default independent axis configurations
defaultIndepConfig = {'transform': None, 'invert': False}

## Functions used to transform axes
axisTransforms = {}

# Pressure: x**(1/pressureExponent)
pressureExponent = 3.0
def forwardPressure(x):
    global pressureExponent
    return np.power(x, 1./pressureExponent)

def inversePressure(y):
    global pressureExponent
    return np.power(y, pressureExponent)

axisTransforms['Pressure'] = {
    'forward': forwardPressure,
    'inverse': inversePressure,
    'locator': mticker.FixedLocator,
}

#ciExponent = np.e
ciExponent = 2.0
def forwardCloudImpact(x):
    global ciExponent
    return np.sign(x) * np.power(np.abs(x), 1./ciExponent)

def inverseCloudImpact(y):
    global ciExponent
    return np.sign(y) * np.power(np.abs(y), ciExponent)

axisTransforms['CloudImpact'] = {
    'forward': forwardCloudImpact,
    'inverse': inverseCloudImpact,
    'locator': mticker.FixedLocator,
}

class axisTransform:
    '''
    Accessor class for axisTransforms dictionary of axis transformation and locator functions
    '''
    def __init__(self, name):
        assert name in axisTransforms, ('ERROR in axisTransform: ',name,' is not a defined transformation')
        assert 'forward' in axisTransforms[name], ('ERROR in axisTransform: ',name,' is missing forward transform')
        assert 'inverse' in axisTransforms[name], ('ERROR in axisTransform: ',name,' is missing inverse transform')
        assert 'locator' in axisTransforms[name], ('ERROR in axisTransform: ',name,' is missing locator')
        self.tname = name

    def forward(self):
        return axisTransforms[self.tname]['forward']

    def inverse(self):
        return axisTransforms[self.tname]['inverse']

    def locator(self):
        return axisTransforms[self.tname]['locator']


class OneCenteredAxisTicks():
  '''
  produces axis locators and formatters for a one-centered axis, e.g., for ratios
  '''
  fineBaseRatios = np.array([2., 3., 5., 8., 10.]) #Fibonacci-like sequence
  coarseBaseRatios = np.array([2., 5., 10.]) #coarser sequence

  nDecades = 1
  def __init__(self, maxLim):
    maxDecade = np.log10(np.max(np.array([maxLim, 1.00001]))-1.0)
    maxDecade = int(np.round(maxDecade))

    allDecades = np.arange(maxDecade-self.nDecades, maxDecade+1)

    nDecimalDigits = str(int(np.max(np.array([1, 1+np.abs(allDecades[-1])]))))
    self.majorFormatter = mticker.FormatStrFormatter('%.'+nDecimalDigits+'f')
    self.minorFormatter = mticker.NullFormatter()

    def baseRatios(decade):
      if decade >= 0:
        return self.fineBaseRatios
      else:
        return self.coarseBaseRatios

    def decadeRatios(decade):
      if decade >= 0:
        return list(baseRatios(decade) * (10.**decade))
      else:
        return list(1. + baseRatios(decade) * (10.**decade))

    #MajorLocator
    ratiosUpper = np.array([])
    majorRatiosUpper = []
    for decade in allDecades:
      majorRatiosUpper += decadeRatios(decade)

    ratiosUpper = np.array(majorRatiosUpper)
    ratiosLower = np.flip(np.divide(1.,ratiosUpper))
    ratios = np.append(ratiosLower, np.array([1.]))
    ratios = np.append(ratios, ratiosUpper)
    self.majorLocator = mticker.FixedLocator(ratios)

    #MinorLocator - fractally divides baseRatios(minorDecade) using a 10x reduction relative to
    #               minorDecade
    minorDecade = allDecades[0]
    minorUpper = []
    lastRatio = 1.
    for thisRatio in decadeRatios(minorDecade):
      diff = thisRatio - lastRatio
      for r in 0.1*baseRatios(minorDecade):
        minorUpper.append(lastRatio + r*diff)
      lastRatio = thisRatio

    ratiosUpper = np.array(minorUpper)
    ratiosLower = np.flip(np.divide(1.,ratiosUpper))
    ratios = np.append(ratiosLower, np.array([1.]))
    ratios = np.append(ratios, ratiosUpper)
    self.minorLocator = mticker.FixedLocator(ratios)

  def setTicks(self, axis):
    axis.set_major_locator(self.majorLocator)
    axis.set_minor_locator(self.minorLocator)
    axis.set_major_formatter(self.majorFormatter)
    axis.set_minor_formatter(self.minorFormatter)


# maximum number of legend entries for line-style plots
maxLegendEntries = 12

# common line width for "ax.plot"
commonLineWidth = 0.80
significanceLineWidth = 0.3
errorbarLineWidth = 0.38
gridLineWidth = 0.05
coastLineWidth = 0.15

titleFontSize = 5.5
cbarLabelFontSize = 4.6
axisLabelFontSize = 4.8
legendLabelFontSize1 = 3.8
legendLabelFontSize2 = 3.2
categoryLabelFontSize = 3.8

###############################################################################
lenWarnSer = 0
nanWarnSer = 0
def plotSeries(fig,
               linesVals, xVals,
               linesLabel,
               title="", indepLabel=None, dataLabel=None,
               indepConfig=defaultIndepConfig,
               sciTicks=False, logScale= False, centralValue=None,
               ny=1, nx=1, nplots=1, iplot=0,
               linesValsMinCI=None, linesValsMaxCI=None,
               dmin=np.NaN, dmax=np.NaN,
               lineAttribOffset=0,
               legend_inside=True,
               interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# xVals            - independent variable on x-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# indepLabel       - label for xVals, optional
# dataLabel        - label for linesVals, optional

# sciTicks         - whether linesVals needs scientific formatting for ticks, optional
# logScale         - y-axis is scaled logarithmically, optional, overrides sciTicks
# centralValue     - central value of linesVals, optional
# indepConfig      - x axis formatting config, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum confidence interval bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum confidence interval bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional
# interiorLabels   - whether to add titles and axis labels for interior subplots, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

    #add lines
    plotVals = np.asarray([])
    nLines = 0
    for iline, lineVals in enumerate(linesVals):
        if np.all(np.isnan(lineVals)):
            global nanWarnSer
            if nanWarnSer==0:
                _logger.warning("skipping all-NaN data")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            nanWarnSer=nanWarnSer+1
            continue
        if len(lineVals)!=len(xVals):
            global lenWarnSer
            if lenWarnSer==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+dataLabel+"; "+linesLabel[iline])
            lenWarnSer=lenWarnSer+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(xVals, lineVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth)
        nLines += 1
        plotVals = np.append(plotVals, lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            lineArr = np.array(lineVals)
            xArr = np.array(xVals)

            # user error bars when there is more than one experiment
            if len(plotVals) > 1:
              deltaError = np.asarray([
                -(linesValsMinCI[iline]-lineArr),
                linesValsMaxCI[iline]-lineArr])

              ax.errorbar(xVals, lineArr, yerr=deltaError,
                          fmt = 'none',
                          ecolor=pColor,
                          elinewidth=errorbarLineWidth,
                          capsize=1.5, capthick=0.2)

            # otherwise use shaded significance region
            else:
              # test statistical significance versus centralValue
              if centralValue is None:
                  isSignificant = np.empty(len(lineVals))
                  isSignificant[:] = np.NaN
                  centralValue_ = 0.0
              else:
                  isSignificant = np.multiply(np.subtract(linesValsMinCI[iline], centralValue),
                                            np.subtract(linesValsMaxCI[iline], centralValue))
                  centralValue_ = centralValue
              isSignificant = np.array([x if np.isfinite(x) else -1.0 for x in isSignificant])

              negsiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] < centralValue_)],dtype=int)
              if len(negsiginds) > 0:
                  ax.plot(xArr[negsiginds], lineArr[negsiginds],
                          color=pColor,
                          ls='',
                          marker='v',
                          markersize=0.2)

              possiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] > centralValue_)],dtype=int)
              if len(possiginds) > 0:
                  ax.plot(xArr[possiginds], lineArr[possiginds],
                          color=pColor,
                          ls='',
                          marker='^',
                          markersize=0.2)

              ax.plot(xVals, linesValsMinCI[iline],
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)
              ax.plot(xVals, linesValsMaxCI[iline],
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)

              ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                              color=pColor,
                              edgecolor=pColor,
                              linewidth=0.0, alpha = 0.1)
              ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                              where=isSignificant > 0.0,
                              color=pColor,
                              edgecolor=pColor,
                              linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add horizontal centralValue line for unbounded quantities
    if centralValue is not None:
        ax.plot([xVals[0], xVals[-1]], [centralValue, centralValue], ls="--", c=".3",
            linewidth=0.7,markersize=0)

    # standardize data-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, plotVals, centralValue)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    #axes settings
    scaleType = 'linear'

    # TODO: encapsulate this log-axis formatting and oneCenteredRatios tick formatting
    #       into a reusable function for all plot types, including color axis in 2D plots
    if logScale:
        finite = np.isfinite(plotVals)
        nonzero = finite
        nonzero[finite] = np.greater(np.abs(plotVals[finite]), 0.)
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if centralValue is None:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0 and vmax > 0.:
                    mindval, maxdval = pu.get_clean_ax_limits(
                        np.log10(np.nanmax(np.array([vmin, vmax*1.e-3]))),
                        np.log10(vmax))
                    mindval = np.power(10., mindval)
                    maxdval = np.power(10., maxdval)
                    scaleType = 'log'
                    maxp = None
                    if np.isfinite(maxdval):
                        maxp = maxdval
                    pint = plotVals.astype(int)
                    isInt = np.all((plotVals - pint) == 0)
                    if isInt:
                        minp = np.nanmax(np.array([1., mindval]))
                    else:
                        minp = mindval
            elif float(centralValue) == 0.:
                scaleType = 'symlog'
                if maxp is not None: minp = -maxp
            else:
                scaleType = 'log'
                maxratio = np.max([vmax/centralValue, centralValue/vmin])
                maxp = maxratio / centralValue
                p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
                maxratio = maxp / centralValue
                minp = centralValue / maxratio
                if minp == maxp:
                    minp = dmin
                    maxp = dmax

    # for non-linear axis scales, limits must be set before scale
    ax.set_ylim(minp, maxp)
    ax.set_yscale(scaleType)

    if centralValue is not None and float(centralValue) == 1. and scaleType == 'log':
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(ax.yaxis)

    if scaleType == 'linear':
        if sciTicks:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    if indepConfig['transform'] is None:
        ax.set_xscale('linear')
        ax.xaxis.set_major_locator(mticker.AutoLocator())
    elif indepConfig['transform'] == 'logit':
        ax.set_xscale('logit', nonpositive='clip')
        ax.xaxis.set_major_locator(mticker.AutoLocator())
    else:
        transform = axisTransform(indepConfig['transform'])
        locator = transform.locator()
        ax.set_xscale('function', functions=(transform.forward(), transform.inverse()))
        ax.xaxis.set_major_locator(locator(xVals))
        ax.tick_params(axis='x', rotation=60.)

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if indepLabel is not None: ax.set_xlabel(indepLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if dataLabel is not None: ax.set_ylabel(dataLabel,fontsize=axisLabelFontSize)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                           framealpha=0.4,ncol=1)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if indepConfig['invert']:
        ax.invert_xaxis()

    ax.grid(linewidth=gridLineWidth)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return

###############################################################################
lenWarnProf = 0
nanWarnProf = 0
def plotProfile(fig,
                linesVals, yVals,
                linesLabel,
                title="", indepLabel=None, dataLabel=None,
                indepConfig=defaultIndepConfig,
                sciTicks=False, logScale=False, centralValue=None,
                ny=1, nx=1, nplots=1, iplot=0,
                linesValsMinCI=None, linesValsMaxCI=None,
                dmin=np.NaN, dmax=np.NaN,
                lineAttribOffset=0,
                legend_inside=True,
                interiorLabels=True,
                marker=None):

# ARGUMENTS
# fig              - matplotlib figure object
# linesVals        - dependent variable (list of arrays)
# yVals            - independent variable on y-axis (array)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# indepLabel       - label for yVals, optional
# dataLabel        - label for linesVals, optional
# indepConfig      - y axis formatting config, optional

# sciTicks         - whether linesVals needs scientific formatting for ticks, optional
# logScale         - x-axis is scaled logarithmically, optional, overrides sciTicks
# centralValue     - central value of linesVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum confidence interval bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum confidence interval bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional
# interiorLabels   - whether to add titles and axis labels for interior subplots, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

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
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(lineVals, yVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth,
                marker=marker)
        nLines += 1
        plotVals = np.append(plotVals,lineVals)

        # Add shaded error regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            lineArr = np.array(lineVals)
            yArr = np.array(yVals)

            # user error bars when there is more than one experiment
            if len(plotVals) > 1:
              deltaError = np.asarray([
                -(linesValsMinCI[iline]-lineArr),
                linesValsMaxCI[iline]-lineArr])

              ax.errorbar(lineArr, yVals, xerr=deltaError,
                          fmt = 'none',
                          ecolor=pColor,
                          elinewidth=errorbarLineWidth,
                          capsize=1.5, capthick=0.2)

            # otherwise use shaded significance region
            else:
              # test statistical significance versus centralValue
              if centralValue is None:
                  isSignificant = np.empty(len(lineVals))
                  isSignificant[:] = np.NaN
                  centralValue_ = 0.0
              else:
                  isSignificant = np.multiply(np.subtract(linesValsMinCI[iline], centralValue),
                                            np.subtract(linesValsMaxCI[iline], centralValue))
                  centralValue_ = centralValue
              isSignificant = np.array([x if np.isfinite(x) else -1.0 for x in isSignificant])

              negsiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] < centralValue_)],dtype=int)
              if len(negsiginds) > 0:
                  ax.plot(lineArr[negsiginds], yArr[negsiginds],
                          color=pColor,
                          ls='',
                          marker='<',
                          markersize=0.2)

              possiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] > centralValue_)],dtype=int)
              if len(possiginds) > 0:
                  ax.plot(lineArr[possiginds], yArr[possiginds],
                          color=pColor,
                          ls='',
                          marker='>',
                          markersize=0.2)

              ax.plot(linesValsMinCI[iline], yVals,
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)
              ax.plot(linesValsMaxCI[iline], yVals,
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)

              ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                               color=pColor,
                               edgecolor=pColor,
                               linewidth=0.0, alpha = 0.1)
              ax.fill_betweenx(yVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                               where=isSignificant > 0.0,
                               color=pColor,
                               edgecolor=pColor,
                               linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add vertical centralValue line for unbounded quantities
    if centralValue is not None:
        ax.plot([centralValue, centralValue], [yVals[0], yVals[-1]], ls="--", c=".3",
            linewidth=0.7,markersize=0)

    # standardize data-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, plotVals, centralValue)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    #axes settings
    scaleType = 'linear'
    if logScale:
        finite = np.isfinite(plotVals)
        nonzero = finite
        nonzero[finite] = np.greater(np.abs(plotVals[finite]), 0.)
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if centralValue is None:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0 and vmax > 0.:
                    mindval, maxdval = pu.get_clean_ax_limits(
                        np.log10(np.nanmax(np.array([vmin, vmax*1.e-3]))),
                        np.log10(vmax))
                    mindval = np.power(10., mindval)
                    maxdval = np.power(10., maxdval)
                    scaleType = 'log'
                    maxp = None
                    if np.isfinite(maxdval):
                        maxp = maxdval
                    pint = plotVals.astype(int)
                    isInt = np.all((plotVals - pint) == 0)
                    if isInt:
                        minp = np.nanmax(np.array([1., mindval]))
                    else:
                        minp = mindval
            elif float(centralValue) == 0.:
                scaleType = 'symlog'
                if maxp is not None: minp = -maxp
            else:
                scaleType = 'log'
                maxratio = np.max([vmax/centralValue, centralValue/vmin])
                maxp = maxratio / centralValue
                p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
                maxratio = maxp / centralValue
                minp = centralValue / maxratio
                if minp == maxp:
                    minp = dmin
                    maxp = dmax

    # for non-linear axis scales, limits must be set before scale
    ax.set_xlim(minp, maxp)
    ax.set_xscale(scaleType)

    if centralValue is not None and float(centralValue) == 1. and scaleType == 'log':
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(ax.xaxis)
        ax.tick_params(axis='x', rotation=90)

    if scaleType == 'linear':
        if sciTicks:
            ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        if (np.isfinite(minp) and
            np.isfinite(maxp)):
            if maxp-minp < 1.0 or \
               maxp-minp > 100.0:
                ax.tick_params(axis='x', rotation=90)
        ax.xaxis.get_offset_text().set_fontsize(3)

    if indepConfig['transform'] is None:
        ax.set_yscale('linear')
        ax.yaxis.set_major_locator(mticker.AutoLocator())
    elif indepConfig['transform'] == 'logit':
        ax.set_yscale('logit', nonpositive='clip')
        ax.yaxis.set_major_locator(mticker.AutoLocator())
    else:
        transform = axisTransform(indepConfig['transform'])
        locator = transform.locator()
        ax.set_yscale('function', functions=(transform.forward(), transform.inverse()))
        ax.yaxis.set_major_locator(locator(yVals))

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if dataLabel is not None: ax.set_xlabel(dataLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if indepLabel is not None: ax.set_ylabel(indepLabel,fontsize=axisLabelFontSize)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                           framealpha=0.4,ncol=1)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)

    if indepConfig['invert']:
        ax.invert_yaxis()

    ax.grid(linewidth=gridLineWidth)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return


###############################################################################
lenWarnTS=0
nanWarnTS=0
def plotTimeSeries(fig,
                   xsDates, linesVals,
                   linesLabel,
                   title="", indepLabel=None, dataLabel=None,
                   sciTicks=False, logScale = False, centralValue=None,
                   ny=1, nx=1, nplots=1, iplot=0,
                   linesValsMinCI=None, linesValsMaxCI=None,
                   dmin=np.NaN, dmax=np.NaN,
                   lineAttribOffset=0,
                   legend_inside=True,
                   interiorLabels=True):

# ARGUMENTS
# fig              - matplotlib figure object
# xsDates          - date x-values (list/array or list of lists/arrays
#                                   of float seconds, datetime.timedelta, datetime.datetime)
# linesVals        - dependent variable (list of arrays)
# linesLabel       - legend label for linesVals (list)

# title            - subplot title, optional
# dataLabel        - label for linesVals, optional
# sciTicks         - whether linesVals needs scientific formatting for ticks, optional
# logScale         - y-axis is scaled logarithmically, optional, overrides sciTicks
# centralValue     - central value of linesVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# linesValsMinCI   - minimum confidence interval bound for linesVals (list of arrays), optional
# linesValsMaxCI   - maximum confidence interval bound for linesVals (list of arrays), optional
# Note: linesValsMinCI and linesValsMaxCI must be specified together

# lineAttribOffset - offset for selecting line attributes, optional
# dmin, dmax       - min/max values of linesVals, optional
# legend_inside    - whether legend should be placed inside the subplot, optional
# interiorLabels   - whether to add titles and axis labels for interior subplots, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

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
        if isinstance(xsDates[0], Iterable):
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
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(xVals, lineVals,
                label=linesLabel[iline],
                color=pColor,
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth)
        nLines += 1
        plotVals = np.append(plotVals, lineVals)

        # Add shaded CI regions if specified
        if linesValsMinCI is not None and \
           linesValsMaxCI is not None:

            lineArr = np.array(lineVals)
            xArr = np.array(xVals)

            # user error bars when there is more than one experiment
            if len(plotVals) > 1:
              deltaError = np.asarray([
                -(linesValsMinCI[iline]-lineArr),
                linesValsMaxCI[iline]-lineArr])

              ax.errorbar(xVals, lineArr, yerr=deltaError,
                          fmt = 'none',
                          ecolor=pColor,
                          elinewidth=errorbarLineWidth,
                          capsize=1.5, capthick=0.2)

            # otherwise use shaded significance region
            else:
              # test statistical significance versus centralValue
              if centralValue is None:
                  isSignificant = np.empty(len(lineVals))
                  isSignificant[:] = np.NaN
                  centralValue_ = 0.0
              else:
                  isSignificant = np.multiply(np.subtract(linesValsMinCI[iline], centralValue),
                                            np.subtract(linesValsMaxCI[iline], centralValue))
                  centralValue_ = centralValue

              isSignificant = np.array([x if np.isfinite(x) else -1.0 for x in isSignificant])

              negsiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] < centralValue_)],dtype=int)
              if len(negsiginds) > 0:
                  ax.plot(xArr[negsiginds], lineArr[negsiginds],
                          color=pColor,
                          ls='',
                          marker='v',
                          markersize=0.8)

              possiginds = np.array([i for i,x in enumerate(isSignificant)
                                     if (x > 0.0 and lineArr[i] > centralValue_)],dtype=int)
              if len(possiginds) > 0:
                  ax.plot(xArr[possiginds], lineArr[possiginds],
                          color=pColor,
                          ls='',
                          marker='^',
                          markersize=0.8)

              ax.plot(xVals, linesValsMinCI[iline],
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)
              ax.plot(xVals, linesValsMaxCI[iline],
                      color=pColor,
                      alpha=0.4,
                      ls='-',
                      linewidth=significanceLineWidth)

              ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                              color=pColor,
                              edgecolor=pColor,
                              linewidth=0.0, alpha = 0.1)
              ax.fill_between(xVals, linesValsMinCI[iline], linesValsMaxCI[iline],
                              where=isSignificant > 0.0,
                              color=pColor,
                              edgecolor=pColor,
                              linewidth=0.2, alpha = 0.3)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add horizontal centralValue line for unbounded quantities
    if centralValue is not None:
        ax.plot([minX, maxX], [centralValue, centralValue], ls="--", c=".3",
            linewidth=0.7,markersize=0)

    #axes settings
    if isinstance(xsDates[0], Iterable):
        pu.format_x_for_dates(ax, xsDates[0])
    else:
        pu.format_x_for_dates(ax, xsDates)

    # standardize data-limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, plotVals, centralValue)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    #axes settings
    scaleType = 'linear'
    if logScale:
        finite = np.isfinite(plotVals)
        nonzero = finite
        nonzero[finite] = np.greater(np.abs(plotVals[finite]), 0.)
        if nonzero.sum() > 0:
            vmin = np.nanmin(np.abs(plotVals[nonzero]))
            vmax = np.nanmax(np.abs(plotVals[nonzero]))
            if centralValue is None:
                # log tick labels look bad for single decade
                if vmax / vmin > 10.0 and vmax > 0.:
                    mindval, maxdval = pu.get_clean_ax_limits(
                        np.log10(np.nanmax(np.array([vmin, vmax*1.e-3]))),
                        np.log10(vmax))
                    mindval = np.power(10., mindval)
                    maxdval = np.power(10., maxdval)
                    scaleType = 'log'
                    maxp = None
                    if np.isfinite(maxdval):
                        maxp = maxdval
                    pint = plotVals.astype(int)
                    isInt = np.all((plotVals - pint) == 0)
                    if isInt:
                        minp = np.nanmax(np.array([1., mindval]))
                    else:
                        minp = mindval
            elif float(centralValue) == 0.:
                scaleType = 'symlog'
                if maxp is not None: minp = -maxp
            else:
                scaleType = 'log'
                maxratio = np.max([vmax/centralValue, centralValue/vmin])
                maxp = maxratio / centralValue
                p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
                maxratio = maxp / centralValue
                minp = centralValue / maxratio
                if minp == maxp:
                    minp = dmin
                    maxp = dmax

    # for non-linear axis scales, limits must be set before scale
    ax.set_ylim(minp, maxp)
    ax.set_yscale(scaleType)

    if centralValue is not None and float(centralValue) == 1. and scaleType == 'log':
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(ax.yaxis)

    if scaleType == 'linear':
        if sciTicks:
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.yaxis.get_offset_text().set_fontsize(3)

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    ax.grid(linewidth=gridLineWidth)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if indepLabel is not None: ax.set_xlabel(indepLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if dataLabel is not None: ax.set_ylabel(dataLabel,fontsize=axisLabelFontSize)

    #legend
    if nLines <= maxLegendEntries:
        if legend_inside:
            #INSIDE AXES
            nlcol = np.int(np.ceil(np.sqrt(nLines)))
            lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                           framealpha=0.4,ncol=nlcol)
            lh.get_frame().set_linewidth(0.0)
        elif ix==nx-1 or iplot==nplots-1:
            #OUTSIDE AXES
            ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                      bbox_to_anchor=(1.02, 1), borderaxespad=0)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return


###############################################################################
def scoreCard(fig,
           xVals, binNames, contourVals,
           title='', xLabel=None, yLabel=None, cLabel=None,
           xConfig=defaultIndepConfig,
           sciTicks=False, logScale=False, centralValue=None,
           ny=1, nx=1, nplots=1, iplot=0,
           contourValsMinCI=None, contourValsMaxCI=None,
           dmin=np.NaN, dmax=np.NaN,
           interiorLabels=True):

# ARGUMENTS
# fig           - matplotlib figure object
# xVals         - x-values (list/array of float, datetime.datetime, or datetime.timedelta)
# binNames      - bin names to display on y axis
# contourVals   - dependent variable (2d numpy array), len(binNames) x len(xVals)

# title         - subplot title, optional
# xLabel        - label for xVals, optional
# yLabel        - label for yVals, optional
# cLabel        - label for dependent variable, optional
# xConfig       - x axis formatting config, optional

# log_y_axis    - whether y-axis is scaled logarithmically, optional
# sciTicks      - whether contourVals needs scientific formatting for ticks, optional
# logScale      - whether contours are spaced logarithmically, optional, overrides sciTicks
# centralValue  - central value of contourVals, optional

# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# contourValsMinCI   - minimum confidence interval bound for contourVals (2d numpy array), optional
# contourValsMaxCI   - maximum confidence interval bound for contourVals (2d numpy array), optional
# Note: contourValsMinCI and contourValsMaxCI must be specified together

# dmin, dmax    - min/max values of contourVals, optional
# interiorLabels- whether to add titles, axis, and colorbar labels for interior subplots, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # standardize color limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, contourVals, centralValue, buffr=0.05)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    # default cAxisType
    cAxisType = 'sequential'

    if centralValue is not None and pu.isfloat(centralValue):
        cAxisType = 'diverging'

    # log contours
    isLogScale = logScale
    finite = np.isfinite(contourVals)
    nonzero = finite
    nonzero[finite] = np.greater(np.abs(contourVals[finite]), 0.)
    if nonzero.sum()==0:
        isLogScale = False

    elif isLogScale:
        if centralValue is None:
            minp_, maxp_ = pu.get_clean_ax_limits(
                np.log10(np.nanmax(np.array([dmin, dmax*1.e-3]))),
                np.log10(dmax), np.log10(contourVals[nonzero]), buffr=0.05)
            minp_ = np.power(10., minp_)
            maxp_ = np.power(10., maxp_)
            if maxp_ / minp_ < 20.:
                isLogScale = False
            else:
                minp = minp_
                maxp = maxp_
                cint = contourVals.astype(int)
                isInt = np.all((contourVals - cint) == 0)
                if isInt:
                    minp = np.nanmax(np.array([1., minp]))
                lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)
        elif float(centralValue) == 0.:
            lognorm = mcolors.SymLogNorm(vmin=minp, vmax=maxp,
                        linthresh=1.e-3*maxp, linscale=1.3, base=10)
        else:
            vmin = np.nanmin(np.abs(contourVals[nonzero]))
            vmax = np.nanmax(np.abs(contourVals[nonzero]))

            maxratio = np.max([vmax/centralValue, centralValue/vmin])
            maxp = maxratio / centralValue
            p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
            maxratio = maxp / centralValue
            minp = centralValue / maxratio

            lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)

    cmap = copy(pstyle.cmaps[cAxisType]['map'])
    #cmap.set_bad(color='gray', alpha=0.5)
    if isLogScale:
        norm = lognorm
    #elif cAxisType == 'diverging':
    #    norm = mcolors.CenteredNorm(centralValue)
    else:
        levels = mticker.MaxNLocator(nbins=pstyle.cmaps[cAxisType]['n']).tick_values(minp, maxp)
        norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # remove rows filled with missing values
    keep = []
    for iy in range(len(binNames)):
      if any(np.isfinite(contourVals[iy,:])):
        keep.append(iy)
    binNames_ = np.array(binNames)[keep]
    contourVals_ = contourVals[keep,:]

    # (x, y) specify corners instead of centers
    yBins = np.arange(len(binNames_))
    xBins = np.arange(len(xVals))
    xBins_edges, yBins_edges = transformXY_to_edges(xBins, yBins)

    pcm = ax.pcolormesh(xBins_edges, yBins_edges, contourVals_, cmap=cmap, norm=norm)

    ax.set_yticks(yBins)
    ax.set_yticklabels(binNames_)

    if isinstance(xVals[0], dt.timedelta):
      xVals_ = pu.TDeltas2Seconds(xVals)
      xTicks = pu.timeDeltaTicks(xVals)
      xBinsTicks = [xVals_.index(t) for t in xTicks]
      ax.set_xticks(xBinsTicks)
      ax.set_xticklabels([pu.timeDeltaTickLabels(t, None) for t in xTicks])
      ax.xaxis.set_tick_params(rotation=90)
      ax.set_aspect('equal')
    else:
      # TODO: clean this up if non-timedelta x values are needed
      ax.set_xticks(xBins)
      ax.set_xticklabels(xVals)

    # Add significance regions if specified
    if (contourValsMinCI is not None and
       np.isfinite(contourValsMinCI).sum() > 0 and
       contourValsMaxCI is not None and
       np.isfinite(contourValsMaxCI).sum() > 0 ):

        contourValsMinCI_ = contourValsMinCI[keep,:]
        contourValsMaxCI_ = contourValsMaxCI[keep,:]

        # test statistical significance versus centralValue
        if centralValue is not None:
            isSignificant = np.multiply(np.subtract(contourValsMinCI_, centralValue),
                                     np.subtract(contourValsMaxCI_, centralValue))
            centralValue_ = centralValue
            isSignificant[~np.isfinite(isSignificant)] = -1.0

            # Create patches indicating significance
            w100 = xBins_edges[1:] - xBins_edges[0:-1]
            h100 = yBins_edges[1:] - yBins_edges[0:-1]
            patches = []

            q = 97.73 # 2 sigma, or get quartile from argument in place of CI
            quartile = np.full_like(isSignificant, q)
            distortion = 6. # best for triangles
            cellFraction = np.power(np.abs(quartile - 50.) / 50., distortion) # distort patch with large power
            cellFraction[isSignificant<=0.] = 0.

            for ix in range(len(xBins)):
              for iy in range(len(yBins)):
                if cellFraction[iy, ix] > 0.:
                  # patch center
                  xc = xBins[ix]
                  yc = yBins[iy]

                  # patch width and height
                  w = w100[ix] * cellFraction[iy, ix]
                  h = h100[iy] * cellFraction[iy, ix]

                  # add triangle patch
                  if contourValsMinCI_[iy, ix] > centralValue:
                    patches.append(
                      mpatches.RegularPolygon((xc, yc+h/8.), 3, radius=w/2.,
                        orientation = np.pi, fill = False)
                    )
                  if contourValsMaxCI_[iy, ix] < centralValue:
                    patches.append(
                      mpatches.RegularPolygon((xc, yc-h/8.), 3, radius=w/2.,
                        orientation = 0., fill = False)
                    )

            # Create patch collection with specified color/alpha
            pc = mcollections.PatchCollection(patches,
                edgecolor = 'black',
                facecolor = 'none',
                snap = False,
                alpha = 0.5,
                linewidth = 0.05,
                linestyles = (0, (5, 1)), # (offset, onoffseq), where onoffseq is an even length tuple of on and off ink lengths in points
            )

            # Add collection to axes
            ax.add_collection(pc)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

    #axes settings
    if xConfig['transform'] is not None:
        if xConfig['transform'] == 'logit':
            ax.set_xscale('logit', nonpositive='clip')
            ax.xaxis.set_major_locator(mticker.AutoLocator())
        else:
            transform = axisTransform(xConfig['transform'])
            locator = transform.locator()
            ax.set_xscale('function', functions=(transform.forward(), transform.inverse()))
            ax.xaxis.set_major_locator(locator(xVals))
            ax.tick_params(axis='x', rotation=90)

    ax.tick_params(axis='x', which='major', labelsize=3)
    ax.tick_params(axis='x', which='minor', labelsize=2)
    ax.tick_params(axis='y', which='major', labelsize=categoryLabelFontSize)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if interiorLabels or iy == ny-1:
        if xLabel is not None: ax.set_xlabel(xLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if yLabel is not None: ax.set_ylabel(yLabel,fontsize=axisLabelFontSize)

    #colorbar
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(contourVals)
    m.set_norm(norm)
    if not isLogScale:
        m.set_clim(minp, maxp)

    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes("right", size="5%", pad=0.05)
    f = ax.get_figure()
    f.add_axes(ax_cb, label='cb')
    cb = plt.colorbar(m, cax=ax_cb)

    #scientific formatting
    if sciTicks and not logScale:
        cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    cb.ax.yaxis.get_offset_text().set_fontsize(3)

    if centralValue is not None and float(centralValue) == 1. and isLogScale:
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(cb.ax.yaxis)

    cb.ax.tick_params(axis='y', which='major', labelsize=3)
    cb.ax.tick_params(axis='y', which='minor', labelsize=2)

    if cLabel is not None: cb.set_label(cLabel,fontsize=cbarLabelFontSize)

    if xConfig['invert']:
        ax.invert_xaxis()

    # optionally add a grid
    #ax.grid(linewidth=gridLineWidth)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return


###############################################################################
def plot2D(fig,
           xVals, yVals, contourVals,
           title='', xLabel=None, yLabel=None, cLabel=None,
           xConfig=defaultIndepConfig,
           yConfig=defaultIndepConfig,
           sciTicks=False, logScale=False, centralValue=None,
           ny=1, nx=1, nplots=1, iplot=0,
           contourValsMinCI=None, contourValsMaxCI=None,
           dmin=np.NaN, dmax=np.NaN,
           interiorLabels=True):

# ARGUMENTS
# fig           - matplotlib figure object
# xVals         - x-values (list/array of float, datetime.datetime, or datetime.timedelta)
# yVals         - y-values
# contourVals   - dependent variable (2d numpy array), len(yVals) x len(xVals)

# title         - subplot title, optional
# xLabel        - label for xVals, optional
# yLabel        - label for yVals, optional
# cLabel        - label for dependent variable, optional
# xConfig       - x axis formatting config, optional
# yConfig       - y axis formatting config, optional

# log_y_axis    - whether y-axis is scaled logarithmically, optional
# sciTicks      - whether contourVals needs scientific formatting for ticks, optional
# logScale      - whether contours are spaced logarithmically, optional, overrides sciTicks
# centralValue  - central value of contourVals, optional

# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# contourValsMinCI   - minimum confidence interval bound for contourVals (2d numpy array), optional
# contourValsMaxCI   - maximum confidence interval bound for contourVals (2d numpy array), optional
# Note: contourValsMinCI and contourValsMaxCI must be specified together

# dmin, dmax    - min/max values of contourVals, optional
# interiorLabels- whether to add titles, axis, and colorbar labels for interior subplots, optional

    ax = fig.add_subplot(ny, nx, iplot+1)

    if (np.isnan(contourVals)).all():
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    xVals_ = pu.TDeltas2Seconds(xVals)

    # standardize color limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, contourVals, centralValue, buffr=0.05)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    # default cAxisType
    cAxisType = 'sequential'

    if centralValue is not None and pu.isfloat(centralValue):
        cAxisType = 'diverging'

    # log contours
    isLogScale = logScale
    finite = np.isfinite(contourVals)
    nonzero = finite
    nonzero[finite] = np.greater(np.abs(contourVals[finite]), 0.)
    if nonzero.sum()==0:
        isLogScale = False

    elif isLogScale:
        if centralValue is None:
            minp_, maxp_ = pu.get_clean_ax_limits(
                np.log10(np.nanmax(np.array([dmin, dmax*1.e-3]))),
                np.log10(dmax), np.log10(contourVals[nonzero]), buffr=0.05)
            minp_ = np.power(10., minp_)
            maxp_ = np.power(10., maxp_)
            if maxp_ / minp_ < 20.:
                isLogScale = False
            else:
                minp = minp_
                maxp = maxp_
                cint = contourVals.astype(int)
                isInt = np.all((contourVals - cint) == 0)
                if isInt:
                    minp = np.nanmax(np.array([1., minp]))
                lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)
        elif float(centralValue) == 0.:
            lognorm = mcolors.SymLogNorm(vmin=minp, vmax=maxp,
                        linthresh=1.e-3*maxp, linscale=1.3, base=10)
        else:
            vmin = np.nanmin(np.abs(contourVals[nonzero]))
            vmax = np.nanmax(np.abs(contourVals[nonzero]))

            maxratio = np.max([vmax/centralValue, centralValue/vmin])
            maxp = maxratio / centralValue
            p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
            maxratio = maxp / centralValue
            minp = centralValue / maxratio

            lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)

    # plot contour
    # option 1: smoothed contours
    #cp = ax.contourf(xVals_, yVals, contourVals, pstyle.cmaps[cAxisType]['n'], cmap=pstyle.cmaps[cAxisType]['map'], extend='both',
    #     vmin=minp, vmax=maxp)

    # option 2: pixel contours
    cmap = copy(pstyle.cmaps[cAxisType]['map'])
    #cmap.set_bad(color='gray', alpha=0.5)
    if isLogScale:
        norm = lognorm
    #elif cAxisType == 'diverging':
    #    norm = mcolors.CenteredNorm(centralValue)
    else:
        levels = mticker.MaxNLocator(nbins=pstyle.cmaps[cAxisType]['n']).tick_values(minp, maxp)
        norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # (x, y) specify corners instead of centers
    xVals_edges, yVals_edges = transformXY_to_edges(xVals_,yVals)

    pcm = ax.pcolormesh(xVals_edges, yVals_edges, contourVals, cmap=cmap, norm=norm)

    # Add significance regions if specified
    if (contourValsMinCI is not None and
       np.isfinite(contourValsMinCI).sum() > 0 and
       contourValsMaxCI is not None and
       np.isfinite(contourValsMaxCI).sum() > 0 ):

        # test statistical significance versus centralValue
        if centralValue is not None:
            isSignificant = np.multiply(np.subtract(contourValsMinCI, centralValue),
                                     np.subtract(contourValsMaxCI, centralValue))
            centralValue_ = centralValue
            isSignificant[~np.isfinite(isSignificant)] = -1.0

            # Create patches indicating significance
            w100 = xVals_edges[1:] - xVals_edges[0:-1]
            h100 = yVals_edges[1:] - yVals_edges[0:-1]
            patches = []

            q = 97.73 # 2 sigma, or get quartile from argument in place of CI
            quartile = np.full_like(isSignificant, q)
            cellFraction = np.power(np.abs(quartile - 50.) / 50., 16.) # distort patch with large power
            cellFraction[isSignificant<=0.] = 0.

            for ix in range(len(xVals)):
              for iy in range(len(yVals)):
                if cellFraction[iy, ix] > 0.:

                  # patch center
                  xc = xVals_[ix]
                  yc = yVals[iy]

                  # patch width and height
                  w = w100[ix] * cellFraction[iy, ix]
                  h = h100[iy] * cellFraction[iy, ix]

                  # add Ellipse patch
                  patches.append(
                    mpatches.Ellipse((xc, yc), w, h, fill = False)
                  )

                  ## origin of Rectangle (bottom left)
                  #x0 = xc - w / 2.
                  #y0 = yc - h / 2.

                  ## add Rectangle patch
                  #patches.append(
                  #  mpatches.Rectangle((x0, y0), w, h, fill = False)
                  #)

            # Create patch collection with specified color/alpha
            pc = mcollections.PatchCollection(patches,
                edgecolor = 'black',
                facecolor = 'none',
                snap = False,
                alpha = 0.5,
                linewidth = 0.05,
                linestyles = (0, (5, 1)), # (offset, onoffseq), where onoffseq is an even length tuple of on and off ink lengths in points
            )

            # Add collection to axes
            ax.add_collection(pc)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

    #axes settings
    pu.format_x_for_dates(ax, xVals)
    if xConfig['transform'] is not None:
        if xConfig['transform'] == 'logit':
            ax.set_xscale('logit', nonpositive='clip')
            ax.xaxis.set_major_locator(mticker.AutoLocator())
        else:
            transform = axisTransform(xConfig['transform'])
            locator = transform.locator()
            ax.set_xscale('function', functions=(transform.forward(), transform.inverse()))
            ax.xaxis.set_major_locator(locator(xVals))
            ax.tick_params(axis='x', rotation=90)

    if yConfig['transform'] is not None:
      if yConfig['transform'] == 'logit':
          ax.set_yscale('logit', nonpositive='clip')
          ax.yaxis.set_major_locator(mticker.AutoLocator())
      else:
          transform = axisTransform(yConfig['transform'])
          locator = transform.locator()
          ax.set_yscale('function', functions=(transform.forward(), transform.inverse()))
          ax.yaxis.set_major_locator(locator(yVals))

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if xLabel is not None: ax.set_xlabel(xLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if yLabel is not None: ax.set_ylabel(yLabel,fontsize=axisLabelFontSize)

    #colorbar
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(contourVals)
    m.set_norm(norm)
    if not isLogScale:
        m.set_clim(minp, maxp)
    cb = plt.colorbar(m, ax=ax)
    #scientific formatting
    if sciTicks and not logScale:
        cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    cb.ax.yaxis.get_offset_text().set_fontsize(3)

    if centralValue is not None and float(centralValue) == 1. and isLogScale:
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(cb.ax.yaxis)

    #cb.ax.tick_params(labelsize=3)
    cb.ax.tick_params(axis='y', which='major', labelsize=3)
    cb.ax.tick_params(axis='y', which='minor', labelsize=2)

    if cLabel is not None: cb.set_label(cLabel,fontsize=cbarLabelFontSize)

    if xConfig['invert']:
        ax.invert_xaxis()

    if yConfig['invert']:
        ax.invert_yaxis()

    # optionally add a grid
    #ax.grid(linewidth=gridLineWidth)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return



###############################################################################
def map2D(fig,
          lonVals, latVals, contourVals,
          title='', cLabel='',
          sciTicks=False, logScale=False, centralValue=None,
          ny=1, nx=1, nplots=1, iplot=0,
          contourValsMinCI=None, contourValsMaxCI=None,
          dmin=np.NaN, dmax=np.NaN,
          interiorLabels=True):

# ARGUMENTS
# fig           - matplotlib figure object
# lonVals       - longitude
# latVals       - latitude
# contourVals   - dependent variable (2d numpy array), len(latVals) x len(lonVals)

# title         - subplot title, optional
# cLabel        - label for dependent variable, optional
# sciTicks      - whether contourVals needs scientific formatting for ticks, optional
# logScale      - whether contours are spaced logarithmically, optional, overrides sciTicks
# centralValue  - central value of contourVals, optional
# ny, nx        - number of subplots in x/y direction, optional
# nplots        - total number of subplots, optional
# iplot         - this subplot index (starting at 0), optional

# contourValsMinCI   - minimum confidence interval bound for contourVals (2d numpy array), optional
# contourValsMaxCI   - maximum confidence interval bound for contourVals (2d numpy array), optional

# dmin, dmax    - min/max values of contourVals, optional
# interiorLabels- whether to add titles, axis, and colorbar labels for interior subplots, optional

    # setup map
    cLon = None

    # Is shortest range of lonVals <= 180 degrees?
    # check two different total longitude ranges on periodic domain

    # 0 < longitude <= 360
    lonVals_360 = deepcopy(lonVals)
    while np.max(lonVals_360) >= 360.0:
        lonVals_360[lonVals_360 >= 360.0] -= 360.0
    while np.min(lonVals_360) < 0.0:
        lonVals_360[lonVals_360 < 0.0] += 360.0

    # -180 < longitude <= 180
    lonVals_180 = deepcopy(lonVals_360)
    lonVals_180[lonVals_180 > 180.0] -= 360.0

    for lon in [lonVals_360, lonVals_180]:
        if np.max(lon) - np.min(lon) <= 180.0:
            cLon = 0.5*(np.max(lon) + np.min(lon))

    if cLon is None:
        # plot entire Earth
        ax = fig.add_subplot(ny, nx, iplot+1, projection=ccrs.Mollweide(0.0))

    else:
        # plot single projected side of Earth
        if cLon > 180.0: cLon-=360.0
        ax = fig.add_subplot(ny, nx, iplot+1, projection=ccrs.Orthographic(cLon))


    gl = ax.gridlines(linewidth=0.01, color='gray', alpha=0.5, linestyle='--')

    # all plots
    gl.left_labels = True
    gl.ylocator = LatitudeLocator()
    gl.yformatter = LatitudeFormatter()
    gl.ylabel_style = {'size': 3, 'color': 'black'}
    gl.xlabel_style = {'size': 3, 'color': 'black'}

    # only global projections
    if cLon is None:
        gl.bottom_labels = True
        gl.xlocator = LongitudeLocator()
        gl.xformatter = LongitudeFormatter()

    # standardize color limits
    mindval, maxdval = pu.get_clean_ax_limits(dmin, dmax, contourVals, centralValue, buffr=0.05)
    minp = None
    maxp = None
    if np.isfinite(mindval): minp = mindval
    if np.isfinite(maxdval): maxp = maxdval

    # default cAxisType
    cAxisType = 'sequential'
    if centralValue is not None and pu.isfloat(centralValue):
      cAxisType = 'diverging'

    # log contours
    isLogScale = logScale
    finite = np.isfinite(contourVals)
    nonzero = finite
    nonzero[finite] = np.greater(np.abs(contourVals[finite]), 0.)
    if nonzero.sum()==0:
        isLogScale = False
    elif isLogScale:
        if centralValue is None:
            minp_, maxp_ = pu.get_clean_ax_limits(
                np.log10(np.nanmax(np.array([dmin, dmax*1.e-3]))),
                np.log10(dmax), np.log10(contourVals[nonzero]), buffr=0.05)
            minp_ = np.power(10., minp_)
            maxp_ = np.power(10., maxp_)
            if maxp_ / minp_ < 20.:
                isLogScale = False
            else:
                minp = minp_
                maxp = maxp_
                cint = contourVals.astype(int)
                isInt = np.all((contourVals - cint) == 0)
                if isInt:
                    minp = np.nanmax(np.array([1., minp]))
                lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)
        elif float(centralValue) == 0.:
            lognorm = mcolors.SymLogNorm(vmin=minp, vmax=maxp,
                        linthresh=1.e-3*maxp, linscale=1.3, base=10)
        else:
            vmin = np.nanmin(np.abs(contourVals[nonzero]))
            vmax = np.nanmax(np.abs(contourVals[nonzero]))

            maxratio = np.max([vmax/centralValue, centralValue/vmin])
            maxp = maxratio / centralValue
            p, maxp = pu.get_clean_ax_limits(centralValue, maxp, centralValue=centralValue)
            maxratio = maxp / centralValue
            minp = centralValue / maxratio

            lognorm = mcolors.LogNorm(vmin=minp, vmax=maxp)


    cmap = copy(pstyle.cmaps[cAxisType]['map'])
    cmap.set_bad(color='gray', alpha=0.)
    if isLogScale:
        norm = lognorm
    #elif cAxisType == 'diverging':
    #    norm = mcolors.CenteredNorm(centralValue)
    else:
        levels = mticker.MaxNLocator(nbins=pstyle.cmaps[cAxisType]['n']).tick_values(minp, maxp)
        norm = mcolors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    ## transform to pcolormesh and cartopy conventions
    # longitude monotonically increasing
    lonSort = np.argsort(lonVals_180)
    lonVals_180 = lonVals_180[lonSort]
    contourVals_ = np.empty_like(contourVals)
    for ii, jj in enumerate(lonSort):
      contourVals_[:,ii] = contourVals[:,jj]

    # lat/lon vectors specify corners instead of centers
    lonVals_edges, latVals_edges = transformXY_to_edges(lonVals_180, latVals)

    pcm = ax.pcolormesh(lonVals_edges, latVals_edges, contourVals_, cmap=cmap, norm=norm,
        transform=ccrs.PlateCarree())

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

    # show full projection extent
    ax.set_global()

    # add coastlines
    ax.coastlines(linewidth=coastLineWidth)

    #colorbar
    m = plt.cm.ScalarMappable(cmap=cmap)
    m.set_array(contourVals)
    m.set_norm(norm)
    if not isLogScale:
        m.set_clim(minp, maxp)
    cb = plt.colorbar(m, ax=ax)
    #scientific formatting
    if sciTicks and not logScale:
        cb.ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    cb.ax.yaxis.get_offset_text().set_fontsize(3)

    if centralValue is not None and float(centralValue) == 1. and isLogScale:
        axTicks = OneCenteredAxisTicks(maxp)
        axTicks.setTicks(cb.ax.yaxis)

    #cb.ax.tick_params(labelsize=3)
    cb.ax.tick_params(axis='y', which='major', labelsize=3)
    cb.ax.tick_params(axis='y', which='minor', labelsize=2)
    cb.set_label(cLabel,fontsize=cbarLabelFontSize)

    # label physical distance to the left and up:
    #label = pstyle.subplotLabel(iplot)
    #trans = mtransforms.ScaledTranslation(-20/72, 7/72, fig.dpi_scale_trans)
    #ax.text(0.0, 1.0, label, transform=ax.transAxes + trans,
    #        fontsize='medium', va='bottom', fontfamily='serif')

    return


###############################################################################
def transformXY_to_edges(xs,ys):
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

    return np.array(xs_pcolormesh_midpoints), np.array(ys_pcolormesh_midpoints)


###############################################################################
lenWarnPDF = 0
nanWarnPDF = 0
def plotPDF(fig,
            countsVals, xVals,
            countLabels,
            title="",
            indepLabel=None,
            ny=1, nx=1, nplots=1, iplot=0,
            lineAttribOffset=1,
            legend_inside=True,
            interiorLabels=True,
            normalized=True,
            standardGaussian=False,
):

# ARGUMENTS
# fig              - matplotlib figure object
# countsVals       - list of arrays, each containing counts across xVals
# xVals            - independent variable on x-axis (array)
# countLabels      - legend label for countsVals (list)

# title            - subplot title, optional
# indepLabel       - label for xVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# lineAttribOffset - offset for selecting line attributes, optional
# legend_inside    - whether legend should be placed inside the subplot, optional
# interiorLabels   - whether to add titles and axis labels for interior subplots, optional
# normalized       - whether to normalize countsVals in order to produce a "PDF"
# standardGaussian - whether to overlay a standard Gaussian distribution

    ax = fig.add_subplot(ny, nx, iplot+1)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

    #add counts
    plotVals = []
    nPDFs = 0
    jhist = -1
    for ihist, countVals in enumerate(countsVals):
        if np.all(np.isnan(countVals)) or np.all(countVals == 0):
            global nanWarnPDF
            if nanWarnPDF==0:
                _logger.warning("skipping all-NaN or all-Zero data")
                _logger.warning(title+"; "+indepLabel+"; "+countLabels[ihist])
            nanWarnPDF=nanWarnPDF+1
            continue
        if len(countVals)!=len(xVals):
            global lenWarnPDF
            if lenWarnPDF==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+indepLabel+"; "+countLabels[ihist])
            lenWarnPDF=lenWarnPDF+1
            continue

        jhist += 1

        # Plot line for each countVals that has non-missing data
        # only plot non-zero counts
        yVals = np.array(countVals, dtype=np.float64)
        positive = np.greater(yVals, 0.)
        xVals_ = xVals[positive]
        yVals = yVals[positive]

        if normalized:
          # assume constant dx between bins
          dx = xVals[1] - xVals[0]
          yVals = np.divide(yVals,np.sum(countVals[positive])*dx)

        ax.plot(xVals_, yVals,
                color=pstyle.color(jhist+lineAttribOffset),
                label=countLabels[ihist],
                ls=pstyle.line(jhist+lineAttribOffset),
                linewidth=commonLineWidth)
        nPDFs = nPDFs + 1
        plotVals.append(yVals)

    if nPDFs == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return

    # add a standard Gaussian pdf
    if standardGaussian:
        from scipy.stats import norm
        ax.plot(xVals, norm.pdf(xVals),
                    color='0.4',
                    ls='-.',
                    linewidth=0.35,
                    label='N(0,1)'
        )
        ax.set_ylim(bottom=5.e-4)

    #axes settings
    ax.set_yscale('log')

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)

    yLabel = 'Count'
    if normalized: yLabel = 'PDF'
    if interiorLabels or iy == ny-1:
        if indepLabel is not None: ax.set_xlabel(indepLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        ax.set_ylabel(yLabel,fontsize=axisLabelFontSize)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=legendLabelFontSize2,frameon=True,
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=legendLabelFontSize2,frameon=False,
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid(linewidth=gridLineWidth)

    return


###############################################################################
lenWarnRamp = 0
nanWarnRamp = 0
def plotCompositeSeries(fig,
               xVals,
               countVals,
               meanVals,
               rmsVals,
               stdVals,
               title="", dataLabel=None,
               indepLabel=None,
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

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

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
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(xVals, lineVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth)
        nLines += 1
        plotVals.append(lineVals)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return None

    #axes settings
    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(plotVals=plotVals)
    if (np.isfinite(mindval) and
        np.isfinite(maxdval)):
        ax.set_ylim(mindval, maxdval)

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if indepLabel is not None: ax.set_xlabel(indepLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if dataLabel is not None: ax.set_ylabel(dataLabel,fontsize=axisLabelFontSize)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid(linewidth=gridLineWidth)

    # Add count on RHS y-axis
    ax2 = ax.twinx()
    color = 'black'
    if interiorLabels or ix == nx-1:
        ax2.set_ylabel('Count',fontsize=axisLabelFontSize,color=color)
    ax2.plot(xVals, countVals,
             color=color,
             label='Count',
             ls=':',
             linewidth=commonLineWidth)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_yscale('log')
    ax2.set_ylim(bottom=100.)
    ax2.tick_params(axis='y', which='major', labelsize=3)
    ax2.tick_params(axis='y', which='minor', labelsize=2)

    return None



###############################################################################
lenWarnRamp = 0
nanWarnRamp = 0
def plotCompositeProfile(fig,
               yVals,
               countVals,
               meanVals,
               rmsVals,
               stdVals,
               title="", dataLabel=None,
               indepLabel=None,
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

# yVals           - independent variable on y-axis (array)

# title            - subplot title, optional
# dataLabel        - label for x-axis, optional
# indepLabel       - label for yVals, optional

# ny, nx           - number of subplots in x/y direction, optional
# nplots           - total number of subplots, optional
# iplot            - this subplot index (starting at 0), optional

# lineAttribOffset - offset for selecting line attributes, optional
# legend_inside    - whether legend should be placed inside the subplot, optional

    ax = fig.add_subplot(ny, nx, iplot+1)
    ix = int(iplot)%int(nx)
    iy = int(iplot)/int(nx)

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

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
        if len(lineVals)!=len(yVals):
            global lenWarnRamp
            if lenWarnRamp==0:
                _logger.warning("skipping data where len(x)!=len(y)")
                _logger.warning(title+"; "+indepLabel+"; "+linesLabel[iline])
            lenWarnRamp=lenWarnRamp+1
            continue

        # Plot line for each lineVals that has non-missing data
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(lineVals, yVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth)
        nLines += 1
        plotVals.append(lineVals)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return None

    #axes settings
    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(plotVals=plotVals)
    if (np.isfinite(mindval) and
        np.isfinite(maxdval)):
        ax.set_xlim(mindval, maxdval)

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if dataLabel is not None: ax.set_xlabel(dataLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if indepLabel is not None: ax.set_ylabel(indepLabel,fontsize=axisLabelFontSize)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid(linewidth=gridLineWidth)

    # Add count on RHS y-axis
    ax2 = ax.twiny()
    color = 'black'
    if interiorLabels or iy == 0:
        ax2.set_xlabel('Count',fontsize=axisLabelFontSize,color=color)
    ax2.plot(countVals, yVals,
             color=color,
             label='Count',
             ls=':',
             linewidth=commonLineWidth)
    ax2.tick_params(axis='x', labelcolor=color)
    ax2.set_xscale('log')
    ax2.set_xlim(left=100.)
    ax2.tick_params(axis='x', which='major', labelsize=3)
    ax2.tick_params(axis='x', which='minor', labelsize=2)

    return None

###############################################################################
lenWarnRamp = 0
nanWarnRamp = 0
def plotfitRampComposite(fig,
               xVals,
               countVals,
               meanVals,
               rmsVals,
               stdVals,
               title="", dataLabel=None,
               indepLabel=None,
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

    # title
    t = title
    if nplots > 1:
      t = pstyle.subplotLabel(iplot)+' '+title
    ax.set_title(t,fontsize=titleFontSize)

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
        pColor = pstyle.color(iline+lineAttribOffset)

        ax.plot(xVals, lineVals,
                color=pColor,
                label=linesLabel[iline],
                ls=pstyle.line(iline+lineAttribOffset),
                linewidth=commonLineWidth)
        nLines += 1
        plotVals.append(lineVals)

    if nLines == 0:
        ax.tick_params(axis='x',labelbottom=False)
        ax.tick_params(axis='y',labelleft=False)
        return None

    fitQuantity = deepcopy(stdVals)
    fitStat = 'STD'

    #fitQuantity = deepcopy(rmsVals)
    #fitStat = 'RMS'

    # Add fit for fitQuantity here using info from countVals

    # determine maximum x index for finite values
    indMaxX4Quant = 0
    for ii, quant in enumerate(fitQuantity):
        if np.isfinite(quant):
            indMaxX4Quant = ii
        elif indMaxX4Quant==0:
            continue
        else:
            break

#    indMaxCount = np.argmax(countVals)

#    # determine maximum x index for linear fit
#    indMaxXFit = indMaxX4Quant
#    maxCount = np.nanmax(countVals+[0])
#    for ii, count in enumerate(countVals):
#        if float(count) >= 0.002*float(maxCount) and ii >= indMaxCount and ii <= indMaxX4Quant:
#            indMaxXFit = ii
#
#    # determine minimum x index for finite values
#    indMinX4Quant = 0
#    for jj, quant in enumerate(np.flip(fitQuantity)):
#        ii = len(fitQuantity)-jj-1
#        if np.isfinite(quant) and ii < indMaxX4Quant:
#            indMinX4Quant = ii
#        elif indMinX4Quant == 0:
#            continue
#        else:
#            break
#
#    # determine minimum x index for linear fit
#    indMinXFit = indMinX4Quant
#    for jj, count in enumerate(np.flip(countVals)):
#        ii = len(countVals)-jj-1
#        if float(count) >= 0.002*float(maxCount) and ii <= indMaxCount and ii < indMaxXFit and ii >= indMinX4Quant:
#            indMinXFit = ii

    validQuant = np.isfinite(fitQuantity)
    validCount = np.greater(countVals,0)
    valid = np.logical_and(validQuant, validCount)

    qFit = np.array(fitQuantity)[valid]
    cFit = np.array(countVals)[valid]
    xFit = np.array(xVals)[valid]

    indMaxCount = np.argmax(cFit)
    maxCount = np.nanmax(cFit+[0])

    # determine maximum x index for linear fit
    indMaxXFit = len(qFit)-1
    for ii, count in enumerate(cFit):
        if float(count) >= 0.002*float(maxCount) and ii >= indMaxCount:
            indMaxXFit = ii

    # determine minimum x index for linear fit
    indMinXFit = 0
    for jj, count in enumerate(np.flip(cFit)):
        ii = len(countVals)-jj-1
        if float(count) >= 0.002*float(maxCount) and ii <= indMaxCount:
            indMinXFit = ii

    # should not need finite test after validQuant test above
    #finite = np.isfinite(qFit[indMinXFit:indMaxXFit+1])
    #inds = np.arange(indMinXFit, indMaxXFit+1)
    #indMaxQ = np.argmax(qFit[indMinXFit:indMaxXFit+1][finite])
    #indMaxQ = inds[finite][indMaxQ]

    inds = np.arange(indMinXFit, indMaxXFit+1)
    indMaxQ = np.argmax(qFit[indMinXFit:indMaxXFit+1])
    indMaxQ = inds[indMaxQ]

    ERRfitDict = None

    # only create fit if there is more than 1 point
    if indMaxQ-indMinXFit+1 > 1:
      ERRfitDict = {}

      weights = cFit[indMinXFit:indMaxQ+1]
      weights = weights / weights.sum()

      #print(title+"; "+indepLabel)
      #print(xFit[indMinXFit:indMaxQ+1])
      #print(qFit[indMinXFit:indMaxQ+1])
      #print(weights)

      p = np.polyfit(xFit[indMinXFit:indMaxQ+1],qFit[indMinXFit:indMaxQ+1],1,
                     w=weights)

      X0 = xFit[indMinXFit]
      ERR0 = X0 * p[0] + p[1]

      # X1 = xFit[indMaxQ]
      # ERR1 = X1 * p[0] + p[1]
      ERR1 = qFit[indMaxQ]
      X1 = (ERR1 - p[1]) / p[0]

      ERRfitDict['bin_utils'] = {
          'X':  [round(X0,2),  round(X1,2)],
          'ERR': [round(ERR0,2), round(ERR1,2)],
      }
      ERRfitDict['YAML'] = {
          'x0': [round(X0,2)],
          'x1': [round(X1,2)],
          'err0': [round(ERR0,2)],
          'err1': [round(ERR1,2)],
      }

      fitX  = np.asarray([0.0]  + ERRfitDict['bin_utils']['X']  + [xVals[indMaxX4Quant]])
      fitERR = np.asarray([ERR0] + ERRfitDict['bin_utils']['ERR'] + [ERR1])

      plotVals.append(fitERR)

      pColor = pstyle.color(linesLabel.index(fitStat)+lineAttribOffset)

      ax.plot(fitX, fitERR,
              color=pColor,
              label='Fit-'+fitStat,
              ls='-.',
              linewidth=commonLineWidth*1.5,
              marker='+',
              ms=1.5
      )

    #axes settings
    # standardize x-limits
    mindval, maxdval = pu.get_clean_ax_limits(plotVals=plotVals)
    if (np.isfinite(mindval) and
        np.isfinite(maxdval)):
        ax.set_ylim(mindval, maxdval)

    ax.tick_params(axis='both', which='major', labelsize=3)
    ax.tick_params(axis='both', which='minor', labelsize=2)

    #handle interior subplot ticks/labels
    if not interiorLabels \
       and (iy < ny-2 or ( iy == ny-2 and (int(nplots)%int(nx)==0 or ix <= (int(nplots)%int(nx) - 1)) )):
        ax.tick_params(axis='x',labelbottom=False)
    if interiorLabels or iy == ny-1:
        if indepLabel is not None: ax.set_xlabel(indepLabel,fontsize=axisLabelFontSize)
    if interiorLabels or ix == 0:
        if dataLabel is not None: ax.set_ylabel(dataLabel,fontsize=axisLabelFontSize)

    #legend
    if legend_inside:
        #INSIDE AXES
        lh = ax.legend(loc='best',fontsize=legendLabelFontSize1,frameon=True,
                       framealpha=0.4,ncol=1)
        lh.get_frame().set_linewidth(0.0)
    elif ix==nx-1 or iplot==nplots-1:
        #OUTSIDE AXES
        ax.legend(loc='upper left',fontsize=legendLabelFontSize1,frameon=False,
                  bbox_to_anchor=(1.02, 1), borderaxespad=0)

    ax.grid(linewidth=gridLineWidth)

    # Add count on RHS y-axis
    ax2 = ax.twinx()
    color = 'black'
    if interiorLabels or ix == nx-1:
        ax2.set_ylabel('Count',fontsize=axisLabelFontSize,color=color)
    #ax2.plot(xVals[:indMaxX4Quant], countVals[:indMaxX4Quant],
    ax2.plot(xVals, countVals,
             color=color,
             label='Count',
             ls=':',
             linewidth=commonLineWidth)
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.set_yscale('log')
    ax2.set_ylim(bottom=100.)
    ax2.tick_params(axis='y', which='major', labelsize=3)
    ax2.tick_params(axis='y', which='minor', labelsize=2)

    return ERRfitDict
