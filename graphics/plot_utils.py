#!/usr/bin/env python3

from collections.abc import Iterable
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import ConciseDateFormatter, DateFormatter, AutoDateLocator
import numpy as np
import os
import pandas as pd


#============================
# figure/plotting definitions
#============================

colorIterator = [
#  [0.6350, 0.0780, 0.1840],
  [0., 0., 0.],
  [0.0000, 0.4470, 0.7410],
  [0.8500, 0.3250, 0.0980],
  [0.9290, 0.6940, 0.1250],
  [0.4940, 0.1840, 0.5560],
  [0.4660, 0.6740, 0.1880],
  [0.3010, 0.7450, 0.9330],
  [0.6350, 0.0780, 0.1840],
]
greyIterator = [
  [0.2, 0.2, 0.2],
  [0.4, 0.4, 0.4],
  [0.6, 0.6, 0.6],
  [0.8, 0.8, 0.8],
]

styleIterator = ['-', '--', '-.', ':']

# fast varying line colors
# produces: 'k-', 'b-', 'g-', ..., 'k--', 'b--', 'g--', ..., 'm:'
defaultPColors = colorIterator*len(styleIterator)
defaultPLineStyles = []
for style in styleIterator:
  defaultPLineStyles += [style]*len(colorIterator)

# fast varying line styles
# produces: 'k-', 'k--', 'k-.', 'k:', ..., 'm-', 'm--', 'm-.', 'm:'
#defaultPLineStyles = styleIterator*len(colorIterator)
#defaultPColors = []
#for color in colorIterator:
#  defaultPColors += [color]*len(styleIterator)

def plotColor(nLines = 1, index = 0, nSpaghetti = None):
  if nSpaghetti is not None and nLines >= nSpaghetti:
    pColors = ['0.45']*nSpaghetti
    for i in list(range(0, nLines - nSpaghetti + 1)):
      pColors += [plotColor(1, i)]
    return pColors[index]
  else:
    return defaultPColors[np.mod(index,len(defaultPColors))]

def plotLineStyle(nLines = 1, index = 0, nSpaghetti = None):
  if nSpaghetti is not None and nLines >= nSpaghetti:
    pLineStyles = ['--']*nSpaghetti
    for i in list(range(0, nLines - nSpaghetti + 1)):
      pLineStyles += [plotLineStyle(1, i)]
    return pLineStyles[index]
  else:
    return defaultPLineStyles[np.mod(index,len(defaultPLineStyles))]

plotSpecs = ['k-*', 'b-*', 'g-*', 'r-*', 'c-*', 'm-*',
             'k--+','b--+','g--+','r--+','c--+','m--+']
plotMarkers = ['*','*','*','*','*','*',
               '+','+','+','+','+','+']


###############################################################################
def setup_fig(nx=1, ny=1, inch_width=1.5, aspect=1.0, ybuffer=True):
#INPUTS
# nx - number of subplots in x direction
# ny - number of subplots in y direction
# inch_width - rough subplot size in inches
# ybuffer - whether to give extra y space for labeling
#
#OUTPUT
# fig - a new figure with standard sizing

    fig = plt.figure()

    if ybuffer:
        fig.set_size_inches(nx*inch_width,aspect*ny*inch_width)
    else:
        fig.set_size_inches(0.9*nx*inch_width,0.9*aspect*ny*inch_width)

    return(fig)


###############################################################################
def finalize_fig(fig, filename='temporary_figure', filetype='png',
                 ybuffer=True, xbuffer=False):
#INPUTS
# fig - plt.figure() type
# filename - name of figure file without extension
# filetype - file extension, e.g., 'png'
# x/ybuffer - whether to give extra x/y space for labeling

    wspace = 0.35
    if xbuffer: wspace = 0.6

    hspace = 0.40
    if ybuffer: hspace = 0.70

    fig.subplots_adjust(wspace=wspace,hspace=hspace)

    if filetype == 'png':
        fig.savefig(filename+'.'+filetype,dpi=200,bbox_inches='tight')
    if filetype == 'pdf':
        fig.savefig(filename+'.'+filetype,bbox_inches='tight')

    plt.close(fig)


###############################################################################
def TDeltas2Seconds(x_):
    if isinstance(x_[0],dt.timedelta):
        x = []
        for xVal in x_:
            x.append(xVal.total_seconds())
        return x
    return x_


###############################################################################
def timeDeltaTicks(x, pos):
    d = dt.timedelta(seconds=x)
    i = '{:d}'
    i02 = '{:02d}'
    vals = {}
    fmts = {}
    prefs = {}
    suffs = {}
    depends = {}

    vals['D'] = d.days
    fmts['D'] = i
    prefs['D'] = ''
    suffs['D'] = 'd '
    depends['D'] = ['D']

    vals['HH'], hrem = divmod(d.seconds, 3600)
    fmts['HH'] = i02
    prefs['HH'] = ''
    suffs['HH'] = ''
    #depends['HH'] = ['HH','MM','SS']

    vals['MM'], vals['SS'] = divmod(hrem, 60)
    fmts['MM'] = i02
    prefs['MM'] = ':'
    suffs['MM'] = ''
    depends['MM'] = ['MM','SS']
    fmts['SS'] = i02
    prefs['SS'] = ':'
    suffs['SS'] = ''
    depends['SS'] = ['SS']

    if vals['MM'] == 0 and vals['SS'] == 0:
        fmts['HH'] = i
        suffs['HH'] = 'h'

    out = ''
    for key in vals.keys():
        include = False
        if key in depends:
            for dep in depends[key]:
                if vals[dep] > 0: include = True
        else:
            include = True
        if include:
            out += prefs[key]+fmts[key].format(vals[key])+suffs[key]

    return out

#DTimeLocator = AutoDateLocator(interval_multiples=True)
DTimeLocator = AutoDateLocator()
DTimeFormatter = ConciseDateFormatter(DTimeLocator) #DateFormatter('%m-%d_%HZ')
TDeltaFormatter = matplotlib.ticker.FuncFormatter(timeDeltaTicks)


###############################################################################
def format_x_for_dates(ax,x):
    if isinstance(x[0],dt.datetime):
        ax.xaxis.set_major_locator(DTimeLocator)
        ax.xaxis.set_major_formatter(DTimeFormatter)
#        ax.xaxis.set_tick_params(rotation=30)
#        ax.set_xlabel('Date',fontsize=4)
        ax.xaxis.get_offset_text().set_fontsize(3)
    if isinstance(x[0],dt.timedelta):
        x = TDeltas2Seconds(x)
        ax.set_xlim(min(x),max(x))
        tstep = 3600*3 #3 hours
        ntick = 500
        while ntick > 8:
            ticks = np.arange(x[0],x[-1]+tstep,tstep)
            ntick = len(ticks)
            tstep = tstep * 2
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(TDeltaFormatter)
        ax.xaxis.set_tick_params(rotation=30)
        ax.set_xlabel('Lead Time',fontsize=4)


###############################################################################
def get_clean_ax_limits(xmin_=np.NaN, xmax_=np.NaN, plotVals=[np.NaN],
                        centralValue=None, buffr=0.2):

    x = np.empty(2)
    for ii, (func, x_) in enumerate(zip(
      [np.nanmin, np.nanmax],
      [xmin_, xmax_],
    )):
      if np.isfinite(x_):
        x[ii] = x_
      else:
        if isinstance(plotVals[0], Iterable):
          x[ii] = func(np.asarray([func(p) for p in plotVals]))
        else:
          x[ii] = func(plotVals)
    xmin = x[0]
    xmax = x[1]
    #print(xmin, centralValue, xmax)

    xmaxabs=np.nanmax([abs(xmin), abs(xmax)])*(1.+buffr)
    if not np.isfinite(xmaxabs) or xmaxabs == 0.0:
        minxval = 0.0
        maxxval = 1.0
    else:
        roundfact = np.round(1. / 10.0 ** np.floor(np.log10(xmaxabs)))
        if not np.isfinite(roundfact) or roundfact <= 0.0: roundfact = 1.0

        if centralValue is None:
            maxxval = np.ceil(    xmax*roundfact ) / roundfact
            minxval = np.floor(   xmin*roundfact ) / roundfact
        elif xmin < 0. and xmax > 0. and centralValue == 0.0:
            maxxval = np.ceil(    (xmaxabs)*roundfact ) / roundfact
            minxval = np.floor( - (xmaxabs)*roundfact ) / roundfact
        else:
            if xmax > centralValue:
                roundfact = np.round(1. / 10.0 ** np.floor(np.log10(xmax - centralValue)))
                if not np.isfinite(roundfact) or roundfact <= 0.0: roundfact = 1.0
                delta = np.ceil(  (xmax - centralValue)*roundfact ) / roundfact
                maxxval = centralValue + delta
                minxval = centralValue - delta
            elif xmin < centralValue:
                roundfact = np.round(1. / 10.0 ** np.floor(np.log10(centralValue - xmin)))
                if not np.isfinite(roundfact) or roundfact <= 0.0: roundfact = 1.0
                delta = np.ceil(  (centralValue - xmin)*roundfact ) / roundfact
                maxxval = centralValue + delta
                minxval = centralValue - delta
            else:
                roundfact = np.round(1. / 10.0 ** np.floor(np.log10(xmax)))
                maxxval = np.ceil(  xmax*roundfact ) / roundfact
                minxval = np.floor( xmin*roundfact ) / roundfact



    return minxval, maxxval


#===========================
# general purpose functions
#===========================

def uniqueMembers(listVar):
#PURPOSE return a list without repeating members
    output = []
    for x in listVar:
        if x not in output:
            output.append(x)
    return output


#########################################################
def isfloat(value):
#PURPOSE determine if value can be converted to float
  try:
    float(value)
    return True
  except ValueError:
    return False


#########################################################
def isint(value):
#PURPOSE determine if value can be converted to int
  try:
    int(value)
    if isfloat(value):
        if float(value).is_integer():
            return True
        else:
            return False
    else:
        return True
  except ValueError:
    return False

#########################################################
def prepends(str1: str, str2: str):
  return str1 == str2[0:min([len(str2), len(str1)])]

#########################################################
def postpends(str1: str, str2: str):
  return str1 == str2[-min([len(str2), len(str1)]):]
