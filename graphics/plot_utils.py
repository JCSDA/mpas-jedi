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

colorIterator = ['k','b','g','r','c','m']

styleIterator = ['-','--','-.',':']

# fast varying line colors
defaultPColors = colorIterator*len(styleIterator)
defaultPLineStyles = []
for style in styleIterator:
  defaultPLineStyles += [style]*len(colorIterator)

# fast varying line styles
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
def get_clean_ax_limits(xmin_=np.NaN,xmax_=np.NaN,plotVals=[np.NaN],
                        signdef=False,symmetric=True):
    if not np.isnan(xmin_) and not np.isnan(xmax_):
        xmin = xmin_
        xmax = xmax_
    else:
        if isinstance(plotVals[0], Iterable):
            xmin = np.nanmin(plotVals[0])
            xmax = np.nanmax(plotVals[0])
            for vals in plotVals[1:]:
                xmin = np.nanmin([np.nanmin(vals),xmin])
                xmax = np.nanmax([np.nanmin(vals),xmax])
        else:
            xmin = np.nanmin(plotVals)
            xmax = np.nanmax(plotVals)

    xmaxabs=np.nanmax([abs(xmin),abs(xmax)])*1.2
    if xmaxabs == 0.0 or np.isnan(xmaxabs):
        minxval = 0.0
        maxxval = 1.0
    else:
        roundfact = np.round(1. / 10.0 ** np.floor(np.log10(xmaxabs)))
        if np.isnan(roundfact) or roundfact <= 0.0: roundfact = 1.0

        if signdef or not symmetric:
            maxxval = np.ceil(    xmax*roundfact ) / roundfact
            minxval = np.floor(   xmin*roundfact ) / roundfact
        else:
            maxxval = np.ceil(    xmaxabs*roundfact ) / roundfact
            minxval = np.floor( - xmaxabs*roundfact ) / roundfact
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
