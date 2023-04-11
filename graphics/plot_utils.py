#!/usr/bin/env python3

from collections.abc import Iterable
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import ConciseDateFormatter, DateFormatter, AutoDateLocator
import numpy as np
import pandas as pd

#===============
# plot utilities
#===============

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
                 ybuffer=True, xbuffer=0.35, hspace_=0.70):
#INPUTS
# fig - plt.figure() type
# filename - name of figure file without extension
# filetype - file extension, e.g., 'png'
# ybuffer - whether to give extra y space for labeling
# xbuffer - percentage of extra x space for labeling

    wspace = xbuffer

    hspace = 0.40
    if ybuffer: hspace = hspace_

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
def timeDeltaTickLabels(x, pos):
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
    depends['HH'] = ['HH','MM','SS']

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

    # default to 0 days when empty
    if out == '': out = '0d'

    return out

#DTimeLocator = AutoDateLocator(interval_multiples=True)
DTimeLocator = AutoDateLocator()
DTimeFormatter = ConciseDateFormatter(DTimeLocator) #DateFormatter('%m-%d_%HZ')
TDeltaFormatter = mpl.ticker.FuncFormatter(timeDeltaTickLabels)


###############################################################################
def timeDeltaTicks(x):
    x = TDeltas2Seconds(x)
    tstep = 3600*3 #3 hours
    ntick = 500
    while ntick > 8:
        ticks = np.arange(x[0],x[-1]+tstep,tstep)
        ntick = len(ticks)
        tstep = tstep * 2
    return ticks


###############################################################################
def format_x_for_dates(ax, x):
    if isinstance(x[0],dt.datetime):
        ax.xaxis.set_major_locator(DTimeLocator)
        ax.xaxis.set_major_formatter(DTimeFormatter)
        ax.xaxis.get_offset_text().set_fontsize(3)
    if isinstance(x[0],dt.timedelta):
        ticks = timeDeltaTicks(x)
        ax.set_xlim(min(ticks), max(ticks))
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(TDeltaFormatter)
        ax.xaxis.set_tick_params(rotation=30)


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
        elif np.abs(centralValue) < 1.e-6:
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
            elif xmin==xmax and xmin==centralValue:
                maxxval = centralValue*3./2.
                minxval = centralValue*2./3.
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
