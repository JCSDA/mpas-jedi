import datetime as dt
import getopt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, AutoDateLocator
import numpy as np
import os
import pandas as pd
import sys

miss_f = -88888.8
miss_i = -88888
miss_s = 'null'
csvSEP = ";"

#=======================
# diagnostic definitions
#=======================

allDiags = ['omb','oma','obs','bak','ana']
nonObsDiags = [diag for diag in allDiags if diag!='obs']

#=====================
# variable definitions
#=====================

# columns: var_name            unit_used   abbr.
varDict = { \
    'air_temperature':        [ 'K',     'T'   ] \
  , 'air_pressure':           [ 'hPa',   'P'   ] \
  , 'altitude':               [ 'm',     'alt' ] \
  , 'bending_angle':          [ '%',     'Bnd' ] \
  , 'brightness_temperature': [ 'K',     'BT'  ] \
  , 'eastward_wind':          [ 'm/s',   'U'   ] \
  , 'northward_wind':         [ 'm/s',   'V'   ] \
  , 'refractivity':           [ '%',     'Ref' ] \
  , 'specific_humidity':      [ 'kg/kg', 'qv'  ] \
  , 'virtual_temperature':    [ 'K',     'Tv'  ] \
  , 'surface_pressure':       [ 'Pa',    'Ps'  ] \
  , 'latitude':               [ 'deg',   'lat' ] \
  , 'longitude':              [ 'deg',   'lon' ] \
    }

#Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
#Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == 'altitude':

#===================
# binning functions
#===================

def equalBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.equal(x,bound))
    return mask

def notEqualBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.not_equal(x,bound))
    return mask

#def lessEqualBound(xVecs, bound):
#    mask = np.full((len(xVecs[0])),False)
#    for x in xVecs:
#        mask = np.logical_or(mask, np.less_equal(x,bound))
#    return mask

def lessBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.less(x,bound))
    return mask

def greatEqualBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.greater_equal(x,bound))
    return mask

def greatBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.greater(x,bound))
    return mask

#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(xVecs,REGION_NAME): 
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = xVecs[0]
#    lats = xVecs[1]
#    nlocs = len(lons)
##    mask = np.full((nlocs),False)
##    for ii in list(range(nlocs)):
##        mask[ii] = isoutside(lons[ii],lats[ii],shp)
#    mask = isoutside(lons,lats,shp)
#    return mask

#===================
# standard bin info
#===================

# altitude, surface to top
altitudeBinBounds = list(np.arange(1000., 50000.0,  2000.))
altitudeMinBounds = []
altitudeMaxBounds = []
altitudeBinStrVals = []
for ibin in list(range(len(altitudeBinBounds)-1)):
    altitudeMinBounds.append(altitudeBinBounds[ibin])
    altitudeMaxBounds.append(altitudeBinBounds[ibin+1])
    binVal = 0.5 * (altitudeBinBounds[ibin+1] + altitudeBinBounds[ibin])
    altitudeBinStrVals.append("{:.0f}".format(binVal))

# pressure, surface to top
pressureBinBounds = [1050., 950., 850., 750., 650., 550., 450., 350., 250., 150., 50., 0.]
pressureMinBounds = []
pressureMaxBounds = []
pressureBinStrVals = []
for ibin in list(range(len(pressureBinBounds)-1)):
    pressureMinBounds.append(pressureBinBounds[ibin+1])
    pressureMaxBounds.append(pressureBinBounds[ibin])
    binVal = 0.5 * (pressureBinBounds[ibin+1] + pressureBinBounds[ibin])
    pressureBinStrVals.append("{:.0f}".format(binVal))

# latitude, north to south
# 'NAMED'
allNamedLatBands = ['NPol','NXTro','NTro','STro','SXTro','SPol']
allNamedLatBandsMinBounds = [60.0, 30.0,  0.0, -30.0, -60.0, -90.0]
allNamedLatBandsMaxBounds = [90.0, 60.0, 30.0,   0.0, -30.0, -60.0]

namedLatBandsStrVals = ['NXTro','NTro','STro','SXTro']
namedLatBandsMinBounds = []
namedLatBandsMaxBounds = []
for latBand in namedLatBandsStrVals:
    iband = allNamedLatBands.index(latBand)
    namedLatBandsMinBounds.append(allNamedLatBandsMinBounds[iband])
    namedLatBandsMaxBounds.append(allNamedLatBandsMaxBounds[iband])

# 'NUMERICAL'; TODO: add these quantities to 2D figures in plot_stats_timeseries.py
dLat = 10.0
latitudeBinBounds = list(np.arange(90.-0.5*dLat, -90.-0.5*dLat,  -dLat))
latitudeMinBounds = []
latitudeMaxBounds = []
latitudeBinStrVals = []
for ibin in list(range(len(latitudeBinBounds)-1)):
    latitudeMinBounds.append(latitudeBinBounds[ibin+1])
    latitudeMaxBounds.append(latitudeBinBounds[ibin])
    binVal = 0.5 * (latitudeBinBounds[ibin+1] + latitudeBinBounds[ibin])
    latitudeBinStrVals.append("{:.0f}".format(binVal))

# general bin definitions
nullBinDesc = { 'filters': [], 'labels': [miss_s] }

#@EffectiveQC* values:
# pass    = 0;   // we like that one!
# missing = 1;   // missing values prevent use of observation
# preQC   = 2;   // observation rejected by pre-processing
# bounds  = 3;   // observation value out of bounds
# domain  = 4;   // observation not within domain of use
# black   = 5;   // observation black listed
# Hfailed = 6;   // H(x) computation failed
# thinned = 7;   // observation removed due to thinning
# diffref = 8;   // metadata too far from reference
# clw     = 9;   // observation removed due to cloud field
# fguess  = 10;  // observation too far from guess
# seaice  = 11;  // observation based sea ice detection, also flags land points
#
# Static list above copied on 3 Oct 2019
# see ufo/src/ufo/filters/QCflags.h for up-to-date list

goodFlags = [0]
goodFlagNames = ['pass']

badFlags     = [1,         2,         3,         4, \
                5,         6,         7,         8, \
                9,         10,        11]
badFlagNames = ['missing', 'preQC',   'bounds',  'domain', \
                'black',   'Hfailed', 'thinned', 'diffref', \
                'clw',     'fguess',  'seaice']

#========================================
# bin dictionary used for all bin groups
#========================================

vNameStr = 'varName'
varQC = vNameStr+'@EffectiveQC'
prsMeta = 'air_pressure@MetaData'
altMeta = 'altitude@MetaData'
lonMeta = 'longitude@MetaData'
latMeta = 'latitude@MetaData'

goodBinNames = goodFlagNames

binGrpDict = { \
#    #name            binVars ... followed by dictionary defitions of each binDesc
    'qc':      { 'variables': [varQC] \
                 , 'good':      { 'filters': [\
                                    {'where': notEqualBound, 'args': [varQC], 'bounds': goodFlags, \
                                     'apply_to': nonObsDiags}], \
                                  'labels': goodBinNames } \
                 , 'bad':       { 'filters': [\
                                    {'where': notEqualBound, 'args': [varQC], 'bounds': badFlags, \
                                     'apply_to': nonObsDiags} \
                                  , {'where': equalBound,    'args': [varQC], 'bounds': badFlags, \
                                     'mask_value': 0.0, 'apply_to': nonObsDiags}] \
                              , 'labels': badFlagNames } \
                 } \
  , 'pressure':  { 'variables': [varQC,prsMeta] \
                 , 'default':   { 'filters': [\
                                    {'where': lessBound,       'args': [prsMeta], 'bounds': pressureMinBounds} \
                                  , {'where': greatEqualBound, 'args': [prsMeta], 'bounds': pressureMaxBounds} \
                                  , {'where': notEqualBound,   'args': [varQC],   'bounds': goodFlags*len(pressureBinStrVals), \
                                     'apply_to': nonObsDiags} ] \
                                , 'labels': pressureBinStrVals } \
                 } \
  , 'altitude':  { 'variables': [varQC,altMeta] \
                 , 'default':   { 'filters': [\
                                    {'where': lessBound,       'args': [altMeta], 'bounds': altitudeMinBounds} \
                                  , {'where': greatEqualBound, 'args': [altMeta], 'bounds': altitudeMaxBounds} \
                                  , {'where': notEqualBound,   'args': [varQC],   'bounds': goodFlags*len(altitudeBinStrVals), \
                                     'apply_to': nonObsDiags} ] \
                                , 'labels': altitudeBinStrVals } \
                 } \
  , 'latband':   { 'variables': [varQC,latMeta] \
                 , 'NAMED':     { 'filters': [\
                                    {'where': lessBound,       'args': [latMeta], 'bounds': namedLatBandsMinBounds} \
                                  , {'where': greatEqualBound, 'args': [latMeta], 'bounds': namedLatBandsMaxBounds} \
                                  , {'where': notEqualBound,   'args': [varQC],   'bounds': goodFlags*len(namedLatBandsStrVals), \
                                     'apply_to': nonObsDiags} ] \
                                , 'labels': namedLatBandsStrVals } \
                 , 'NUMERICAL': { 'filters': [\
                                    {'where': lessBound,       'args': [latMeta], 'bounds': latitudeMinBounds} \
                                  , {'where': greatEqualBound, 'args': [latMeta], 'bounds': latitudeMaxBounds} \
                                  , {'where': notEqualBound,   'args': [varQC],   'bounds': goodFlags*len(latitudeBinStrVals), \
                                     'apply_to': nonObsDiags} ] \
                                , 'labels': latitudeBinStrVals } \
                 } \
  , 'box':       { 'variables': [varQC,lonMeta,latMeta] \
                 , 'AFRICA':    nullBinDesc \
                 , 'ATLANTIC':  nullBinDesc \
                 , 'AUSTRALIA': nullBinDesc \
                 , 'CONUS':     { 'filters': [\
                                    {'where': lessBound,     'args': [lonMeta], 'bounds': [234.0]} \
                                  , {'where': greatBound,    'args': [lonMeta], 'bounds': [294.0]} \
                                  , {'where': lessBound,     'args': [latMeta], 'bounds': [ 25.0]} \
                                  , {'where': greatBound,    'args': [latMeta], 'bounds': [ 50.0]} \
                                  , {'where': notEqualBound, 'args': [varQC],   'bounds': goodFlags, \
                                     'apply_to': nonObsDiags} ] \
                               , 'labels': ['CONUS'] } \
                 , 'EUROPE':    nullBinDesc \
                 , 'E_EUROPE':  nullBinDesc \
                 , 'NAMERICA':  nullBinDesc \
                 , 'PACIFIC':   nullBinDesc \
                 , 'SAMERICA':  nullBinDesc \
                 , 'SE_ASIA':   nullBinDesc \
                 , 'S_ASIA':    nullBinDesc \
                  } \
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g., 
#  , 'polygon':   { 'variables': [varQC,lonMeta,latMeta] \
#                 , 'CONUS':    { 'filters': [ \
#                                  {'where': outsideRegion, 'args': [lonMeta,latMeta], 'bounds': ['CONUS']}], \
#                               , 'labels': ['CONUS'] } \
#                 , 'EUROPE':   { 'filters': [ \
#                                 {'where': outsideRegion, 'args': [lonMeta,latMeta], 'bounds': ['EUROPE']}], \
#                               , 'labels': ['EUROPE'] } \
#                  } \
}


#=====================
# ObsSpace definitions
#=====================

nullBinGrps = [[miss_s,[]]]
profPressBinGrps = [ ['qc',['good','bad']], \
                     ['latband',['NAMED','NUMERICAL']], \
                     ['box',['CONUS']], \
                     ['pressure',['default']] \
                   ]
profAltBinGrps   = [ ['qc',['good','bad']], \
                     ['latband',['NAMED','NUMERICAL']], \
                     ['box',['CONUS']], \
                     ['altitude',['default']] \
                   ]
radianceBinGrps  = [ ['qc',['good','bad']], \
                     ['latband',['NAMED','NUMERICAL']], \
                     ['box',['CONUS']] \
                   ]

nullObsSpaceInfo = {'ObsSpaceGrp': miss_s,    'process': 0, 'binGrps': nullBinGrps }

profile_s = 'profile'
radiance_s = 'radiance'

# columns: ObsSpace name (YAML)    ObsSpaceGrp              process?                  binGrps
ObsSpaceDict_base = { \
    'sondes':                {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps } \
  , 'aircraft':              {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps } \
  , 'satwind':               {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profPressBinGrps } \
  , 'gnssroref':             {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   } \
  , 'gnssrobndropp1d':       {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   } \
  , 'gnssro':                {'ObsSpaceGrp': profile_s,  'process': True, 'binGrps': profAltBinGrps   } \
  , 'airs_aqua':             {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [1,6,7] } \
  , 'amsua_n15':             {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [5,6,7,8,9] } \
  , 'amsua_n18':             {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [5,6,7,8,9] } \
  , 'amsua_n19':             {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [5,6,7,9] } \
  , 'amsua_metop-a':         {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [5,6,9] } \
  , 'amsua_metop-b':         {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [] } \
  , 'amsua_aqua':            {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [8,9] } \
  , 'amsua_n19--ch1-3,15':   {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [1,2,3,15] } \
  , 'amsua_n19--ch4-7,9-14': {'ObsSpaceGrp': radiance_s, 'process': True, 'binGrps': radianceBinGrps, \
                              'channels': [4,5,6,7,9,10,11,12,13,14] } \
  , 'cris-fsr_npp':          {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': [24,26,28,32,37,39] } \
  , 'hirs4_metop-a':         {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,16) } \
  , 'iasi_metop-a':          {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': [16,29,32,35,38,41,44] } \
  , 'mhs_n19':               {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,6) } \
  , 'seviri_m08':            {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': [5] } \
  , 'sndrd1_g15':            {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,16) } \
  , 'sndrd2_g15':            {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,16) } \
  , 'sndrd3_g15':            {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,16) } \
  , 'sndrd4_g15':            {'ObsSpaceGrp': radiance_s, 'process': False, 'binGrps': radianceBinGrps, \
                              'channels': range(1,16) } \
    }

#============================
# figure/plotting definitions
#============================

plotMarkers = ['k-*', 'b-*', 'g-*', 'r-*', 'c-*', 'm-*', \
               'k--+','b--+','g--+','r--+','c--+','m--+']

def setup_fig(nx, ny, inch_size, aspect, ybuffer):
#INPUTS
# nx - number of subplots in x direction
# ny - number of subplots in y direction
# inch_size - rough subplot size in inches
# ybuffer - whether to give extra y space for labeling
#
#OUTPUT
# fig - a new figure with standard sizing

    fig = plt.figure()

    if ybuffer:
        fig.set_size_inches(nx*inch_size,aspect*ny*inch_size)
    else:
        fig.set_size_inches(0.9*nx*inch_size,0.9*aspect*ny*inch_size)

    return(fig)

def finalize_fig(fig, filename, filetype, ybuffer):
#INPUTS
# fig - plt.figure() type
# filename - name of figure file without extension
# filetype - file extension, e.g., 'png'
# ybuffer - whether to give extra y space for labeling

    if ybuffer:
        fig.subplots_adjust(wspace=0.35,hspace=0.70)
    else:
        fig.subplots_adjust(wspace=0.35,hspace=0.40)

    if filetype == 'png':
        fig.savefig(filename+'.png', dpi=200,bbox_inches='tight')
    plt.close(fig)

def TDeltas2Seconds(x_):
    if isinstance(x_[0],dt.timedelta):
        x = []
        for xVal in x_:
            x.append(xVal.total_seconds())
        return x
    return x_

def format_x_for_dates(ax,x):
    if isinstance(x[0],dt.datetime):
        ax.xaxis.set_major_locator(DTimeLocator)
        ax.xaxis.set_major_formatter(DTimeFormatter)
        ax.xaxis.set_tick_params(rotation=30)
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

def timeTicks(x, pos):
    d = dt.timedelta(seconds=x)
    if d.seconds > 0:
       return str(d)
    else:
       return '{:d}'.format(d.days)+'d'

DTimeLocator = AutoDateLocator(interval_multiples=True)
DTimeFormatter = DateFormatter('%m-%d_%HZ')
TDeltaFormatter = matplotlib.ticker.FuncFormatter(timeTicks)


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

def isfloat(value):
#PURPOSE determine if value can be converted to float
  try:
    float(value)
    return True
  except ValueError:
    return False

def isint(value):
#PURPOSE determine if value can be converted to int
  try:
    int(value)
    return True
  except ValueError:
    return False


#========================
# aggregation definitions
#========================

#Ordered list of statistics available in ASCII stats_* files that can be aggregated
aggregatableFileStats = ['Count','Mean','RMS','STD','Min','Max']
allFileStats = aggregatableFileStats

def aggStatsSeries(x):
#PURPOSE: aggregate DataFrame group of statistics
# INPUT: x - group of statistics that include all aggregatableFileStats
# OUTPUT: pandas Series of aggregated values for aggregatableFileStats
    y = {}

    y['Count'] = x['Count'].sum()

    y['Mean'] = (x['Mean'] * x['Count'].astype(float)).sum() / x['Count'].astype(float).sum()

    y['RMS'] = np.sqrt((x['RMS']**2 * x['Count'].astype(float)).sum() / x['Count'].astype(float).sum())

    y['STD'] = np.sqrt( ( ((x['STD'] ** 2 + x['Mean'] ** 2) \
                          * x['Count'].astype(float)).sum()  \
                          / x['Count'].astype(float).sum() ) \
                        - y['Mean'] ** 2 )
    # Pooled variance Formula as described here:
    #  https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-or-more-groups-given-known-group-varianc

    y['Min'] = x['Min'].min()

    y['Max'] = x['Max'].max()

    return pd.Series(y, index=aggregatableFileStats)


#================================
# parallel processing definitions
#================================

def proc_print(nproc,myproc,msg):
    # addendum to normal print()
    # add prefix of processor rank (myproc) when #of processors (nproc) > 0
    prefix = ""
    if nproc > 1:
        prefix = "p="+"{:d}".format(myproc)+": "
    print(prefix+msg)

def print_par_args_usage(thisPyScript):
    print ('Either 0 or 2 arguments are required.')
    print ('usage: python '+thisPyScript+'')
    print ('   or  python '+thisPyScript+' -n <NP> -i <proc>')
    print ('   or  python '+thisPyScript+' --nprocs <NP> --iproc <proc>')
    print ('')
    print ('  where <NP> is the total number of processors')
    print ('  and <proc> is the processor rank, starting at 1')
    print ('')
    print('   E.g., using GNU parallel:')
    print ('')
    print('   parallel -j<NP> --plus "python '+thisPyScript+' -n {##} -i {#} >& diags{#}.log" ::: `seq <NP>`')
    os._exit(-1)

def par_args(argv):
    opts, args = getopt.getopt(argv[1:], 'n:i:', ['nprocs', 'iproc'])
    if len(opts) == 0:
        return [1,0]
    elif len(opts) == 2:
        nproc = int(opts[0][1])
        myproc = int(opts[1][1]) - 1
        if myproc < 0:
           print_par_args_usage(argv[0])
        return [nproc, myproc]
    else:
        print_par_args_usage(argv[0])

#================================
#================================

def main():
    print ('This is not a runnable program.')
    os._exit(0)

if __name__ == '__main__': main()
