import os, sys, getopt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, AutoDateLocator
#                             HOURLY, rrulewrapper, RRuleLocator
import datetime as dt
import numpy as np
import pandas as pd

miss_f = -88888.8
miss_i = -88888
miss_s = 'null'
csvSEP = ";"

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
  , 'latitude':               [ 'deg',   'lat' ] \
  , 'longitude':              [ 'deg',   'lon' ] \
    }

profile_s = 'profile'
radiance_s = 'radiance'

#===================
# binning functions
#===================

def belowBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.less(x,bound))
    return mask

def aboveBound(xVecs, bound):
    mask = np.full((len(xVecs[0])),False)
    for x in xVecs:
        mask = np.logical_or(mask, np.greater_equal(x,bound))
    return mask

#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g.,
#def outsideRegion(LONLAT,REGION_NAME): 
#    Note: depending on shape file definitions, LON may need to be -180 to 180 instead of 0 to 360
#
#    shp = READ_SHAPEFILE(REGION_NAME)
#    lons = x[0]
#    lats = x[1]
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
nullBinDesc = { 'filters': [], 'values': [miss_s] }
allBins = 'ALL'

#========================================
# bin dictionary used for all bin groups
#========================================

binGrpDict = { \
#    #name            binVars ... followed by dictionary defitions of each binDesc
    'pressure':  { 'variables': ['air_pressure@MetaData'] \
                 , 'default': { 'filters': [\
                                  {'where': belowBound, 'args': [0], 'bounds': pressureMinBounds} \
                                , {'where': aboveBound, 'args': [0], 'bounds': pressureMaxBounds} ] \
                              , 'values': pressureBinStrVals } \
                 } \
  , 'altitude':  { 'variables': ['altitude@MetaData'] \
                 , 'default': { 'filters': [\
                                  {'where': belowBound, 'args': [0], 'bounds': altitudeMinBounds} \
                                , {'where': aboveBound, 'args': [0], 'bounds': altitudeMaxBounds} ] \
                              , 'values': altitudeBinStrVals } \
                 } \
  , 'latband':   { 'variables': ['latitude@MetaData'] \
                 , 'NAMED':     { 'filters': [\
                                  {'where': belowBound, 'args': [0], 'bounds': namedLatBandsMinBounds} \
                                , {'where': aboveBound, 'args': [0], 'bounds': namedLatBandsMaxBounds} ] \
                              , 'values': namedLatBandsStrVals } \
                 , 'NUMERICAL': { 'filters': [\
                                  {'where': belowBound, 'args': [0], 'bounds': latitudeMinBounds} \
                                , {'where': aboveBound, 'args': [0], 'bounds': latitudeMaxBounds} ] \
                              , 'values': latitudeBinStrVals } \
                 } \
  , 'box':       { 'variables': ['longitude@MetaData','latitude@MetaData'] \
                 , 'AFRICA':    nullBinDesc \
                 , 'ATLANTIC':  nullBinDesc \
                 , 'AUSTRALIA': nullBinDesc \
                 , 'CONUS':    { 'filters': [\
                                   {'where': belowBound, 'args': [0], 'bounds': [234.0]} \
                                 , {'where': aboveBound, 'args': [0], 'bounds': [294.0]} \
                                 , {'where': belowBound, 'args': [1], 'bounds': [ 25.0]} \
                                 , {'where': aboveBound, 'args': [1], 'bounds': [ 50.0]} ] \
                               , 'values': ['CONUS'] } \
                 , 'EUROPE':    nullBinDesc \
                 , 'E_EUROPE':  nullBinDesc \
                 , 'NAMERICA':  nullBinDesc \
                 , 'PACIFIC':   nullBinDesc \
                 , 'SAMERICA':  nullBinDesc \
                 , 'SE_ASIA':   nullBinDesc \
                 , 'S_ASIA':    nullBinDesc \
                  } \
#TODO: use shapefiles/polygons to describe geographic regions instead of lat/lon boxes, e.g., 
#  , 'polygon':   { 'variables': ['longitude@MetaData','latitude@MetaData'] \
#                 , 'CONUS':  {'filters': [ \
#                                {'where': outsideRegion, 'args': [0,1], 'bounds': ['CONUS']}], \
#                              'values': ['CONUS'] } \
#                 , 'EUROPE': {'filters': [ \
#                                {'where': outsideRegion, 'args': [0,1], 'bounds': ['EUROPE']}], \
#                              'values': ['EUROPE'] } \
#                  } \
}


#=====================
# ObsSpace definitions
#=====================

nullBinKeys = [[miss_s,[]]]
profPressBinKeys = [ ['pressure',['default']], ['latband',['NAMED','NUMERICAL']], ['box',['CONUS']] ]
profAltBinKeys =   [ ['altitude',['default']], ['latband',['NAMED','NUMERICAL']], ['box',['CONUS']] ]
radianceBinKeys =  [ ['latband',['NAMED','NUMERICAL']], ['box',['CONUS']] ]

# columns: ObsSpace name       ObsSpaceGrp    binGrp         process?
ObsSpaceDict_base = { \
    'sonde':                 [ profile_s,    1, profPressBinKeys ] \
  , 'aircraft':              [ profile_s,    1, profPressBinKeys ] \
  , 'satwind':               [ profile_s,    1, profPressBinKeys ] \
  , 'gnssroref':             [ profile_s,    1, profAltBinKeys   ] \
  , 'gnssrobndropp1d':       [ profile_s,    1, profAltBinKeys   ] \
  , 'amsua_n15':             [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_n18':             [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_n19':             [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_metop-a':         [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_metop-b':         [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_aqua':            [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_n19--ch1-3,15':   [ radiance_s,   1, radianceBinKeys  ] \
  , 'amsua_n19--ch4-7,9-14': [ radiance_s,   1, radianceBinKeys  ] \
    }
    #Note, refractivity: we plot RMSE of OMB/O and OMA/O; refractivity unit: N-unit
    #Note, bending_angle: we plot RMSE of OMB/O and OMA/O; bendibinVar == 'altitude':



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
        fig.subplots_adjust(wspace=0.35,hspace=0.65)
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
