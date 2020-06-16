#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd


#==========================================================
# utilities for statistics, aggregation, and bootstrapping
#==========================================================

#============
# statistics
#============

#Ordered list of statistics available in ASCII stats_* files...
# (1) that can be aggregated
aggregatableFileStats = ['Count','Mean','MS','RMS','STD','Min','Max']

allFileStats = aggregatableFileStats

# (2) that can be sampled with bootstrap
sampleableAggStats = ['Count','Mean','MS','RMS']

fileStatAttributes = ['DiagSpaceGrp','varName','varUnits','diagName',
                      'binMethod','binVar','binVal','binUnits']

statsFilePrefix = 'stats_'

###############################################################################
def calcStats(array_f):
    # array_f (float(:)) - 1-d array of float values for which statistics should be calculated

    #Only include non-NaN values in statistics
    STATS = {}
    STATS['Count']  = len(array_f)-np.isnan(array_f).sum()
    if STATS['Count'] > 0:
        STATS['Mean'] = np.nanmean(array_f)
        STATS['MS']   = np.nanmean(np.square(array_f))
        STATS['RMS']  = np.sqrt(STATS['MS'])
        STATS['STD']  = np.nanstd(array_f)
        STATS['Min']  = np.nanmin(array_f)
        STATS['Max']  = np.nanmax(array_f)
    else:
        for stat in allFileStats:
            if stat != 'Count':
                STATS[stat] = np.NaN

    return STATS


###############################################################################
def write_stats_nc(statSpace,statsDict):

    statsFile = statsFilePrefix+statSpace+'.nc'
    if os.path.exists(statsFile):
        os.remove(statsFile)

    ncid = Dataset(statsFile, 'w', format='NETCDF4')
    ncid.description = statSpace+" diagnostic statistics"

    nrows = len(statsDict[fileStatAttributes[0]])
    ncid.createDimension('nrows', nrows)

    for attribName in fileStatAttributes:
        attribHandle = ncid.createVariable(attribName,str,'nrows')
        attribHandle[:] = np.array(statsDict[attribName], dtype=object)

    for statName in allFileStats:
        statHandle = ncid.createVariable(statName,'f4','nrows')  #'f8'
        statHandle[:] = statsDict[statName]

    ncid.close()


###############################################################################
def read_stats_nc(statsFile):

    ncid = Dataset(statsFile, 'r')
    ncid.set_auto_mask(False)

    statsDict = {}
    for attribName in fileStatAttributes:
        statsDict[attribName] = np.asarray(ncid.variables[attribName][:])

    for statName in allFileStats:
        statsDict[statName] = np.asarray(ncid.variables[statName][:])

    ncid.close()

    return statsDict


#===================
# basic aggregation
#===================

###############################################################################
def aggStatsDF(x):
#PURPOSE: aggregate DataFrame containing aggregatableFileStats
# INPUT: x - pandas DataFrame containing the subpopulation stats
# OUTPUT: y - dictionary formatted as a pandas Series containing aggregated stats

    #converting to dictionary first speeds up the memory access
    y = aggStatsDict(x.to_dict('list'))

    return pd.Series(y, index=aggregatableFileStats)

def aggStatsDict(x_): #, stats = aggregatableFileStats):
#PURPOSE: aggregate Dictionary containing aggregatableFileStats
# INPUT: x - dictionary containing the subpopulation stats
# OUTPUT: y - dictionary containing aggregated stats

    # convert lists to np.array to enable math functions
    x = {}
#    for stat in stats:
    for stat in aggregatableFileStats:

        x[stat] = np.array(x_[stat])

    y = {}

    y['Count'] = x['Count'].sum()

    for stat in aggregatableFileStats:
        if stat != 'Count': y[stat] = np.NaN

    if y['Count'] > 0:
        y['Mean'] = ( np.nansum( np.multiply(x['Mean'], x['Count'].astype(float)) )
                        / ( np.sum(x['Count']).astype(float) ) )

        y['MS'] = ( np.nansum( np.multiply(x['MS'], x['Count'].astype(float)) )
                       / np.sum(x['Count']).astype(float) )

        y['RMS'] = np.sqrt( y['MS'] )

        y['STD'] = ( np.sqrt( (
                       np.nansum( np.multiply(
                           (np.square(x['STD']) + np.square(x['Mean'])),
                           x['Count'].astype(float) ) )
                              / np.sum(x['Count']).astype(float) )
                            - np.square(y['Mean']) ) )
    # Pooled variance Formula as described here:
    #  https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-or-more-groups-given-known-group-varianc
        y['Min'] = np.nanmin(x['Min'])

        y['Max'] = np.nanmax(x['Max'])

    return y


#============================================
# bootstrapping for confidence intervals (CI)
#============================================

cimean = 'VALUE'
cimin = 'LO'
cimax = 'HI'
ciTraits = [cimean, cimin, cimax]

################################################################################
def identityFunc(x):
    return x


################################################################################
def rmsFunc(x, axis=0):
    return np.sqrt(np.nanmean(np.square(x),axis=axis))


################################################################################
def bootStrapVector(X, alpha=0.05, n_samples=8000, weights=None):
# PURPOSE: compute bootstrap confidence intervals on vector of data
#
#INPUTs:
# X       - array of values
# alpha   - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n       - number of bootstrap samples, optional
# weights - apply to X when performing mean

#OUTPUT: STATS - dictionary object containing mean, CI min and CI max

    Ndata = len(X) ; # number of "data points" (could be aggregated over a time series, or over space too)

    if type(n_samples) is list:
        nsSamples = n_samples
    elif (type(n_samples) is np.array
        or type(n_samples) is np.ndarray):
        nsSamples = list(n_samples)
    else:
        nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    if weights is None:
        iResample = np.random.choice(Ndata, (Ndata,max_samples) ) #
    else:
        iResample = np.random.choice(Ndata, (Ndata,max_samples), p = weights ) #

    XResample = X[ iResample ]
    Expect = np.nanmean(XResample,axis=0)

    STATS = {}
    for trait in ciTraits:
        STATS[trait] = []

    for nSamples in nsSamples:
        sampleVals = np.sort( Expect[0:nSamples] )
        nonNaNSamples = len(sampleVals)-np.isnan(sampleVals).sum()

        iMid = np.around( 0.5 * float(nonNaNSamples) ).astype(int)
        iLeft = np.around( 0.5 * alpha * float(nonNaNSamples) ).astype(int)
        iRight = np.around( (1 - 0.5 * alpha) * float(nonNaNSamples) ).astype(int)

        STATS[cimean].append(sampleVals[iMid])
        STATS[cimin].append(sampleVals[iLeft])
        STATS[cimax].append(sampleVals[iRight])

    return STATS


################################################################################
def bootStrapVectorRMSFunc(X, Y, statFunc=np.subtract,
                           alpha=0.05, n_samples=8000):
# PURPOSE: compute bootstrap confidence intervals on RMS of vector of data
#
#INPUTs:
# X, Y     - arrays of values for whole population
# statFunc - function f(RMS(X),RMS(Y)), optional, default is np.subtract
# alpha    - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n        - number of bootstrap samples, optional

#OUTPUT: STATS - dictionary object containing mean, CI min and CI max

    Ndata = len(X) ; # number of "data points" (could be aggregated over a time series, or over space too)

    if type(n_samples) is list:
        nsSamples = n_samples
    elif (type(n_samples) is np.array
        or type(n_samples) is np.ndarray):
        nsSamples = list(n_samples)
    else:
        nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    iResample = np.random.choice(Ndata, (Ndata,max_samples) )

    XResample = rmsFunc(X[ iResample ], axis=0)
    YResample = rmsFunc(Y[ iResample ], axis=0)

    Expect = statFunc(XResample,YResample)

    STATS = {}
    for trait in ciTraits:
        STATS[trait] = []

    for nSamples in nsSamples:
        sampleVals = np.sort( Expect[0:nSamples] )
        nonNaNSamples = len(sampleVals)-np.isnan(sampleVals).sum()

        iMid = np.around( 0.5 * float(nonNaNSamples) ).astype(int)
        iLeft = np.around( 0.5 * alpha * float(nonNaNSamples) ).astype(int)
        iRight = np.around( (1 - 0.5 * alpha) * float(nonNaNSamples) ).astype(int)

        STATS[cimean].append(sampleVals[iMid])
        STATS[cimin].append(sampleVals[iLeft])
        STATS[cimax].append(sampleVals[iRight])

    return STATS


################################################################################
def bootStrapAggRMSFunc(X, Y, Ns, statFunc=np.subtract,
                           alpha=0.05, n_samples=8000):
# PURPOSE: compute bootstrap confidence intervals on aggregated RMS of vector of RMS of subpopulations
#
# X, Y     - arrays of RMS for multiple subpopulations
# statFunc - function f(agg(X),agg(Y)), optional, default is np.subtract
# alpha    - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n        - number of bootstrap samples, optional

#OUTPUT: STATS - dictionary object containing mean, CI min and CI max

    if type(n_samples) is list:
        nsSamples = n_samples
    elif (type(n_samples) is np.array
        or type(n_samples) is np.ndarray):
        nsSamples = list(n_samples)
    else:
        nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    # remove zero-size clusters
    X_ = []
    Y_ = []
    Ns_ = []
    for i, n in enumerate(Ns):
        if n > 0.0:
            X_.append(X[i])
            Y_.append(Y[i])
            Ns_.append(n)

    X_ = np.asarray(X_)
    Y_ = np.asarray(Y_)
    Ns_ = np.asarray(Ns_)

    Ndata = len(X_) ; # number of "data points"

    iResample = np.random.choice(Ndata, (Ndata,max_samples) )

    NsResample = Ns_[ iResample ]
    XResample = np.nansum(np.multiply(np.square(X_[ iResample ]), NsResample),axis=0)
    YResample = np.nansum(np.multiply(np.square(Y_[ iResample ]), NsResample),axis=0)

    NaggResample = np.nansum(NsResample,axis=0)
    XaggResample = np.sqrt(np.divide(XResample, NaggResample))
    YaggResample = np.sqrt(np.divide(YResample, NaggResample))

    Expect = statFunc(XaggResample,YaggResample)

    STATS = {}
    for trait in ciTraits:
        STATS[trait] = []

    for nSamples in nsSamples:
        sampleVals = np.sort( Expect[0:nSamples] )
        nonNaNSamples = len(sampleVals)-np.isnan(sampleVals).sum()

        iMid = np.around( 0.5 * float(nonNaNSamples) ).astype(int)
        iLeft = np.around( 0.5 * alpha * float(nonNaNSamples) ).astype(int)
        iRight = np.around( (1 - 0.5 * alpha) * float(nonNaNSamples) ).astype(int)

        STATS[cimean].append(sampleVals[iMid])
        STATS[cimin].append(sampleVals[iLeft])
        STATS[cimax].append(sampleVals[iRight])

    return STATS


################################################################################
def bootStrapVectorFunc(X, Y, alpha=0.05,
                        n_samples=5000,
                        vecFuncs=[identityFunc],
                        bootFuncs=[bootStrapVector],
                        statFunc=np.subtract):
# PURPOSE: compute bootstrap confidence intervals on
#          E[statFunc(vecFunc(X),vecFunc(Y))]
#          using vectors of data, X and Y, that contain the whole population
#INPUTS:
# X, Y      - vectors of values for the whole population
# alpha     - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples - number of bootstrap samples, optional, default==5000
#             can either be a scalar or a list of values
# vecFuncs  - function to apply independently to X and Y, e.g., np.square, optional
# statFunc  - function f(vecFunc(X),vecFunc(Y)), optional, default is np.subtract

#OUTPUTS:
# STATS - dictionary object containing median, low CI, and high CI of bootstrapped stats

    N = len(X) ; # number of "data points" (could be aggregated over a time series, or over space too)

    if type(n_samples) is list:
        nsSamples = n_samples
    elif (type(n_samples) is np.array
        or type(n_samples) is np.ndarray):
        nsSamples = list(n_samples)
    else:
        nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    STATS = {}
    for ifunc, func in enumerate(vecFuncs):
        bootFunc = bootFuncs[ifunc]
        if bootFunc is bootStrapVector:
            ciVals = bootFunc(
                         statFunc(func(X), func(Y)),
                         n_samples=nsSamples)
        elif bootFunc is bootStrapVectorRMSFunc:
            ciVals = bootFunc(
                         X, Y, statFunc,
                         n_samples=nsSamples)

        STATS[ifunc] = {}
        for trait in ciTraits:
            STATS[ifunc][trait] = ciVals[trait]

    return STATS


################################################################################
def bootStrapClusterFunc(X, Y, alpha=0.05,
                         n_samples=5000,
                         statFunc=np.subtract,
                         statNames=sampleableAggStats):
# PURPOSE: compute bootstrap confidence intervals on
#          E[calcStats(X)] or E[statFunc(calcStats(X),calcStats(Y))]
#          using Stats from subgroups of a whole population
#INPUTS:
# X, Y      - pandas DataFrames containing aggregatable stats for all subgroups
# alpha     - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples - number of bootstrap samples, optional, default==5000
#             can either be a scalar or a list of values
# statFunc  - function, optional, default is np.subtract
# statNames - the stats for which CI's will be produced

#OUTPUTS:
# STATS   - dictionary object containing median, low CI, and high CI of bootstrapped stats

    nClust = len(X) ; # number of "data points" (could be aggregated over a time series, or over space too)

    if type(n_samples) is list:
        nsSamples = n_samples
    elif (type(n_samples) is np.array
        or type(n_samples) is np.ndarray):
        nsSamples = list(n_samples)
    else:
        nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    STATS = {}
    for stat in statNames:
        if stat == 'Count': continue
        Ns = X.loc[:,'Count'].to_numpy().astype(float)
        X_ = X.loc[:,stat].to_numpy()
        Y_ = Y.loc[:,stat].to_numpy()

        if any(Ns > 0.0):
            if stat == 'Mean' or stat == 'MS':
                weights = np.divide(Ns,np.nansum(Ns))
                ciVals = bootStrapVector(
                             statFunc(X_,Y_),
                             weights=weights,
                             n_samples=nsSamples)
            elif stat == 'RMS':
                ciVals = bootStrapAggRMSFunc(
                             X_, Y_, Ns, statFunc,
                             n_samples=nsSamples)
            else:
                print("\n\nERROR: stat not implemented: ", stat)
                os._exit(1)
        else:
            ciVals = {}
            for trait in ciTraits: ciVals[trait] = [np.NaN]


        STATS[stat] = {}
        for trait in ciTraits:
            STATS[stat][trait] = ciVals[trait]

    return STATS

