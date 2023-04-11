#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy
import logging
from netCDF4 import Dataset
import numpy as np
import os
os.environ['NUMEXPR_NUM_THREADS'] = '1'
import pandas as pd

_logger = logging.getLogger(__name__)

#==========================================================
# utilities for statistics, aggregation, and bootstrapping
#==========================================================

# single-precision for floating point storage
rKindStore = np.float32
rKindNC = 'f4'

# double-precision for floating point storage
#rKindStore = np.float64
#rKindNC = 'f8'

# double-precision for floating point calculation
rKindCompute = np.float64

# double-precision for robustness
iKind = np.int64

#============
# statistics
#============

# ddof - Delta Degrees of Freedom used in standard deviation calculations
# Note: must be identical across experiments for the two stages:
# + statistics generation
# + plot generation
#ddof = int(1)

#Ordered list of statistics available in ASCII stats_* files...
# (1) that can be aggregated
aggregatableFileStats = [
  'Count', 'Mean', 'MS', 'RMS', 'STD',
  'm3', 'm4', 'Skew', 'ExcessKurtosis',
  'Min', 'Max',
]
nonCountAggregatableFileStats = deepcopy(aggregatableFileStats)
nonCountAggregatableFileStats.remove('Count')

sigmaStatistics = ['RMS']
posSemiDefiniteStats = ['Count', 'MS', 'RMS', 'STD']

allFileStats = aggregatableFileStats
nonCountFileStats = deepcopy(allFileStats)
nonCountFileStats.remove('Count')

statDtypes = {
  'Count': iKind,
  'Mean': rKindStore,
  'MS': rKindStore,
  'RMS': rKindStore,
  'STD': rKindStore,
  'm3': rKindStore,
  'm4': rKindStore,
  'Skew': rKindStore,
  'ExcessKurtosis': rKindStore,
  'Min': rKindStore,
  'Max': rKindStore,
}

statDtypesCompute = {
  'Count': iKind,
  'Mean': rKindCompute,
  'MS': rKindCompute,
  'RMS': rKindCompute,
  'STD': rKindCompute,
  'm3': rKindCompute,
  'm4': rKindCompute,
  'Skew': rKindCompute,
  'ExcessKurtosis': rKindCompute,
  'Min': rKindStore,
  'Max': rKindStore,
}

# (2) that can be sampled with bootstrap
sampleableAggStats = ['Count','Mean','MS','RMS']

bootStrapStats = []
for s in sampleableAggStats:
    if s not in ['Count', 'Mean']: bootStrapStats.append(s)

# pandas indices
fileStatAttributes = ['DiagSpaceGrp','varName','varUnits','diagName',
                      'binMethod','binVar','binVal','binUnits']


###############################################################################
def calcStats(x_):
# INPUTS:
# x[np.array(float)] - 1-d array of floating point values for which
#  statistics will be calculated

# RETURNS:
# y[dict] - population statistics, accounting for NaN values

  #Only include non-NaN values in statistics
  y = {}
  y['Count']  = np.isfinite(x_).sum()

  for statName in nonCountFileStats:
    y[statName] = statDtypes[statName](np.NaN)

  if y['Count'] > 0:
    x = np.asarray(x_).flatten()
    x = x[np.isfinite(x)]
    y['Mean'] = np.nanmean(x)
    y['MS']   = np.nanmean(np.square(x))

    # raw moments
    m1 = y['Mean']
    m2 = y['MS']
    m3 = np.nanmean(np.power(x, 3))
    m4 = np.nanmean(np.power(x, 4))
    y['m3'] = m3
    y['m4'] = m4

    # central moments
    c2 = centralMoment2(m1, m2)
    c3 = centralMoment3(m1, m2, m3)
    c4 = centralMoment4(m1, m2, m3, m4)

    # RMS
    y['RMS']  = np.sqrt(m2)

    if np.isfinite(c2) and c2 >= 0.:
      # STD
      y['STD']  = np.sqrt(c2)

      if c2 > 0:
        # higher-order standardized moments, which are normalized by a power of c2
        y['Skew'] = c3 / np.power(c2, 1.5)
        y['ExcessKurtosis'] = c4 / np.power(c2, 2) - 3.

    #y['Skew'] = np.nanmean(np.power(x - y['Mean'], 3)) / np.power(y['STD'], 3)
    #y['ExcessKurtosis'] = np.nanmean(np.power(x - y['Mean'], 4)) / np.power(y['STD'], 4) - 3.

    y['Min']  = np.nanmin(x)
    y['Max']  = np.nanmax(x)

#  else:
#    for stat in nonCountFileStats:
#      y[stat] = np.NaN

  return y


# central moment 2 in terms of raw moments 1, 2
def centralMoment2(m1, m2):
  return m2 - m1**2

# central moment 3 in terms of raw moments 1, 2, 3
def centralMoment3(m1, m2, m3):
  # m3 - 3*m1*(m2 - m1**2) - m1**3
  # m3 - 3*m1*m2 + 3*m1*m1**2 - m1**3
  # m3 - 3*m1*m2 + 3*m1**3 - m1**3
  a = np.array([m3, -3*m1*m2, 2*(m1**3)])
  a = np.sort(a)
  return np.sum(a)

# central moment 4 in terms of raw moments 1, 2, 3, 4
def centralMoment4(m1, m2, m3, m4):
  # m4 - 4*m1*m2 + 6*(m1**2)*m2 - 3*(m1**4)
  a = np.array([m4, -4*m1*m2, 6*(m1**2)*m2, -3*(m1**4)])
  a = np.sort(a)
  return np.sum(a)

###############################################################################
class BinnedStatisticsFile:
  prefix = 'stats_'
  hdf = '.h5'
  nc = '.nc'
  def __init__(self, statSpace=None, appIdentifier=None, DiagSpace=None, directory='.'):

    if statSpace is not None:
      self.__statSpace = statSpace

    elif DiagSpace is not None:
      self.__statSpace = DiagSpace
      if appIdentifier is not None and appIdentifier != '':
        self.__statSpace = appIdentifier+"_"+self.__statSpace

    else:
      _logger.error("\n\nERROR: missing statSpace constructor info")

    self.__directory = directory
    self.__file = self.fileName()

  def fileName(self):
    return self.__directory+'/'+self.prefix+self.__statSpace

  @staticmethod
  def naturalName(s):
    return s.replace('-', '__')

  def write(self, statsDict):
    '''
    Convert the statsDict to a pandas.DataFrame, then write to hdf using the df's native method
    The dataset is hierarchical.  Storing as an hdf5 file is very efficient.
    '''
    file = self.__file+self.hdf
    if os.path.exists(file):
      os.remove(file)

    # convert dict to df
    df = pd.DataFrame.from_dict(statsDict)

    # separate index from data
    df.set_index(fileStatAttributes, inplace=True)
    df.sort_index(inplace=True)

    # write to hdf
    df.to_hdf(file, key=self.naturalName(self.__statSpace), complevel=6)

  def read(self):
    file = self.__file+self.hdf
    if os.path.exists(file):
      df = pd.read_hdf(file, key=self.naturalName(self.__statSpace))

      # convert index back to columns
      df.reset_index(inplace=True)

      # convert df to dict
      statsDict = df.to_dict('list')

      return statsDict

    # can still read old netcdf files if they are available
    # TODO: fully deprecate old nc format
    file = self.__file+self.nc
    if os.path.exists(file):
      _logger.warning(
      '''
      WARNING: BinnedStatisticsFile.read; nc format will be deprecated soon in
      favor of hdf. Please consider replacing old nc files with h5 files by
      re-running the appropriate stats file generator.
      ''')
      ncid = Dataset(file, 'r')
      ncid.set_auto_mask(False)
      statsDict = {}
      for attribName in fileStatAttributes:
        statsDict[attribName] = np.asarray(ncid.variables[attribName][:])
      for statName in allFileStats:
        statsDict[statName] = np.asarray(ncid.variables[statName][:], statDtypes[statName])
      ncid.close()

      return statsDict


#===================
# basic aggregation
#===================

###############################################################################
def aggStatsDF(x):
#PURPOSE: aggregate DataFrame containing aggregatableFileStats
# INPUT: x[pd.DataFrame] - subpopulation statistics
# RETURNS: y[pd.Series] - aggregated statistics

  #converting to dictionary first speeds up the memory access
  y = aggStatsDict(x.to_dict('list'))

  return pd.Series(y, index=aggregatableFileStats)

def aggStatsDict(x_): #, stats = aggregatableFileStats):
#PURPOSE: aggregate Dictionary containing aggregatableFileStats
# INPUT: x[dict] - subpopulation statistics
# RETURNS: y[dict] - aggregated statistics

  counts = np.array(x_['Count'])
  mask = np.isfinite(counts)
  mask[mask] = np.greater(counts[mask], 0)

  # convert lists to masked np.array to enable math functions
  x = {}
  for stat in aggregatableFileStats:
    x[stat] = np.array(x_[stat])[mask].astype(statDtypesCompute[stat])

  y = {}
  for stat in nonCountAggregatableFileStats:
    y[stat] = np.NaN

  ## Count
  Count = np.int64(np.nansum(x['Count']))
  y['Count'] = Count

  if np.isfinite(Count) and Count > 0:
    # Note: arrays are sorted in ascending order before summation to avoid precision loss
    # all raw (non-central) moments are aggregated the same way
    raw = {}
    xCount = x['Count'].astype(rKindCompute)
    for rawMoment in ['Mean', 'MS', 'm3', 'm4']:
      v1_ = np.multiply(x[rawMoment].astype(rKindCompute), xCount)
      v1_ = np.sort(v1_)
      raw[rawMoment] = np.nansum(v1_) / rKindCompute(Count)
      y[rawMoment] = rKindStore(raw[rawMoment])

    # store for easy re-use
    m1 = raw['Mean']
    m2 = raw['MS']
    m3 = raw['m3']
    m4 = raw['m4']
    c2 = centralMoment2(m1, m2)
    c3 = centralMoment3(m1, m2, m3)
    c4 = centralMoment4(m1, m2, m3, m4)

    ## RMS
    y['RMS'] = rKindStore(np.sqrt(m2))

#    # Pooled variance Formula as described here:
#    #  https://en.wikipedia.org/wiki/Pooled_variance#Sample-based_statistics
#    #  https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-or-more-groups-given-known-group-varianc
#    v1_ = np.multiply(np.square(x['STD']), np.subtract(x['Count'].astype(rKindCompute), ddof))
#    v2_ = np.multiply(np.square(x['Mean']), x['Count'].astype(rKindCompute))
#    v12 = np.nansum(np.sort(np.append(v1_, v2_)))
#    v3 = np.square(m1) * rKindCompute(Count)
#    if v12 >= v3:
#      y['STD'] = rKindStore(np.sqrt((v12 - v3) / rKindCompute(Count - ddof)))

    if np.isfinite(c2) and c2 >= 0.:
      ## STD
      y['STD'] = rKindStore(np.sqrt(c2))

      ## higher-order standardized moments, which are normalized by a power of c2
      if c2 > 0:
        y['Skew'] = rKindStore(c3 / np.power(c2, 1.5))
        y['ExcessKurtosis'] = rKindStore(c4 / np.power(c2, 2) - 3.)

    ## Min
    y['Min'] = np.nanmin(x['Min'])

    ## Max
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
    return np.sqrt(np.nanmean(np.square(x), axis=axis))


################################################################################
def bootStrapVector(X, alpha=0.05, n_samples=8000, weights=None):
# PURPOSE: compute bootstrap confidence intervals on vector of data
#
# INPUTS:
# X[np.array] - values to be bootstrapped.  For linear quantities (e.g., Mean, MS),
#   these can also be subpopulation quantities, so long as the weights are set accordingly.
#   See bootStrapClusterFunc for an example.
# alpha[float] - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n[int or list(int) or array(int)] - number of bootstrap samples, optional
# weights[np.array] - same length as X, adjusts re-sampling rates for each value, optional

# RETURNS:
# ciVals - dictionary object containing mean, CI min and CI max

  if type(n_samples) is list:
    nsSamples = n_samples
    ciVals = defaultdict(list)
  elif (type(n_samples) is np.array
    or type(n_samples) is np.ndarray):
    nsSamples = list(n_samples)
    ciVals = defaultdict(list)
  else:
    nsSamples = [n_samples]
    ciVals = {}

  max_samples = np.max(nsSamples)

  ## number of "data points" must by > 0
  #  could be aggregated over a time series, space, or any other binning characteristic
  Ndata = len(X)

  if weights is None:
    iResample = np.random.choice(Ndata, (Ndata, max_samples))
  else:
    iResample = np.random.choice(Ndata, (Ndata, max_samples), p = weights)

  XResample = X[iResample]
  Expect = np.nanmean(XResample, axis=0)

  for nSamples in nsSamples:
    sampleVals = np.sort(Expect[0:nSamples])
    nonNaNSamples = np.isfinite(sampleVals).sum()

    iMid = int(np.around(0.5 * rKindCompute(nonNaNSamples)))
    iLeft = int(np.around(0.5 * alpha * rKindCompute(nonNaNSamples)))
    iRight = int(np.around((1 - 0.5 * alpha) * rKindCompute(nonNaNSamples)))

    if len(nsSamples) == 1:
      ciVals[cimean] = rKindStore(sampleVals[iMid])
      ciVals[cimin] = rKindStore(sampleVals[iLeft])
      ciVals[cimax] = rKindStore(sampleVals[iRight])
    else:
      ciVals[cimean].append(rKindStore(sampleVals[iMid]))
      ciVals[cimin].append(rKindStore(sampleVals[iLeft]))
      ciVals[cimax].append(rKindStore(sampleVals[iRight]))

  return ciVals


################################################################################
def bootStrapVectorRMSFunc(X, Y, statFunc=np.subtract,
                           alpha=0.05, n_samples=8000):
# PURPOSE: compute bootstrap confidence intervals on RMS of vector of data
#
# INPUTS:
# X, Y - intrinsic values for whole population
# statFunc - function f(RMS(X),RMS(Y)), optional, default is np.subtract
# alpha[float] - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples[int or list(int) or array(int)] - number of bootstrap samples, optional

# RETURNS:
# ciVals - dictionary object containing mean, CI min and CI max

  if type(n_samples) is list:
    nsSamples = n_samples
    ciVals = defaultdict(list)
  elif (type(n_samples) is np.array
    or type(n_samples) is np.ndarray):
    nsSamples = list(n_samples)
    ciVals = defaultdict(list)
  else:
    nsSamples = [n_samples]
    ciVals = {}

  max_samples = np.max(nsSamples)

  ## number of "data points" must by > 0
  #  could be aggregated over a time series, space, or any other binning characteristic
  Ndata = len(X)

  iResample = np.random.choice(Ndata, (Ndata, max_samples))

  XResample = rmsFunc(X[iResample], axis=0)
  YResample = rmsFunc(Y[iResample], axis=0)

  Expect = statFunc(XResample,YResample)

  for nSamples in nsSamples:
    sampleVals = np.sort(Expect[0:nSamples])
    nonNaNSamples = np.isfinite(sampleVals).sum()

    iMid = int(np.around(0.5 * rKindCompute(nonNaNSamples)))
    iLeft = int(np.around(0.5 * alpha * rKindCompute(nonNaNSamples)))
    iRight = int(np.around((1 - 0.5 * alpha) * rKindCompute(nonNaNSamples)))

    if len(nsSamples) == 1:
      ciVals[cimean] = rKindStore(sampleVals[iMid])
      ciVals[cimin] = rKindStore(sampleVals[iLeft])
      ciVals[cimax] = rKindStore(sampleVals[iRight])
    else:
      ciVals[cimean].append(rKindStore(sampleVals[iMid]))
      ciVals[cimin].append(rKindStore(sampleVals[iLeft]))
      ciVals[cimax].append(rKindStore(sampleVals[iRight]))

  return ciVals


################################################################################
def bootStrapAggRMSFunc(X, Y, Ns, statFunc=np.subtract,
                           alpha=0.05, n_samples=8000):
# PURPOSE: compute bootstrap confidence intervals on function of aggregated RMS
#          of two arrays of subpopulation RMSs
#
# INPUTS:
# X[np.array] - one set of RMS for multiple subpopulations
# Y[np.array] - second set of RMS for the same subpopulations
# Ns[np.array] - counts of data associated for the same subpopulations
# statFunc - function f(agg(X),agg(Y)), optional, default is np.subtract
# alpha[float] - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples[int or list(int) or array(int)] - number of bootstrap samples, optional

# RETURNS:
# ciVals - dictionary object containing mean, CI min and CI max

  if type(n_samples) is list:
    nsSamples = n_samples
    ciVals = defaultdict(list)
  elif (type(n_samples) is np.array
    or type(n_samples) is np.ndarray):
    nsSamples = list(n_samples)
    ciVals = defaultdict(list)
  else:
    nsSamples = [n_samples]
    ciVals = {}

  max_samples = np.max(nsSamples)

  # remove invalid clusters
  mask = np.logical_and(np.isfinite(X), np.isfinite(Y))
  mask = np.logical_and(mask, Ns > 0.)
  if mask.sum() == 0:
    _logger.error("\n\nERROR: no valid subpopulation data")
    os._exit(1)

  X_ = np.asarray(X[mask], dtype=rKindCompute)
  Y_ = np.asarray(Y[mask], dtype=rKindCompute)
  Ns_ = np.asarray(Ns[mask], dtype=rKindCompute)

  Ndata = len(X_) ; # number of "data points"

  iResample = np.random.choice(Ndata, (Ndata, max_samples))

  NsResample = Ns_[iResample]
  XResample = np.nansum(np.multiply(np.square(X_[iResample]), NsResample), axis=0)
  YResample = np.nansum(np.multiply(np.square(Y_[iResample]), NsResample), axis=0)

  NaggResample = np.nansum(NsResample, axis=0)
  XaggResample = np.sqrt(np.divide(XResample, NaggResample))
  YaggResample = np.sqrt(np.divide(YResample, NaggResample))

  Expect = statFunc(XaggResample, YaggResample)

  for nSamples in nsSamples:
    sampleVals = np.sort(Expect[0:nSamples])
    nonNaNSamples = np.isfinite(sampleVals).sum()

    iMid = int(np.around(0.5 * rKindCompute(nonNaNSamples)))
    iLeft = int(np.around(0.5 * alpha * rKindCompute(nonNaNSamples)))
    iRight = int(np.around((1 - 0.5 * alpha) * rKindCompute(nonNaNSamples)))

    if len(nsSamples) == 1:
      ciVals[cimean] = rKindStore(sampleVals[iMid])
      ciVals[cimin] = rKindStore(sampleVals[iLeft])
      ciVals[cimax] = rKindStore(sampleVals[iRight])
    else:
      ciVals[cimean].append(rKindStore(sampleVals[iMid]))
      ciVals[cimin].append(rKindStore(sampleVals[iLeft]))
      ciVals[cimax].append(rKindStore(sampleVals[iRight]))

  return ciVals


################################################################################
def bootStrapVectorFunc(X, Y, alpha=0.05,
                        n_samples=5000,
                        vecFuncs=[identityFunc],
                        bootFuncs=[bootStrapVector],
                        statFunc=np.subtract):
# PURPOSE: compute bootstrap confidence intervals (bootFunc) on
#          E[statFunc(vecFunc(X),vecFunc(Y))] using two vectors of data, X and Y,
#          which are unique realizations of the same whole population of data
# INPUTS:
# X[np.array] - one set of values for the whole population
# Y[np.array] - second set of values for the whole population
# alpha[float] - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples[int or list(int) or array(int)] - number of bootstrap samples, optional
# vecFuncs[list(function)] - function to apply independently to X and Y, e.g., np.square, optional
# bootFuncs[list(function)] - function to use for bootstrapping each vecFuncs member (same length), optional
# statFunc[function] - function f(vecFunc(X),vecFunc(Y)), optional, default is np.subtract

# RETURNS:
# funcCIVals - dictionary object containing median, low CI, and high CI of bootstrapped stats

  if type(n_samples) is list:
    nsSamples = n_samples
  elif (type(n_samples) is np.array
    or type(n_samples) is np.ndarray):
    nsSamples = list(n_samples)
  else:
      nsSamples = [n_samples]
  max_samples = np.max(nsSamples)

  mask = np.logical_and(np.isfinite(X), np.isfinite(Y))
  X_ = X[mask]
  Y_ = Y[mask]

  funcCIVals = {}
  for ifunc, func in enumerate(vecFuncs):
    bootFunc = bootFuncs[ifunc]
    if bootFunc is bootStrapVector:
      ciVals = bootFunc(
                 statFunc(func(X_), func(Y_)),
                 n_samples=nsSamples)
    elif bootFunc is bootStrapVectorRMSFunc:
      ciVals = bootFunc(
                 X_, Y_, statFunc,
                 n_samples=nsSamples)

    funcCIVals[ifunc] = {}
    for trait in ciTraits:
      funcCIVals[ifunc][trait] = ciVals[trait]

  return funcCIVals


################################################################################
def bootStrapClusterFunc(X, Y, alpha=0.05,
                         n_samples=5000,
                         statFunc=np.subtract,
                         statNames=sampleableAggStats):
# PURPOSE: compute bootstrap confidence intervals on
#          E[statFunc(X, Y)], where X and Y contain statistics from
#          subgroups of a whole population
# INPUTS:
# X[pd.DataFrame] - one set of aggregatableFileStats for all subgroups
# Y[pd.DataFrame] - second set of aggregatableFileStats for same subgroups
# alpha[float] - confidence interval (CI) percentile (e.g., 0.05 for 95%), optional
# n_samples[int or list(int) or array(int)] - number of bootstrap samples, optional
# statFunc[function] - function f(agg(X),agg(Y)), optional, default is np.subtract
# statNames[list(str)] - the statistics for which CI's will be produced

# RETURNS:
# statCIVals - nested dictionary object containing median, low CI, and high CI
#              of all statNames

  ## initialize output
  statCIVals = {}
  for stat in statNames:
    statCIVals[stat] = {}
    for trait in ciTraits:
      statCIVals[stat][trait] = np.NaN

  ## number of "data points" must by > 0
  #  could be aggregated over a time series, space, or any other binning characteristic
  nClust = len(X)
  if nClust > 0:

    if type(n_samples) is list:
      nsSamples = n_samples
    elif (type(n_samples) is np.array
      or type(n_samples) is np.ndarray):
      nsSamples = list(n_samples)
    else:
      nsSamples = [n_samples]
    max_samples = np.max(nsSamples)

    for stat in statNames:
      if stat == 'Count': continue
      Ns = X.loc[:,'Count'].to_numpy(dtype=rKindCompute)
      X_ = X.loc[:,stat].to_numpy(dtype=rKindCompute)
      Y_ = Y.loc[:,stat].to_numpy(dtype=rKindCompute)

      # remove zero-size and non-finite clusters
      mask = np.logical_and(np.isfinite(X_), np.isfinite(Y_))
      mask = np.logical_and(mask, Ns > 0.)
      if mask.sum() == 0: continue
      Ns = Ns[mask]
      X_ = X_[mask]
      Y_ = Y_[mask]

      if stat == 'Mean' or stat == 'MS':
        weights = np.divide(Ns,np.nansum(Ns))
        ciVals = bootStrapVector(
                   statFunc(X_,Y_),
                   weights=weights,
                   n_samples=nsSamples,
                   alpha=alpha,
        )
      elif stat == 'RMS':
        ciVals = bootStrapAggRMSFunc(
                   X_, Y_, Ns, statFunc,
                   n_samples=nsSamples,
                   alpha=alpha,
        )
      else:
        _logger.error("\n\nERROR: stat not implemented: ", stat)
        os._exit(1)

      for trait in ciTraits:
        statCIVals[stat][trait] = ciVals[trait]

  return statCIVals
