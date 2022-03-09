#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy


# allCIErrParams is a dictionary of lookup tabls for symmetric cloud impact (SCI)
# The structure of allCIErrParams keys is as follows:
#allCIErrParams[
#    ({int: MPAS forecast resoultion},{str: bias correction type})][
#    {str: instrument name}][({int: channel number},{str: SCI method})]
allCIErrParams = {}

# Individual CIErrParams dictionaries below are generated
# automatically when AGGREGATEFC_StatsComposite figures are
# created with plot_stats_timeseries.py

# More keys can be added as they become available

## 120 km 6-hr MPAS Forecast from MPAS-JEDI analysis 2018041500 to 2018051412
## /glade/scratch/wuyl/test2/pandac/test_120km/ben_FC
## VARIATIONAL BIAS CORRECTION
CIErrParams = defaultdict(dict)
#EMPTY
allCIErrParams[('cloud_impact',120,'varbc')] = deepcopy(CIErrParams)


## CONSTANT BIAS CORRECTION
CIErrParams = defaultdict(dict)
CIErrParams['abi_g16'][ (7, 'Okamoto') ]   =  {'X': [1, 29.37], 'ERR': [3.48, 18.7]}
CIErrParams['abi_g16'][ (8, 'Okamoto') ]   =  {'X': [1, 12.54], 'ERR': [1.98, 16.04]}
CIErrParams['abi_g16'][ (9, 'Okamoto') ]   =  {'X': [1, 15.15], 'ERR': [2.05, 18.2]}
CIErrParams['abi_g16'][ (10, 'Okamoto') ]   =  {'X': [1, 17.18], 'ERR': [2.07, 19.14]}
CIErrParams['abi_g16'][ (11, 'Okamoto') ]   =  {'X': [1, 25.68], 'ERR': [2.24, 26.2]}
CIErrParams['abi_g16'][ (13, 'Okamoto') ]   =  {'X': [1, 28.04], 'ERR': [2.41, 28.7]}
CIErrParams['abi_g16'][ (14, 'Okamoto') ]   =  {'X': [1, 28.0], 'ERR': [2.44, 29.02]}
CIErrParams['abi_g16'][ (15, 'Okamoto') ]   =  {'X': [1, 27.13], 'ERR': [2.39, 28.69]}
CIErrParams['abi_g16'][ (16, 'Okamoto') ]   =  {'X': [1, 23.71], 'ERR': [2.18, 24.68]}
#For UFO YAML config:
#  x0: [1, 1, 1, 1, 1, 1, 1, 1, 1]
#  x1: [29.37, 12.54, 15.15, 17.18, 25.68, 28.04, 28.0, 27.13, 23.71]
#  err0: [3.48, 1.98, 2.05, 2.07, 2.24, 2.41, 2.44, 2.39, 2.18]
#  err1: [18.7, 16.04, 18.2, 19.14, 26.2, 28.7, 29.02, 28.69, 24.68]
allCIErrParams[('cloud_impact',120,'constant')] = deepcopy(CIErrParams)


## NO BIAS CORRECTION
CIErrParams = defaultdict(dict)
#For binning_params(bench+IR): - fits to RMS (not STD)
CIErrParams['abi_g16'][ (7, 'Okamoto') ]   =  {'X': [1, 35.91], 'ERR': [2.17, 35.92]}
CIErrParams['abi_g16'][ (8, 'Okamoto') ]   =  {'X': [1, 10.31], 'ERR': [1.79, 14.92]}
CIErrParams['abi_g16'][ (9, 'Okamoto') ]   =  {'X': [1, 12.33], 'ERR': [1.9, 16.85]}
CIErrParams['abi_g16'][ (10, 'Okamoto') ]   =  {'X': [1, 15.91], 'ERR': [1.81, 18.01]}
CIErrParams['abi_g16'][ (11, 'Okamoto') ]   =  {'X': [1, 30.81], 'ERR': [2.25, 26.86]}
CIErrParams['abi_g16'][ (12, 'Okamoto') ]   =  {'X': [11, 15.11], 'ERR': [19.81, 22.77]}
CIErrParams['abi_g16'][ (13, 'Okamoto') ]   =  {'X': [1, 28.76], 'ERR': [2.36, 25.54]}
CIErrParams['abi_g16'][ (14, 'Okamoto') ]   =  {'X': [1, 26.97], 'ERR': [2.34, 24.96]}
CIErrParams['abi_g16'][ (15, 'Okamoto') ]   =  {'X': [1, 25.5], 'ERR': [2.14, 24.59]}
CIErrParams['abi_g16'][ (16, 'Okamoto') ]   =  {'X': [1, 22.41], 'ERR': [1.94, 20.95]}
#For UFO YAML config(bench+IR):
#  x0: [1, 1, 1, 1, 1, 11, 1, 1, 1, 1]
#  x1: [35.91, 10.31, 12.33, 15.91, 30.81, 15.11, 28.76, 26.97, 25.5, 22.41]
#  err0: [2.17, 1.79, 1.9, 1.81, 2.25, 19.81, 2.36, 2.34, 2.14, 1.94]
#  err1: [35.92, 14.92, 16.85, 18.01, 26.86, 22.77, 25.54, 24.96, 24.59, 20.95]
#For binning_params(bench+IR): - fits to RMS (not STD)
CIErrParams['ahi_himawari8'][ (7, 'Okamoto') ]   =  {'X': [1, 39.38], 'ERR': [2.16, 41.63]}
CIErrParams['ahi_himawari8'][ (8, 'Okamoto') ]   =  {'X': [1, 9.75], 'ERR': [1.73, 12.58]}
CIErrParams['ahi_himawari8'][ (9, 'Okamoto') ]   =  {'X': [1, 13.7], 'ERR': [1.87, 16.16]}
CIErrParams['ahi_himawari8'][ (10, 'Okamoto') ]   =  {'X': [1, 19.44], 'ERR': [1.82, 18.94]}
CIErrParams['ahi_himawari8'][ (11, 'Okamoto') ]   =  {'X': [1, 39.48], 'ERR': [2.44, 31.21]}
CIErrParams['ahi_himawari8'][ (12, 'Okamoto') ]   =  {'X': [11, 15.69], 'ERR': [20.14, 22.76]}
CIErrParams['ahi_himawari8'][ (13, 'Okamoto') ]   =  {'X': [1, 34.36], 'ERR': [2.48, 27.38]}
CIErrParams['ahi_himawari8'][ (14, 'Okamoto') ]   =  {'X': [1, 34.35], 'ERR': [2.55, 27.58]}
CIErrParams['ahi_himawari8'][ (15, 'Okamoto') ]   =  {'X': [1, 31.21], 'ERR': [2.4, 25.15]}
CIErrParams['ahi_himawari8'][ (16, 'Okamoto') ]   =  {'X': [1, 26.07], 'ERR': [2.02, 21.14]}
#For UFO YAML config(bench+IR):
#  x0: [1, 1, 1, 1, 1, 11, 1, 1, 1, 1]
#  x1: [39.38, 9.75, 13.7, 19.44, 39.48, 15.69, 34.36, 34.35, 31.21, 26.07]
#  err0: [2.16, 1.73, 1.87, 1.82, 2.44, 20.14, 2.48, 2.55, 2.4, 2.02]
#  err1: [41.63, 12.58, 16.16, 18.94, 31.21, 22.76, 27.38, 27.58, 25.15, 21.14]

allCIErrParams[('cloud_impact',120,None)] = deepcopy(CIErrParams)


## 30 km 6-hr MPAS Forecasts from GFSANA 2018041500 to 2018042200
## /glade/scratch/wuyl/test2/pandac/test_30km_cld/FC
## VARIATIONAL BIAS CORRECTION
CIErrParams = defaultdict(dict)
#EMPTY
allCIErrParams[('cloud_impact',30,'varbc')] = deepcopy(CIErrParams)


## CONSTANT BIAS CORRECTION
CIErrParams = defaultdict(dict)
#EMPTY
allCIErrParams[('cloud_impact',30,'constant')] = deepcopy(CIErrParams)


## NO BIAS CORRECTION
CIErrParams = defaultdict(dict)
# 18FEB2022, JJG, 30 days clrama
#For binning_params(clrama):^[[0m
CIErrParams['abi_g16'][(7, 'Okamoto')]   = {'X': [0.5, 20.08], 'ERR': [1.06, 21.36]} #original x0 == 0
CIErrParams['abi_g16'][(8, 'Okamoto')]   = {'X': [0.5, 10.94], 'ERR': [0.46, 16.67]} #original x0 == 0
CIErrParams['abi_g16'][(9, 'Okamoto')]   = {'X': [0.5, 13.29], 'ERR': [0.48, 19.32]} #original x0 == 0
CIErrParams['abi_g16'][(10, 'Okamoto')]   = {'X': [0.5, 14.86], 'ERR': [0.6, 19.59]} #original x0 == 0
CIErrParams['abi_g16'][(11, 'Okamoto')]   = {'X': [0.5, 26.21], 'ERR': [0.59, 28.49]} #original x0 == 0
CIErrParams['abi_g16'][(12, 'Okamoto')]   = {'X': [3, 30.57], 'ERR': [5.0, 22.47]} #original err0 == -4.13
CIErrParams['abi_g16'][(13, 'Okamoto')]   = {'X': [0.5, 27.91], 'ERR': [0.75, 31.07]} #original x0 == 0
CIErrParams['abi_g16'][(14, 'Okamoto')]   = {'X': [0.5, 27.32], 'ERR': [0.84, 31.41]} #original x0 == 0
CIErrParams['abi_g16'][(15, 'Okamoto')]   = {'X': [0.5, 27.1], 'ERR': [0.75, 30.96]} #original x0 == 0
CIErrParams['abi_g16'][(16, 'Okamoto')]   = {'X': [0.5, 23.11], 'ERR': [0.14, 25.79]} #original x0 == 0
#For UFO YAML config(clrama):
#  x0: [0.5, 0.5, 0.5, 0.5, 0.5, 3, 0.5, 0.5, 0.5, 0.5]
#  x1: [20.08, 10.94, 13.29, 14.86, 26.21, 30.57, 27.91, 27.32, 27.1, 23.11]
#  err0: [1.06, 0.46, 0.48, 0.6, 0.59, 5.0, 0.75, 0.84, 0.75, 0.14]
#  err1: [21.36, 16.67, 19.32, 19.59, 28.49, 22.47, 31.07, 31.41, 30.96, 25.79]

#For binning_params(clrama):
CIErrParams['ahi_himawari8'][(7, 'Okamoto')]   = {'X': [0.5, 21.37], 'ERR': [1.11, 22.29]} #original x0 == 0
CIErrParams['ahi_himawari8'][(8, 'Okamoto')]   = {'X': [0.5, 10.34], 'ERR': [0.49, 15.1]} #original x0 == 0
CIErrParams['ahi_himawari8'][(9, 'Okamoto')]   = {'X': [0.5, 12.76], 'ERR': [0.54, 17.65]} #original x0 == 0
CIErrParams['ahi_himawari8'][(10, 'Okamoto')]   = {'X': [0.5, 15.85], 'ERR': [0.69, 19.03]} #original x0 == 0
CIErrParams['ahi_himawari8'][(11, 'Okamoto')]   = {'X': [0.5, 24.59], 'ERR': [0.78, 26.15]} #original x0 == 0
CIErrParams['ahi_himawari8'][(12, 'Okamoto')]   = {'X': [2, 28.76], 'ERR': [5.0, 20.92]} #original err0 == -4.94
CIErrParams['ahi_himawari8'][(13, 'Okamoto')]   = {'X': [0.5, 25.29], 'ERR': [0.92, 28.34]} #original x0 == 0
CIErrParams['ahi_himawari8'][(14, 'Okamoto')]   = {'X': [0.5, 24.66], 'ERR': [0.93, 28.61]} #original x0 == 0
CIErrParams['ahi_himawari8'][(15, 'Okamoto')]   = {'X': [0.5, 24.43], 'ERR': [0.86, 27.92]} #original x0 == 0
CIErrParams['ahi_himawari8'][(16, 'Okamoto')]   = {'X': [0.5, 21.54], 'ERR': [0.43, 23.03]} #original x0 == 0
#For UFO YAML config(clrama):
#  x0: [0.5, 0.5, 0.5, 0.5, 0.5, 2, 0.5, 0.5, 0.5, 0.5]
#  x1: [21.37, 10.34, 12.76, 15.85, 24.59, 28.76, 25.29, 24.66, 24.43, 21.54]
#  err0: [1.11, 0.49, 0.54, 0.69, 0.78, 5.0, 0.92, 0.93, 0.86, 0.43]
#  err1: [22.29, 15.1, 17.65, 19.03, 26.15, 20.92, 28.34, 28.61, 27.92, 23.03]

allCIErrParams[('cloud_impact',30,None)] = deepcopy(CIErrParams)

scaleInflation = 1.0
ABEIParams = {
  # Minamide and Zhang, 2018 values derived from ABI radiances
  'abi_g16': {
    (8): {'LambdaOverACI': 0.015*scaleInflation},
    (9): {'LambdaOverACI': 0.012*scaleInflation},
    (10): {'LambdaOverACI': 0.009*scaleInflation},
  },
  'ahi_himawari8': {
    (8): {'LambdaOverACI': 0.015*scaleInflation},
    (9): {'LambdaOverACI': 0.012*scaleInflation},
    (10): {'LambdaOverACI': 0.009*scaleInflation},
  },
}
