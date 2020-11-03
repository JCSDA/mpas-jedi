#!/usr/bin/env python3

import datetime as dt
import os

########################################################################
'''
This module is used to configure statistical analyses.  Those analyses
can be intialized either by directly executing AnalyzeStats.py or
by submitting a series of jobs for multiple DiagSpaces using
SpawnAnalyzeStats.py and AnalyzeStats.csh.

Command-line examples:
----------------------
+ Carry out analyses for all DiagSpaces that contain "amsua"

    python AnalyzeStats.py -d amsua

+ Use 12 processes to carry out analyses for all
  DiagSpaces with anGroup == "conv"

    python AnalyzeStats.py -n 12 -g conv

+ Use 12 processes to carry out analyses for all
  DiagSpaces that contain "abi", and use 30 processes
  for reading the StatisticsDatabase

    python AnalyzeStats.py -n 12 -r 30 -d abi

+ Get info about more options

    python AnalyzeStats.py --help

Job-submission examples:
------------------------
+ Spawn one job for each DiagSpace that is enabled in config
  using the anGroupConfig specified therein

    python SpawnAnalyzeStats.py

+ Specify that statistics files come from a JEDI hofx application

    python SpawnAnalyzeStats.py -app hofx

+ Spawn one job for each DiagSpace that contains "amsua"

    python SpawnAnalyzeStats.py -d amsua

+ Choose a unique job account number

    python SpawnAnalyzeStats.py -a NMMM0043

+ Get info about more options

    python SpawnAnalyzeStats.py --help

'''
########################################################################
#Select the statistics to be analyzed
# options: see su.allFileStats
# 'analysisStatistics' from individual diagnostics will override this setting
analysisStatistics = ['Count','Mean','RMS','STD']

#Select diagnostic groupings
# + each entry in the dictionary should be of the format:
#   diagnosticGroup: diagnosticNames
# + diagnosticNames is a list of diagnostics, such as
#   [diagnosticName1, diagnosticName2, etc...]
# + diagnosticGroup is up to the user, but it is best to have it match
#   the strings used for axis labeling in Analyses
# + the default behavior is to plot each diagnostic on an independent axis,
#   which will still be done for any analysis type that does not use
#   diagnosticGroupings or has maxDiagnosticsPerAnalysis < len(diagnosticNames)
diagnosticGroupings = {}
diagnosticGroupings['omm'] = ['omb', 'oma']
diagnosticGroupings['rltv_omm'] = ['rltv_omb', 'rltv_oma']

## Configure the StatisticsDatabase.StatsDB objects
dbConf = {}

## hasFCLenDir whether directory structure includes forecast length
## note: overridden when fcTDeltaLast > fcTDeltaFirst
dbConf['hasFCLenDir'] = False

## expDirectory is the top-level directory that contains data
#  from all experiments.  The environment variable, EXP_DIR, can be used
#  to specify its value externally if desired.
user = 'guerrett'
dbConf['expDirectory'] = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac')

## cntrlExpIndex is the index of the control experiment (used for DiffCI analyses)
#  in the expNames list
dbConf['cntrlExpIndex'] = 0

## expLongNames is a list of directories within expDirectory that contains
#  data from individual experiments
dbConf['expLongNames'] = []

## expNames is a list of experiment names used for database lookups and
#  analysis labels, e.g., on figures.  Make these names concise and
#  exclude spaces.
dbConf['expNames'] = []

## DAMethods is a list of secondary labels that are within the database
#  intermediate file names.  They may be the actual DA method used (e.g., 3denvar),
#  an indication of the workflow component (e.g., omb or omm), or a customized
#  label for each experiment.  DAMethods is only important for file naming and is
#  not used in the analyses
dbConf['DAMethods'] = []
## Note: refer to the StatisticsDatabase.StatsDB class for more details

## -------------------------------------------
## Append to expLongNames, expNames, DAMethods
## -------------------------------------------
## EDA w/ baseline config
#nEnsDAMembers = 20
#for mem in list(range(1,nEnsDAMembers+1)):
#  member = '{:03d}'.format(mem)
#  dbConf['expLongNames'].append(
#    'guerrett_eda_3denvar_conv_clramsua_NMEM20_120km/Verification/bg/mem'+member)
#  dbConf['expNames'].append('EDA'+member)
#  dbConf['DAMethods'].append('omm')
#dbConf['cntrlExpIndex'] = nEnsDAMembers

## 3denvar baseline (conventional + clear-sky AMSUA)
#dbConf['expLongNames'].append('guerrett_3denvar_conv_clramsua_NMEM1_120km/Verification/bg')
dbConf['expLongNames'].append('guerrett_3denvar_conv_clramsua_NMEM1_120km/Verification/fc/mean')
dbConf['expNames'].append('deterministic')
dbConf['DAMethods'].append('omm')

## eda-3denvar mean
dbConf['expLongNames'].append('guerrett_eda_3denvar_conv_clramsua_NMEM20_120km/Verification/fc/mean')
dbConf['expNames'].append('eda')
dbConf['DAMethods'].append('omm')

## eda-3denvar RTPP0.5 mean
dbConf['expLongNames'].append('guerrett_eda_3denvar_conv_clramsua_NMEM20_RTPP0.5_120km/Verification/fc/mean')
dbConf['expNames'].append('eda-RTPP0.5')
dbConf['DAMethods'].append('omm')

## -------------------------------
## Cycle times and forecast length
## -------------------------------
#First and Last CYCLE dates and increment
dbConf['firstCycleDTime'] = dt.datetime(2018,4,15,0,0,0)
dbConf['lastCycleDTime'] = dt.datetime(2018,4,18,0,0,0)
dbConf['cyTimeInc'] = dt.timedelta(hours=12)

#First and Last FORECAST durations and increment
dbConf['fcTDeltaFirst'] = dt.timedelta(days=0)
dbConf['fcTDeltaLast'] = dt.timedelta(days=0,hours=6)
dbConf['fcTimeInc'] = dt.timedelta(hours=6)

## fcDirFormats is used to declare the directory string format
#  for forecast lengths. Can include any combination of substrings
#  from the following examples:
#  "%D" (only number of days)
#  "%D_%HH:%MM:%SS"
#  "%MIN:%SEC"
#  "%s" (total seconds)
#  "%m", "%MIN" (total full minutes)
#  "%h" (total full hours)
#
# By default, all experiments use the same format string. Override
# that behavior by redefining the fcDirFormats list below.
commonFCDirFormat = "%hhr"
dbConf['fcDirFormats'] = [commonFCDirFormat]*len(dbConf['expNames'])

## statsFileSubDirs is the subdirectory within the date directory(ies)
#  that contains the statstics files for constructing the StatsDB object
# examples:
#TODO: make StatsDB search for the correct subdirectory for each experiment
# 'diagnostic_stats'
# 'diagnostic_stats/obs'
# 'diagnostic_stats/model'
commonStatsFileSubDir = 'diagnostic_stats/obs'
dbConf['statsFileSubDirs'] = [commonStatsFileSubDir]*len(dbConf['expNames'])

########################################################################
## Configure the analysisTypes to apply to the statistics
#  - below are recommendations for single/multiple forecast lengths
#  - analysisTypes can be mixed and matched as desired,
#    however some of them require nCY, nFC, or nExp > 1
#  - see the individual classes for more details (Analyses.py)
analysisTypes = []
if dbConf['fcTDeltaFirst'] == dbConf['fcTDeltaLast']:
    ## gross error analysisTypes for single forecast length
    ## -------------------------------------------------------
    ## recommended
    analysisTypes.append('CYAxisExpLines')

    ## potentially useful
    analysisTypes.append('CYAxisBinValLines')
    analysisTypes.append('CYandBinValAxes2D')

else:
    ## gross error analysisTypes for multiple forecast lengths
    ## -------------------------------------------------------
    ## recommended
    analysisTypes.append('FCAxisExpLines')
    if len(dbConf['expNames']) > 1: analysisTypes.append('FCAxisExpLinesDiffCI')

    ## potentially useful
    analysisTypes.append('CYandBinValAxes2D')
    analysisTypes.append('FCandBinValAxes2D')
    analysisTypes.append('CYAxisExpLines')
    analysisTypes.append('CYAxisFCLines')

## used to dissect gross errors in more detail
analysisTypes.append('BinValAxisProfile')
if len(dbConf['expNames']) > 1: analysisTypes.append('BinValAxisProfileDiffCI')
##Note: BinValAxisProfile* analyses work for any forecast duration

## useful for prescribing/evaluating R statistics
#analysisTypes.append('BinValAxisPDF')
#analysisTypes.append('BinValAxisStatsComposite')
#analysisTypes.append('GrossValues')

