#!/usr/bin/env python3

import datetime as dt
import os
import StatisticsDatabase as sdb


########################################################################
'''
This module is used to configure statistical analyses.  Those analyses
can be intialized either by directly executing analyze_stats.py or
by submitting a series of jobs for multiple DiagSpaces using 
spawn_analyze_stats_jobs.py and analyze_stats_job.csh.

Command-line examples:
----------------------
+ Carry out analyses for all DiagSpaces that contain "amsua"

    python analyze_stats.py -d amsua

+ Use 12 processes to carry out analyses for all
  DiagSpaces with anGroup == "conv" 

    python analyze_stats.py -n 12 -g conv

+ Use 12 processes to carry out analyses for all
  DiagSpaces that contain "abi", and use 30 processes
  for reading the StatisticsDatabase

    python analyze_stats.py -n 12 -r 30 -d abi

+ Get info about more options

    python analyze_stats.py --help

Job-submission examples:
------------------------
+ Spawn one job for each DiagSpace that is enabled in config
  using the anGroupConfig specified therein

    python spawn_analyze_stats_jobs.py

+ Spawn one job for each DiagSpace that contains "amsua"

    python spawn_analyze_stats_jobs.py -d amsua

+ Choose a unique job account number

    python spawn_analyze_stats_jobs.py -a NMMM0043

+ Get info about more options

    python spawn_analyze_stats_jobs.py --help

'''
########################################################################
## Configure the sdb.StatsDB objects
dbConf = {}

## How many forecast lengths are included?
## one: sdb.singleFCLen
## more: sdb.multiFCLen
## - use these as templates
dbConf['plotGroup'] = sdb.multiFCLen

## expDirectory is the top-level directory that contains data
#  from all experiments.  The environment variable, EXP_DIR, can be used
#  to specify its value externally if desired.
user = 'guerrett'
dbConf['expDirectory'] = os.getenv('EXP_DIR','/glade/scratch/'+user+'/pandac/')

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

## Note: refer to the sdb.StatsDB class for more details

## cntrlExpIndex is the index of the control experiment (used for DiffCI analyses)
#  in the expNames list
dbConf['cntrlExpIndex'] = 0

if dbConf['plotGroup'] == sdb.singleFCLen:
# Warm start from PANDA-C baseline, April 2020 (Yali)
    dbConf['expLongNames'].append('OLDCODE/120km_3denvar_conv_clramsua/VF/fc-6hr-MPAS_conv_clramsua_YW')
    dbConf['expNames'].append('CNTRL')
    dbConf['DAMethods'].append('omb-6hr-MPAS_conv_clramsua')

## Clear-sky ABI + PANDA-C baseline
    dbConf['expLongNames'].append('OLDCODE/120km_3denvar_conv_clramsua_clrabi/VF/bg')
    dbConf['expNames'].append('CLRABI')
    dbConf['DAMethods'].append('omm')

## All-sky ABI + PANDA-C baseline
    dbConf['expLongNames'].append('OLDCODE/120km_3denvar_conv_clramsua_cldabi/VF/bg')
    dbConf['expNames'].append('ALLABI')
    dbConf['DAMethods'].append('omm')

    #First and Last CYCLE dates
    dbConf['firstCycleDTime'] = dt.datetime(2018,4,15,0,0,0)
    dbConf['lastCycleDTime'] = dt.datetime(2018,5,14,12,0,0)
    dbConf['cyTimeInc'] = dt.timedelta(hours=6)

    #First and Last FORECAST durations
    dbConf['fcTDeltaFirst'] = dt.timedelta(days=0)
    dbConf['fcTDeltaLast'] = dt.timedelta(days=0)
    dbConf['fcTimeInc'] = dt.timedelta(hours=24)

if dbConf['plotGroup'] == sdb.multiFCLen:
## PANDA-C baseline, April 2020 (Yali)
    dbConf['expLongNames'].append('OLDCODE/120km_omf_conv_clramsua/VF/fc-ANA_conv_clramsua_YW')
    dbConf['expNames'].append('CNTRL')
    dbConf['DAMethods'].append('omm')

## Clear-sky ABI + PANDA-C baseline
    dbConf['expLongNames'].append('OLDCODE/120km_3denvar_conv_clramsua_clrabi/VF/fc')
    dbConf['expNames'].append('CLRABI')
    dbConf['DAMethods'].append('omm')

## All-sky ABI + PANDA-C baseline
    dbConf['expLongNames'].append('OLDCODE/120km_3denvar_conv_clramsua_cldabi/VF/fc')
    dbConf['expNames'].append('ALLABI')
    dbConf['DAMethods'].append('omm')

    #First and Last CYCLE dates
    dbConf['firstCycleDTime'] = dt.datetime(2018,4,15,0,0,0)
    dbConf['lastCycleDTime'] = dt.datetime(2018,5,11,12,0,0)
    dbConf['cyTimeInc'] = dt.timedelta(hours=12)

    #First and Last FORECAST durations
    dbConf['fcTDeltaFirst'] = dt.timedelta(hours=0)
    dbConf['fcTDeltaLast'] = dt.timedelta(days=3, hours=0)
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


########################################################################
## Configure the analysisTypes to apply to the statistics
#  - below are recommendations for single/multiple forecast lengths
#  - analysisTypes can be mixed and matched as desired, 
#    however some of them require nCY, nFC, or nExp > 1
#  - see the individual classes for more details (AnalyzeStatistics.py)
analysisTypes = []
if dbConf['fcTDeltaFirst'] == dbConf['fcTDeltaLast']:
    ## gross error analysisTypes for single forecast length
    ## -------------------------------------------------------
    ## recommended
    analysisTypes.append('CYAxisExpLines')

    ## potentially useful
    # analysisTypes.append('CYAxisBinValLines')
    # analysisTypes.append('CYandBinValAxes2D')

else:
    ## gross error analysisTypes for multiple forecast lengths
    ## -------------------------------------------------------
    ## recommended
    analysisTypes.append('FCAxisExpLines')
    if len(dbConf['expNames']) > 1: analysisTypes.append('FCAxisExpLinesDiffCI')

    ## potentially useful
    # analysisTypes.append('FCandBinValAxes2D')
    # analysisTypes.append('CYAxisExpLines')
    # analysisTypes.append('CYAxisFCLines')

## used to dissect gross errors in more detail
analysisTypes.append('BinValAxisProfile')
if len(dbConf['expNames']) > 1: analysisTypes.append('BinValAxisProfileDiffCI')
##Note: BinValAxisProfile* analyses work for any forecast duration

## useful for prescribing/evaluating R statistics
#analysisTypes.append('BinValAxisPDF')
#analysisTypes.append('BinValAxisStatsComposite')
#analysisTypes.append('GrossValues')

