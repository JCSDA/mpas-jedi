#!/usr/bin/env python3

from collections import OrderedDict
import datetime as dt
import os
import var_utils as vu

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

    ./SpawnAnalyzeStats.py

+ Specify that statistics files come from a JEDI hofx application

    ./SpawnAnalyzeStats.py -app hofx

+ Spawn one job for each DiagSpace that contains "amsua"

    ./SpawnAnalyzeStats.py -d amsua

+ Spawn one job for the MPAS model DiagSpace

    ./SpawnAnalyzeStats.py -d mpas

+ Spawn one job for each DiagSpace in the typical MPAS-Workflow variational application

    ./SpawnAnalyzeStats.py -nout 2 -d amsua_,sonde,airc,sfc,gnssrobndropp1d,satwind

+ Spawn one job for each DiagSpace in the typical MPAS-Workflow hofx application

    ./SpawnAnalyzeStats.py -app hofx -d mhs,amsua,abi_,ahi_,sonde,airc,sfc,gnssrobndropp1d,gnssrorefncep,satwind

+ Choose a unique job account number

    ./SpawnAnalyzeStats.py -a NMMM0043

+ Get info about more options

    ./SpawnAnalyzeStats.py --help

'''

'''
readjust configuration given the provided run time values for:
  control the control experiment for comparison graphs
  expArray an array of colon delimited short and long array names,
    e.g. ["exp1:exp1_long_name", "exp2:exp2_long_name",...]
  verifySpace the verification space (model or obs)
  verifyType the verification type (forecast or omb/oma)
  firstCycle a datetime object of the first cycle to graph
  lastCycle a datetime object of the last cycle to graph
'''
def adjust_experiments(control, expArray, verifySpace, verifyType, firstCycle, lastCycle):
  print('control:', control)
  print('experiments:', expArray)
  print(verifySpace, verifyType, firstCycle, lastCycle)

  verificationChanged = False
  exps = None
  if control:
    dbConf['cntrlExpName'] = control

  if verifySpace:
    VerificationSpace = verifySpace
    verificationChanged = True
  if verifyType:
    VerificationType = verifyType
    verificationChanged = True
  if firstCycle:
    dbConf['firstCycleDTime'] = firstCycle
  if lastCycle:
    dbConf['lastCycleDTime'] = lastCycle

  ## if the verification space/type changed, reconfigure
  if verificationChanged:
    if zeroDurationForecast:
      if VerificationType == 'omb/oma' and VerificationSpace == 'obs':
        obsAppIdentifier = 'da'
        deterministicVerifyDir = OMBOMAVerification
        # ensemble verification for omb/oma not currently supported by MPAS-Workflow

      if VerificationType == 'forecast':
        # single forecast duration omf (fcTDeltaFirst==0 and fcTDeltaLast == 0)
        deterministicVerifyDir = ShortRangeFCVerification
        ensembleVerifyDir = ShortRangeFCVerification
        obsAppIdentifier = 'hofx'
    else:
      # multiple extended forecast omf (0 < fcTDeltaLast <= extended forecast length)
      # override VerificationType; only 'forecast' is available
      VerificationType = 'forecast'
      deterministicVerifyDir = ExtendedFCVerification
      ensembleVerifyDir = ExtendedFCVerification
      obsAppIdentifier = 'hofx'

    if workflowType == 'MPAS-Workflow':
      commonStatsFileSubDir = statsFileSubDirBase+'/'+VerificationSpace

    # set up the experiments array
    if expArray is not None:
      experiments = OrderedDict()
      for exp in expArray:
        shortLong = exp.split(':')
        experiments[shortLong[0]] = shortLong[1] + deterministicVerifyDir

    ## expNames is a list of experiment names used for database lookups and figure labels, e.g., legend
    #  entries.  expNames must have the same length as expLongNames.  Make these brief names concise
    #  and exclude spaces. It is recommended to only use characters from [A-Z], [a-z], [0-9], and
    #  ['-', '+', '=', '_', ',', ';', '.'], although others are allowed.
    dbConf['expNames'] = []

    ## expLongNames is a list of directories within expDirectory. Each list member is associated with a
    #  unique experiment.  Add to the list to analyze for experiments simultaneously.
    dbConf['expLongNames'] = []

    for key, value in experiments.items():
      dbConf['expNames'].append(key)
      dbConf['expLongNames'].append(value)

    if VerificationSpace == 'model':
      dbConf['appIdentifiers'] = ['']*len(dbConf['expNames'])
    else:
      dbConf['appIdentifiers'] = [obsAppIdentifier]*len(dbConf['expNames'])

    assert dbConf['cntrlExpName'] in dbConf['expNames'], (
      'cntrlExpName must be one of the available expNames')

    ## statsFileSubDirs is the final subdirectory within the date directory(ies)
    #  that contains the statstics files for constructing the StatsDB object
    dbConf['statsFileSubDirs'] = [commonStatsFileSubDir]*len(dbConf['expNames'])

    if VerificationSpace == 'model':
      dbConf['appIdentifiers'] = ['']*len(dbConf['expNames'])
    else:
      dbConf['appIdentifiers'] = [obsAppIdentifier]*len(dbConf['expNames'])

    commonFCDirFormat = "%hhr"
    dbConf['fcDirFormats'] = [commonFCDirFormat]*len(dbConf['expNames'])


## ================================================================================================
## ================================================================================================
## dbConf: configureation for all StatisticsDatabase.StatsDB objects
dbConf = {}

## ================================================================================================
## -------------------------------
## Cycle times and forecast length
## -------------------------------
# first and last valid Cycle date-times (DTime) (e.g., analysis time, forecast initial time)
dbConf['firstCycleDTime'] = dt.datetime(2018,4,15,0,0,0)
dbConf['lastCycleDTime'] = dt.datetime(2018,5,14,18,0,0)

# time increment (TimeInc) between valid Cycle (cy) date-times
dbConf['cyTimeInc'] = dt.timedelta(hours=6)

# first and last forecast (fc) durations (TDelta)
dbConf['fcTDeltaFirst'] = dt.timedelta(days=0)
dbConf['fcTDeltaLast'] = dt.timedelta(days=0,hours=0)

if (dbConf['fcTDeltaFirst'].total_seconds() == 0 and
    dbConf['fcTDeltaLast'].total_seconds() == 0):
  zeroDurationForecast = True
  #note: a different directory structure is assumed for zeroDurationForecast
else:
  zeroDurationForecast = False

# time increment (TimeInc) between forecast (fc) durations
# E.g., interval between each execution of verification on forecast output files
dbConf['fcTimeInc'] = dt.timedelta(hours=12)


## ================================================================================================
## -------------------------------------------------------
## Verification data sub-directory configuration
## -------------------------------------------------------

# VerificationType and VerificationSpace are super-settings that automatically control the directory
# structure where the statistics files are located.  The directory names are derived from
# MPAS-Workflow.  If the user has a different directory structure, it is recommended to override the
# values for
# + OMBOMAVerification, ShortRangeFCVerification, ExtendedFCVerification, and
#   deterministicFCMemberDir
# OR
# + deterministicVerifyDir, ensembleVerifyDir, obsAppIdentifier, and commonStatsFileSubDir

## VerificationType
# OPTIONS: 'omb/oma', 'forecast'
# 'omb/oma' - calculated from a da application, only available when
#             VerificationSpace=='obs'
# 'forecast' - single- or multi-duration forecasts either in observation or model space
VerificationType = 'forecast'

## VerificationSpace
# OPTIONS: 'obs', 'model'
# 'obs' - observation space
# 'model' - compare to analyses in model space, only available when VerificationType=='forecast'
VerificationSpace = 'obs'

## workflowType
# allows easy selection of directory structures from between workflow designs
# OPTIONS: 'MPAS-Workflow', 'liuz,jban'
workflowType = 'MPAS-Workflow'

## -----------------------------------------------------------
## Automated sub-directory selection from above configurations
## -----------------------------------------------------------

## OMBOMAVerification, ShortRangeFCVerification, ExtendedFCVerification, and
#  deterministicFCMemberDir are used as presets for deterministicVerifyDir and ensembleVerifyDir

# MPAS-Workflow settings
if workflowType == 'MPAS-Workflow':
  deterministicMemberDir = '/mean'
  OMBOMAVerification = '/Verification/da'+deterministicMemberDir
  ShortRangeFCVerification = '/Verification/bg'+deterministicMemberDir
  ExtendedFCVerification = '/Verification/fc'+deterministicMemberDir

# liuz, jban settings
if workflowType == 'liuz,jban':
  OMBOMAVerification = '/DADIAG'
  ShortRangeFCVerification = '/FC1DIAG'
  ExtendedFCVerification = '/FC2DIAG'

## deterministicVerifyDir, ensembleVerifyDir, and obsAppIdentifier are assumed to be
# uniform across all experiments

## deterministicVerifyDir, set automatically from above settings
#  subdirectory within a deterministic experiment directory where date subdirectories are located
deterministicVerifyDir = 'NotSupported'

## ensembleVerifyDir, set automatically from above settings
#  subdirectory within an ensemble experiment directory where date subdirectories are located
ensembleVerifyDir = 'NotSupported'

if zeroDurationForecast:
  if VerificationType == 'omb/oma' and VerificationSpace == 'obs':
    obsAppIdentifier = 'da'
    deterministicVerifyDir = OMBOMAVerification
    # ensemble verification for omb/oma not currently supported by MPAS-Workflow

  if VerificationType == 'forecast':
    # single forecast duration omf (fcTDeltaFirst==0 and fcTDeltaLast == 0)
    deterministicVerifyDir = ShortRangeFCVerification
    ensembleVerifyDir = ShortRangeFCVerification
    obsAppIdentifier = 'hofx'

else:
  # multiple extended forecast omf (0 < fcTDeltaLast <= extended forecast length)
  # override VerificationType; only 'forecast' is available
  VerificationType = 'forecast'
  deterministicVerifyDir = ExtendedFCVerification
  ensembleVerifyDir = ExtendedFCVerification
  obsAppIdentifier = 'hofx'


## ================================================================================================
## ----------------------------------------
## Experiment-specific sub-directory naming
## ----------------------------------------
# The overall directory structure is defined in StatisticsDatabase.StatsDB.  For each experiment key
# in the `experiments` dictionary below, its statistics files are located in
# if zeroDurationForecast:
#   dbConf['expDirectory']+'/'+experiments[key]+'/YYYYMMDDHH/'+statsFileSubDir
# else:
#   dbConf['expDirectory']+'/'+experiments[key]+'/YYYYMMDDHH/FCDuration/'+statsFileSubDir

# where

# YYYYMMDDHH is the valid date of the forecast initial state, DA background, or DA analysis
# FCDuration is the forecast length, formatted used commonFCDirFormat below

# statsFileSubDir is a user-optional subdirectory where the statistics are stored, and
# by defult, it is set to commonStatsFileSubDir below

## expDirectory is the top-level directory that contains data
#  from all experiments.  The environment variable, EXP_DIR, can be used
#  to specify its value externally if desired.
user = os.environ['USER']
dbConf['expDirectory'] = os.getenv('EXP_DIR','/glade/derecho/scratch/'+user+'/pandac')

## hasFCLenDir whether directory structure includes forecast length
#  overridden to True within StatisticsDatabase.StatsDB when fcTDeltaLast > fcTDeltaFirst
dbConf['hasFCLenDir'] = False

## cntrlExpName is the experiments key of the control experiment, which is used for DiffCI analyses
dbConf['cntrlExpName'] = 'clrama'

## experiments - dictionary with key, value pairs as follows
#  + the key is a short name for the experiment (see expNames below)
#  + the value is the directory where the verification statistics files are located
#  + if using MPAS-Workflow, users only need to add one new `experiments` entry per experiment and
#    select their desired VerificationType and VerificationSpace above

experiments = OrderedDict()

experiments['clrama'] = \
    'guerrett_3dhybrid-60-60-iter_gnssrorefncep_O30kmI60km_ensB-SE80+RTPP70_VarBC_RefNCEP_2ndDoaDob' + \
    deterministicVerifyDir


## Additional examples for adding experiments
## ------------------------------------------

## -------------------------------------------------------------------------
## (1) Single-state (deterministic background) experiment
#  This is an instance of the 2018 120km 3denvar clrama experiment,
#  which assimilates conventional obs and clear-sky AMSUA brightness
#  temperature
#  date simulated: JUNE 2021 (update as needed)

#experiments['clrama'] = \
#  'guerrett_3denvar_OIE120km_ioda-v2'+deterministicVerifyDir

## -------------------------------------------------------------------------
## (2) Add another single-state experiment
#  This will be compared to the 'clrama' experiment above

#experiments['variant-1'] = \
#  'guerrett_3denvar_OIE120km_VARIANT1'+deterministicVerifyDir

## -------------------------------------------------------------------------
## (3) Add an ensemble-state experiment
#  omf will be y - h(Mean(xf)), where Mean(xf) is the mean of an ensemble of
#  6-hr forecasts from NCEP's GEFS analyses.
#  + the same procedure can be used for MPAS-Workflow EDA experiment output
#  + add the 'mean' ensemble member

#experiments['gefs-cold'] = \
#  'guerrett_eda_3denvar_NMEM20_OIE120km_GEFSVerify'+ensembleVerifyDir

#  + add the 5th ensemble member

#experiments['gefs-cold-5'] = \
#  'guerrett_eda_3denvar_NMEM20_OIE120km_GEFSVerify/Verification/bg/mem005'

## -------------------------------------------------------------------------
## (4) Make spaghetti plots for a 20-member ensemble of cold-start
#  6-hr forecasts from GEFS analyses
#  + the default nSpaghetti arguments to plotColor and plotLineStyle in
#    plot_utiles need to be set equal to nEnsDAMembers below (the number of
#    ensemble members).  For this to work correctly, this ensemble of
#    experiments must be first, and all other experiments should be placed
#    after this one

#nEnsDAMembers = 20

#  + treat each ensemble member as a different experiment

#thisExpDir = 'guerrett_eda_3denvar_NMEM20_OIE120km_GEFSVerify'
#for mem in list(range(1,nEnsDAMembers+1)):
#  member = '{:03d}'.format(mem)
#  experiments['gefs-cold-'+str(mem)] = thisExpDir+'/Verification/bg/mem'+member

#  + the same procedure can be used for MPAS-Workflow EDA experiment output
#  + this kind of figure is very expensive to generate for large nEnsDAMembers

## ================================================================================================
## The settings below are automatically populated from the above configurations

## expNames is a list of experiment names used for database lookups and figure labels, e.g., legend
#  entries.  expNames must have the same length as expLongNames.  Make these brief names concise
#  and exclude spaces. It is recommended to only use characters from [A-Z], [a-z], [0-9], and
#  ['-', '+', '=', '_', ',', ';', '.'], although others are allowed.
dbConf['expNames'] = []

## expLongNames is a list of directories within expDirectory. Each list member is associated with a
#  unique experiment.  Add to the list to analyze for experiments simultaneously.
dbConf['expLongNames'] = []

for key, value in experiments.items():
  dbConf['expNames'].append(key)
  dbConf['expLongNames'].append(value)

assert dbConf['cntrlExpName'] in dbConf['expNames'], (
  'cntrlExpName must be one of the available expNames')


## commonStatsFileSubDir, set automatically from above settings
# subdirectory within each date subdirectory where statistics files are located, assumed to be
# uniform across all experiments
statsFileSubDirBase = 'diagnostic_stats'

# MPAS-Workflow settings
if workflowType == 'MPAS-Workflow':
  commonStatsFileSubDir = statsFileSubDirBase+'/'+VerificationSpace

# liuz, jban settings
if workflowType == 'liuz,jban':
  commonStatsFileSubDir = statsFileSubDirBase

## statsFileSubDirs is the final subdirectory within the date directory(ies)
#  that contains the statstics files for constructing the StatsDB object
dbConf['statsFileSubDirs'] = [commonStatsFileSubDir]*len(dbConf['expNames'])

## appIdentifiers is a list of secondary labels that are within the database
#  intermediate file names.  Typical values in MPAS-Workflow applications are
#  'hofx' for OMF, 'da' for OMB/OMA, and '' (empty) for model-space
#  verification.  It is optional for the user to use whatever string they choose
#  as part of the file name to distinguish the information contained within from
#  other files, e.g., an experiment characteristic.  appIdentifiers is only important
#  for file naming and is not used in the analyses
if VerificationSpace == 'model':
  dbConf['appIdentifiers'] = ['']*len(dbConf['expNames'])
else:
  dbConf['appIdentifiers'] = [obsAppIdentifier]*len(dbConf['expNames'])

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


## ================================================================================================
#Select the statistics to analyze
# options: see stat_utils.allFileStats
# 'selectedStatistics' for individual diagnostics (see diag_utils) will override this setting
selectedStatistics = ['Count', 'Mean', 'RMS', 'STD']
#selectedStatistics = ['Skew', 'ExcessKurtosis']

#Select diagnostic groupings (if applicable)
# + each entry in the dictionary should be of the format:
#   diagnosticGroup: diagnosticNames
# + diagnosticNames is a list of diagnostics, such as
#   [diagnosticName1, diagnosticName2, etc...]
# + diagnosticGroup (the dict key) is up to the user, but it is best to have it match
#   the strings used for axis labeling in Analyses
# + the default behavior is to plot each diagnostic on an independent set of axes,
#   which will still be done for any analysis type that does not use
#   diagnosticGroupings or has maxDiagnosticsPerAnalysis < len(diagnosticNames),
#   e.g., CYandBinValAxes2D, FCandBinValAxes2D, and BinValAxes2D
diagnosticGroupings = {}
# OMM
diagnosticGroupings['omm'] = ['omb', 'oma']
diagnosticGroupings['rltv_omm'] = ['rltv_omb', 'rltv_oma']
diagnosticGroupings['omm_nobc'] = ['omb_nobc', 'oma_nobc']

# ObsError
diagnosticGroupings['ObsError'] = ['sigmaob', 'doadob', 'ApproxDoaDob']
diagnosticGroupings['RelativeObsError'] = ['rltv_sigmaob', 'rltv_doadob', 'ApproxRelativeDoaDob']
diagnosticGroupings['sigmao'] = ['sigmaof', 'ideal-sigmaof']
#diagnosticGroupings['rltv_sigmao'] = ['rltv_omf', 'ideal-rltv_sigmaof']
diagnosticGroupings['ErrorRatios'] = ['OENIb', 'OENIa', 'InnovationRatio']

# R, HBH^T, obs consistency
diagnosticGroupings['ddT'] = ['doadob', 'dabdob']
diagnosticGroupings['Relative-ddT'] = ['rltv_doadob', 'rltv_dabdob']
#SpreadDiagnostics = [vu.EddT, vu.HBHT, vu.R]
SpreadDiagnostics = [vu.EddT, vu.HBHTplusR]
diagnosticGroupings['dy'] = [vu.DiagnosticVars[dy] for dy in SpreadDiagnostics]
#diagnosticGroupings['dy'] = ['omf', 'sigmaof', 'sigmaf']
diagnosticGroupings['CRy'] = ['CRyb', 'CRya']

# model spread and consistency
#diagnosticGroupings['sigmax'] = ['sigmaxb', 'sigmaxa', 'sigmainf']
#diagnosticGroupings['dx'] = ['mmgfsan', 'sigmaxb']
#diagnosticGroupings['CRx'] = ['CRxb', 'CRxa']
#diagnosticGroupings['SRx'] = ['SRx-eda', 'SRx-rtpp']


########################################################################
## Configure the analysisTypes to apply to the statistics
#  - below are recommendations for single/multiple forecast lengths
#  - analysisTypes can be mixed and matched as desired,
#    however some of them require nCY, nFC, or nExp > 1
#  - see the individual classes for more details (Analyses.py)
analysisTypes = []
if zeroDurationForecast:
  print('Generating CY-type figures')
  ## gross error analysisTypes for single forecast length
  ## -------------------------------------------------------
  ## recommended
  analysisTypes.append('CYAxisExpLines')
  analysisTypes.append('CYandBinValAxes2D')
  analysisTypes.append('BinValAxes2D')
  if len(dbConf['expNames']) > 1: analysisTypes.append('BinValAxisProfileDiffCI')

  ## useful for prescribing/evaluating R statistics
  #analysisTypes.append('BinValAxisPDF')
  #analysisTypes.append('BinValAxisPDFMultiExp')
  analysisTypes.append('BinValAxisStatsComposite')
  #analysisTypes.append('GrossValues')
  analysisTypes.append('BinValAxisProfile')

  ## potentially useful
  #analysisTypes.append('CYAxisBinValLines')

else:
  print('Generating FC-type figures')
  ## gross error analysisTypes for multiple forecast lengths
  ## -------------------------------------------------------
  ## recommended
  if len(dbConf['expNames']) > 1: analysisTypes.append('FCScoreCard')
  if len(dbConf['expNames']) > 1: analysisTypes.append('BinValAxisProfileDiffCI')

  ## potentially useful
  #analysisTypes.append('FCandBinValAxes2D')
  #analysisTypes.append('FCAxisExpLines')
  #if len(dbConf['expNames']) > 1: analysisTypes.append('FCAxisExpLinesDiffCI')
  #analysisTypes.append('CYandBinValAxes2D')
  #analysisTypes.append('CYAxisExpLines')
  #analysisTypes.append('CYAxisFCLines')
  #analysisTypes.append('BinValAxisProfile') 
  #analysisTypes.append('BinValAxes2D')
  #analysisTypes.append('BinValAxisStatsComposite')


