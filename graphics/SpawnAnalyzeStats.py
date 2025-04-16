#!/usr/bin/env python3

from SpawnAnalyzeStatsArgs import args

import argparse
import config as conf
from predefined_configs import outerIter
from copy import deepcopy
import JobScript as js
import os
import subprocess
from pathlib import Path
import re
import textwrap

import AnalyzeStats as mainScript
from analyze_config import analysisTypes

from analysis.AnalysisBase import anWorkingDir

jobenv = 'csh'
jobbody = ['''
date

#
# set environment:
# ================
source /etc/profile.d/z00_modules.csh
setenv HDF5_DISABLE_VERSION_CHECK 1
setenv NUMEXPR_MAX_THREADS 1

setenv pySourceDir PYSOURCE

set mainScript = '''+mainScript.__name__+'''
ln -sf ${pySourceDir}/${mainScript}.py ./

#
# make plots:
# ===========
module purge
#module load python
#source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh default
#module load ncarenv/1.3
#module load gnu/10.1.0
#module load ncarcompilers/0.5.0
#module load netcdf/4.8.1
module load conda/latest
conda activate npl-2023a
setenv PYTHONDONTWRITEBYTECODE 1 # avoid __pycache__ creation
module list

set success = 1
set try = 0
while ($success != 0 && $try < 5)
  @ try++
  echo "try=$try"

  python ${mainScript}.py -n NPWORK -r NPREAD -d DIAGSPACE -app JEDIAPP -nout NOUTER -a ANALYSISTYPE EXPARGS >& an.log
  grep 'TypeError: super() takes at least 1 argument (0 given)' an.log
  if ( $status == 0 ) then
    sleep 2
    deactivate
    module purge
    #module load python
    module load conda/latest
    conda activate npl-2023a
    module list
    #source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh default
  else
    set success = 0
  endif
end

grep 'Finished main() successfully' an.log
if ( $status != 0 ) then
  touch ./FAIL
  exit 1
endif

rm ./*.py

date

exit''']

def main():
    '''
    Main function that sequentially
    () collates command-line agruments (SpawnAnalyzeStatsArgs) and static
       configuration module (config)
    () loops over selected DiagSpaces, and for each
       - spawns a job that executes AnalyzeStats on multiple processors
    '''
    DiagSpaceConfig = deepcopy(conf.DiagSpaceConfig)
    for key in sorted(DiagSpaceConfig):
        if not DiagSpaceConfig[key]['process']: del DiagSpaceConfig[key]

    ## process DiagSpace command-line selection
    selectDiagSpaces = None
    if args.diagSpaces:
        selectDiagSpaces = str(args.diagSpaces).split(',')

    ## remove DiagSpaces that are not selected
    for key in sorted(DiagSpaceConfig):
        if selectDiagSpaces is not None:
            match = False
            for space in selectDiagSpaces:
                if space in key: match = True
            if not match: del DiagSpaceConfig[key]

    ## set python source directory
    if args.scriptdir:
        ## from command-line argument, if provided
        scriptDir = Path(args.scriptdir)
    else:
        ## otherwise, use current directory
        scriptDir = Path(os.getcwd())

    jobConf = {}

    ## get job configuration command-line arguments
    if args.account: jobConf['account'] = args.account
    if args.queue: jobConf['queue'] = args.queue
    if args.memory: jobConf['memory'] = args.memory
    if args.nppernode: jobConf['nppernode'] = args.nppernode

    ## get experiment command-line arguments
    expArgs = ""
    if args.controlExperiment: expArgs += " -c " + args.controlExperiment
    if args.experiments: expArgs += " -e " + args.experiments
    if args.verifySpace: expArgs += " -p " + args.verifySpace
    if args.verifyType: expArgs += " -t " + args.verifyType
    if args.firstCycle: expArgs += " -f " + args.firstCycle
    if args.lastCycle: expArgs += " -l " + args.lastCycle

    ## if experiments were provided, adjust the analysis types
    if args.experiments:
      exps = args.experiments.split(',')
      if len(exps) > 1:analysisTypes.append('BinValAxisProfileDiffCI')

    ## check for command line option to create a sync job
    if args.waitForSync:
      sync = True
    else:
      sync = False

    jobConf['env'] = jobenv

    ## submit a job for each selected DiagSpace
    job_procs = []
    for DiagSpace, dsConf in DiagSpaceConfig.items():
      for analysisType in analysisTypes:
        myJobConf = deepcopy(jobConf)

        anGroup = dsConf['anGrp']
        grpAtt = conf.anGroupConfig[anGroup]
        npwork = grpAtt['npwork']
        myJobConf['nnode'] = npwork
        npread = grpAtt['npread']
        if 'nppernode' not in jobConf:
          myJobConf['nppernode'] = max(npwork, npread)
        myJobConf['walltime'] = grpAtt['analyze walltime']

        myJobConf['name'] = 'AnStat_'+DiagSpace+'_'+analysisType
        myJobConf['filename'] = 'AnStat_'+DiagSpace
        myJobConf['path'] = './'+anWorkingDir(DiagSpace, analysisType)

        myJobConf['script'] = []
        substitutions = {
            'DIAGSPACE': DiagSpace,
            'ANALYSISTYPE': analysisType,
            'PYSOURCE': str(scriptDir),
            'NPWORK': str(npwork),
            'NPREAD': str(npread),
            'JEDIAPP': args.jediAppName,
            'NOUTER': str(args.nOuterIter),
            'EXPARGS': expArgs,
        }
        for line in jobbody:
            newline = line
            for key, val in substitutions.items():
                newline = re.sub(key, val, newline)
            myJobConf['script'] += [newline+'\n']

        job = js.JobScriptFactory(myJobConf)
        job.create()
        job_procs.append(job.submit(sync))

    if len(job_procs) > 0:

      # wait for all jobs to complete, order doesn't matter
      for job_proc in job_procs:
        if sync == True:
          print("waiting for subprocess:", job_proc)
        job_proc.wait()

if __name__ == '__main__': main()

