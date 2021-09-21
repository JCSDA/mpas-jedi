#!/usr/bin/env python3

from SpawnAnalyzeStatsArgs import args
import AnalyzeStats as mainScript
from Analyses import anWorkingDir
import argparse
import config as conf
from predefined_configs import outerIter
from copy import deepcopy
import JobScript as js
import os
from pathlib import Path
import re
import textwrap

jobenv = 'csh'
jobbody = ['''
date

#
# set environment:
# ================
module load python/3.7.5
source /glade/u/apps/ch/opt/usr/bin/npl/ncar_pylib.csh

setenv pySourceDir PYSOURCE

set mainScript = '''+mainScript.__name__+'''

#
# link dependencies:
# ==================
set pyDepends = ( \\
  ${mainScript} \\''']
for dep in mainScript.depends_on:
    jobbody += ['  '+dep+' \\']
jobbody += [''')

foreach pySource ($pyDepends)
  ln -sf ${pySourceDir}/${pySource}.py ./
end

#
# make plots:
# ===========
python ${mainScript}.py -n NPWORK -r NPREAD -d DIAGSPACE -app JEDIAPP -nout NOUTER >& an.log
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

    jobConf['env'] = jobenv

    ## submit a job for each selected DiagSpace
    for DiagSpace, dsConf in DiagSpaceConfig.items():
        myJobConf = deepcopy(jobConf)

        anGroup = dsConf['anGrp']
        grpAtt = conf.anGroupConfig[anGroup]
        npwork = grpAtt['npwork']
        npread = grpAtt['npread']
        myJobConf['nppernode'] = max(npwork, npread)
        myJobConf['walltime'] = grpAtt['analyze_walltime']

        myJobConf['name'] = 'AnStat_'+DiagSpace
        myJobConf['path'] = './'+anWorkingDir(DiagSpace)

        myJobConf['script'] = []
        substitutions = {
            'DIAGSPACE': DiagSpace,
#            'ANALYSISGROUP': anGroup,
            'PYSOURCE': str(scriptDir),
            'NPWORK': str(npwork),
            'NPREAD': str(npread),
            'JEDIAPP': args.jediAppName,
            'NOUTER': str(args.nOuterIter),
        }
        for line in jobbody:
            newline = line
            for key, val in substitutions.items():
                newline = re.sub(key, val, newline)
            myJobConf['script'] += [newline+'\n']

        job = js.JobScriptFactory(myJobConf)
        job.create()
        job.submit()

if __name__ == '__main__': main()

