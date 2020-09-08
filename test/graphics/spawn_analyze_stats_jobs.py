#!/usr/bin/env python3

#from analyze_config import analysisTypes
import analyze_stats as mainScript
from AnalyzeStatistics import anWorkingDir
import argparse
import config as conf
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
python ${mainScript}.py -n NPWORK -r NPREAD -d DIAGSPACE >& an.log
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
    () collates command-line agruments and static configuration modules
    () loops over selected DiagSpaces, and for each
       - initializes StatisticsDatabase object (internal multiprocessing)
       - analyzes the statistics for all selectd analysisTypes (internal multiprocessing)
    '''

    DiagSpaceConfig = deepcopy(conf.DiagSpaceConfig)
    for key in sorted(DiagSpaceConfig):
        if not DiagSpaceConfig[key]['process']: del DiagSpaceConfig[key]

    optionalDS = [key for key in sorted(DiagSpaceConfig)]

    # Parse command line
    ap = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument('-d', '--diagSpaces',
                    help=textwrap.dedent(
                       '''
                       Comma-separated list of DiagSpaces with non-conservative matching
                          e.g., amsua selects amsua_metop-a, amsua_n19, etc.
                          default behavior is to select all DiagSpaces in config
                          Current options:
                       ''')+str(optionalDS))
    ap.add_argument('-a', '--account',
                    help='JobScript account number')
    ap.add_argument('-q', '--queue',
                    help='JobScript submission queue')
    ap.add_argument('-m', '--memory',
                    help='JobScript requested memory per node')
    ap.add_argument('-s', '--scriptdir',
                    help='Location of scripts')

    MyArgs = ap.parse_args()

    ## process DiagSpace command-line selection
    selectDiagSpaces = None
    if MyArgs.diagSpaces:
        selectDiagSpaces = str(MyArgs.diagSpaces).split(',')

    ## remove DiagSpaces that are not selected
    for key in sorted(DiagSpaceConfig):
        if selectDiagSpaces is not None:
            match = False
            for space in selectDiagSpaces:
                if space in key: match = True
            if not match: del DiagSpaceConfig[key]

    ## set python source directory
    if MyArgs.scriptdir:
        ## from command-line argument, if provided
        scriptDir = Path(MyArgs.scriptdir)
    else:
        ## otherwise, use current directory
        scriptDir = Path(os.getcwd())

    jobConf = {}

    ## get job configuration command-line arguments
    if MyArgs.account: jobConf['account'] = MyArgs.account
    if MyArgs.queue: jobConf['queue'] = MyArgs.queue
    if MyArgs.memory: jobConf['memory'] = MyArgs.memory

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

