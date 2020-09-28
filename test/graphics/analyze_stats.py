#!/usr/bin/env python3

import AnalyzeStatistics as AS
import analyze_config as anconf
import argparse
import config as conf
from copy import deepcopy
import diag_utils as du
import logging
import logsetup
import multiprocessing as mp
import os
import StatisticsDatabase as sdb
import textwrap

_logger = logging.getLogger(__name__)

depends_on = [
  'analyze_config',
  'AnalyzeStatistics',
  'basic_plot_functions',
  'binning_configs',
  'binning_params',
  'binning_utils',
  'config',
  'diag_utils',
  'plot_utils',
  'stat_utils',
  'StatisticsDatabase',
  'var_utils',
]

def main():
    '''
    Main function that sequentially
    () processes command-line agruments
    () collates those with static configuration modules (config and analysis_config)
    () loops over selected DiagSpaces, and for each
       - initializes StatisticsDatabase object (internal multiprocessing)
       - analyzes the statistics for all selected anconf.analysisTypes (internal multiprocessing)
    See analysis_config for more information
    '''
    _logger.info('Starting main()')

    DiagSpaceConfig = deepcopy(conf.DiagSpaceConfig)
    for key in sorted(DiagSpaceConfig):
        if not DiagSpaceConfig[key]['process']: del DiagSpaceConfig[key]

    optionalDS = [key for key in sorted(DiagSpaceConfig)]

    # Parse command line
    ap = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("-n", "--npan",
                    help="Number of tasks/processors for analyses")
    ap.add_argument("-r", "--npread",
                    help="Number of tasks/processors for reading StatsDB objects, defaults to npan")
    ap.add_argument("-d", "--diagSpace",
                    help=textwrap.dedent(
                       '''
                       Select DiagSpaces with non-conservative matching
                          e.g., amsua selects amsua_metop-a, amsua_n19, etc.
                          default behavior is to select all DiagSpaces in config
                          Current options:
                       ''')+str(optionalDS))
    ap.add_argument("-g", "--anGrp",
                    help="Select a group of DiagSpaces (overridden by --diagSpace option)")
    ap.add_argument("-a", "--analysisType",
                    help="Select a single analysisType")

    MyArgs = ap.parse_args()

    ## process processor selection
    #  note: these scripts only work on a single node
    if MyArgs.npan:
        npan = min(int(MyArgs.npan), mp.cpu_count())
    else:
        npan = 1
    if MyArgs.npread:
        npread = min(int(MyArgs.npread), mp.cpu_count())
    else:
        npread = npan

    ## process DiagSpace command-line selections
    anGrp = None
    selectDiagSpace = None
    if MyArgs.diagSpace:
        selectDiagSpace = MyArgs.diagSpace
    elif MyArgs.anGrp:
        anGrp = MyArgs.anGrp

    ## remove DiagSpaces that are not selected
    for key in sorted(DiagSpaceConfig):
        if ((selectDiagSpace is not None and selectDiagSpace not in key) or
            (anGrp is not None and DiagSpaceConfig[key]['anGrp'] != anGrp)):
            del DiagSpaceConfig[key]

    ## process analysisType command-line selection or use defaults from analyze_config
    if MyArgs.analysisType:
        analysisTypes = [MyArgs.analysisType]
    else:
        analysisTypes = anconf.analysisTypes

    ## loop over selected DiagSpaces
    for DiagSpaceName in sorted(DiagSpaceConfig):
        ## setup DiagSpaceName configuration
        myDBConf = deepcopy(anconf.dbConf)
        myDBConf['DiagSpaceName'] = DiagSpaceName

        availableDiagNames = []
        for diag in DiagSpaceConfig[DiagSpaceName]['diagNames']:
            if du.availableDiagnostics[diag].get('analyze',True):
                availableDiagNames.append(diag)

        myDBConf['diagnosticConfigs'] = du.diagnosticConfigs(
            availableDiagNames, DiagSpaceName,
            analysisStatistics = anconf.analysisStatistics)

        # Construct statistical database for each DiagSpace
        db = sdb.StatsDB(myDBConf)

        if db.available:
            _logger.info('')
            _logger.info('Analyzing StatsDB for '+DiagSpaceName)

            ## Initialize the database
            db.read(npread)

            analyses = AS.Analyses(db, analysisTypes, npan)

            analyses.analyze()

    _logger.info('Finished main() successfully')

if __name__ == '__main__': main()
