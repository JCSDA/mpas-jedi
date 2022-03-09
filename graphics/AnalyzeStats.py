#!/usr/bin/env python3

from AnalyzeStatsArgs import args

from Analyses import Analyses
import analyze_config as anconf
import config as conf
from copy import deepcopy
import diag_utils as du
import logging
import logsetup
import multiprocessing as mp
import os
import StatisticsDatabase as sdb

_logger = logging.getLogger(__name__)

depends_on = [
  'AnalyzeStatsArgs',
  'analyze_config',
  'Analyses',
  'basic_plot_functions',
  'predefined_configs',
  'binning_params',
  'binning_utils',
  'config',
  'diag_utils',
  'fit2D',
  'plot_utils',
  'stat_utils',
  'StatisticsDatabase',
  'var_utils',
]

def main():
    '''
    Main function that sequentially
    () collates command-line arguments from AnalyzeStatsArgs with static
       configuration modules (config and analysis_config)
    () loops over selected DiagSpaces, and for each
       - initializes StatisticsDatabase object (internal multiprocessing)
       - analyzes the statistics for all selected anconf.analysisTypes (internal multiprocessing)
    See analysis_config for more information
    '''
    _logger.info('Starting main()')

    DiagSpaceConfig = deepcopy(conf.DiagSpaceConfig)
    for key in sorted(DiagSpaceConfig):
        if not DiagSpaceConfig[key]['process']: del DiagSpaceConfig[key]

    ## process processor selection
    #  note: these scripts only work on a single node
    if args.npan:
        npan = min(int(args.npan), mp.cpu_count())
    else:
        npan = 1
    if args.npread:
        npread = min(int(args.npread), mp.cpu_count())
    else:
        npread = npan

    ## process DiagSpace command-line selections
    anGrp = None
    selectDiagSpace = None
    if args.diagSpace:
        selectDiagSpace = args.diagSpace
    elif args.anGrp:
        anGrp = args.anGrp

    ## remove DiagSpaces that are not selected
    for key in sorted(DiagSpaceConfig):
        if ((selectDiagSpace is not None and selectDiagSpace not in key) or
            (anGrp is not None and DiagSpaceConfig[key]['anGrp'] != anGrp)):
            del DiagSpaceConfig[key]

    ## process analysisType command-line selection or use defaults from analyze_config
    if args.analysisType:
        analysisTypes = [args.analysisType]
    else:
        analysisTypes = anconf.analysisTypes

    ## loop over selected DiagSpaces
    for DiagSpaceName in sorted(DiagSpaceConfig):
        ## setup DiagSpaceName configuration
        myDBConf = deepcopy(anconf.dbConf)
        myDBConf['DiagSpaceName'] = DiagSpaceName

        myDBConf['diagnosticConfigs'] = du.diagnosticConfigs(
            DiagSpaceConfig[DiagSpaceName]['diagNames'], DiagSpaceName,
            selectedStatistics = anconf.selectedStatistics)
        for diag in list(myDBConf['diagnosticConfigs'].keys()):
            if not myDBConf['diagnosticConfigs'][diag]['analyze']:
                del myDBConf['diagnosticConfigs'][diag]

        # Construct statistical database for each DiagSpace
        db = sdb.StatsDB(myDBConf)

        if db.available:
            _logger.info('')
            _logger.info('Analyzing StatsDB for '+DiagSpaceName)

            ## Initialize the database
            db.read(npread)

            analyses = Analyses(db, analysisTypes, anconf.diagnosticGroupings, npan)

            analyses.analyze()

    _logger.info('Finished main() successfully')

if __name__ == '__main__': main()
