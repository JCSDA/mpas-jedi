#!/usr/bin/env python3

import argparse
import modelsp_utils as mu
import os
from ProcessArgs import ProcessArgs

class DiagnoseObsStatisticsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()

  @staticmethod
  def add_arguments(parser):

    parser.add_argument('date', default=os.getcwd().split('/')[-3],
                        type=str, nargs = '?',
                        help='Valid date (YYYYMMDDHH)')
    parser.add_argument('-r', '--referenceState', default = mu.GFSANA_DIR+'/x1.40962.init',
                        type = str,
                        help='Path/prefix of reference analysis state')
    parser.add_argument('-p', '--meanState', default = '../restart',
                        type = str,
                        help='Path/prefix of arbitrary MPAS state')
    parser.add_argument('-rd', '--referenceDiagnostics', default = 'None',
                        type = str,
                        help='Path/prefix of reference analysis state')
    parser.add_argument('-pd', '--meanDiagnostics', default = '../diag',
                        type = str,
                        help='Path/prefix of arbitrary MPAS state')

    parser.add_argument('-n', '--nprocs', default = 1,
                        type = int,
                        help='Number of tasks/processors for multiprocessing')
    parser.add_argument("-m", "--nMembers", default = 0,
                        type = int,
                        help="number of ensemble members; must be >1 to produce ensemble diagnostic stats")
    parser.add_argument("-b", "--bgEnsemblePath", default = '../../../../../../CyclingDA/YYYYMMDDHH/bg/mem{:03d}/bg',
                        type = str,
                        help="Path/prefix of background ensemble member states; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")
    parser.add_argument("-a", "--anEnsemblePath", default = '../../../../../../CyclingInflation/RTPP/YYYYMMDDHH/an/mem{:03d}/an',
                        type = str,
                        help="Path/prefix of analysis ensemble member states; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")
    parser.add_argument("-i", "--inflatedEnsemblePath", default = '../../../../../../CyclingDA/YYYYMMDDHH/an/mem{:03d}/an',
                        type = str,
                        help="Path/prefix of inflated ensemble member states; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")

#    parser.add_argument("-f", "--fcEnsemblePath", default = '../../../../../../CyclingFC/YYYYMMDDHH/mem{:03d}/mpasout',
#                        type = str,
#                        help="Path/prefix of forecast ensemble member states; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")

processor = DiagnoseObsStatisticsArgs()

processor.processArgs()

args = processor.args
