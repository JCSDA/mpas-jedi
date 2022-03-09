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

    parser.add_argument('date', default=os.getcwd().split('/')[-3], type=str, nargs = '?',
                        help='Valid date (YYYYMMDDHH)')
    parser.add_argument('-r', '--referenceState', default = mu.GFSANA_DIR+'/x1.40962.init',
                        type = str,
                        help='Path/prefix of reference analysis state')
    parser.add_argument('-m', '--mpasState', default = '../restart', type = str,
                        help='Path/prefix of arbitrary MPAS state')
    parser.add_argument('-n', '--nprocs', default = 1, type = int,
                        help='Number of tasks/processors for multiprocessing')

processor = DiagnoseObsStatisticsArgs()

processor.processArgs()

args = processor.args
