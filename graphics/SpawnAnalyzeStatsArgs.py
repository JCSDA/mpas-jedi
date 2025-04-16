#!/usr/bin/env python3

import jediApplicationArgs
from ProcessArgs import ProcessArgs
import textwrap

class SpawnAnalyzeStatsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()
    self.argProcessors += [jediApplicationArgs]

  @staticmethod
  def add_arguments(parser):
    parser.add_argument('-d', '--diagSpaces',
      help=textwrap.dedent(
        '''
        Comma-separated list of DiagSpaces with non-conservative matching
           e.g., amsua selects amsua_metop-a, amsua_n19, etc.
           default behavior is to select all DiagSpaces in config
        '''))
    parser.add_argument('-a', '--account', type = str,
      help='JobScript account number')
    parser.add_argument('-q', '--queue', type = str,
      help='JobScript submission queue, one of \n'
           '    "casper" job will run on a shared casper htc node \n'
           '    "gpudev" job will run on a shared casper gpu node \n'
           '    "develop" job will run on a shared derecho node \n'
           '    "main" job will run on a dedicated derecho node (not recommended) \n'
           'Jobs on shared nodes only get billed for the cpus they use, \n'
           '  e.g. -n 1 will only use one cpu and be billed for 1 * wall-clock hours. \n'
           'Jobs with lower cpu and memory requirements often get run sooner.')
    parser.add_argument('-m', '--memory', type = int,
      help='JobScript requested memory per node in GB')
    parser.add_argument('-n', '--nppernode', type = int,
      help='JobScript requested number of processors per node')
    parser.add_argument('-s', '--scriptdir', type = str,
      help='Location of scripts')
    parser.add_argument("-c", "--controlExperiment",
      help="Base experiment for making comparison graphs")
    parser.add_argument("-e", "--experiments",
      help="Comma separated list of <short:long_name> experiments to graph")
    parser.add_argument("-p", "--verifySpace",
      help="verificationSpace (model or obs)")
    parser.add_argument("-t", "--verifyType",
      help="verificationType (forecast or omb/oma)")
    parser.add_argument("-f", "--firstCycle",
      help="first Cycle date/time, e.g. 20180414T18")
    parser.add_argument("-l", "--lastCycle",
      help="last Cycle date/time, e.g. 20180415T06")
    parser.add_argument("-w", "--waitForSync", action="store_true",
      help="create a sync job which won't run until all graph jobs complete")

processor = SpawnAnalyzeStatsArgs()

processor.processArgs()

args = processor.args
