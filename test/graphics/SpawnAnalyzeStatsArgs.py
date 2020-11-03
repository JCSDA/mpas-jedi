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
    parser.add_argument('-a', '--account',
      help='JobScript account number')
    parser.add_argument('-q', '--queue',
      help='JobScript submission queue')
    parser.add_argument('-m', '--memory',
      help='JobScript requested memory per node')
    parser.add_argument('-s', '--scriptdir',
      help='Location of scripts')

processor = SpawnAnalyzeStatsArgs()

processor.processArgs()

args = processor.args
