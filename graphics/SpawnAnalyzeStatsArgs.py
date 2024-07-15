#!/usr/bin/env python3

import traceback
import jediApplicationArgs
from ProcessArgs import ProcessArgs
import textwrap

class SpawnAnalyzeStatsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()
    self.argProcessors += [jediApplicationArgs]
    self.processArgs()

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
      help='JobScript submission queue')
    parser.add_argument('-m', '--memory', type = str,
      help='JobScript requested memory per node')
    parser.add_argument('-s', '--scriptdir', type = str,
      help='Location of scripts')

