#!/usr/bin/env python3

import jediApplicationArgs
from ProcessArgs import ProcessArgs
import textwrap

class analyzeStatsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()
    self.argProcessors += [jediApplicationArgs]

  @staticmethod
  def add_arguments(parser):
    parser.add_argument("-n", "--npan",
      help="Number of tasks/processors for analyses")
    parser.add_argument("-r", "--npread",
      help="Number of tasks/processors for reading StatsDB objects, defaults to npan")
    parser.add_argument("-d", "--diagSpace",
      help=textwrap.dedent(
        '''
        Select DiagSpaces with non-conservative matching
           e.g., amsua selects amsua_metop-a, amsua_n19, etc.
           default behavior is to select all DiagSpaces in config
        '''))
    parser.add_argument("-g", "--anGrp",
      help="Select a group of DiagSpaces (overridden by --diagSpace option)")
    parser.add_argument("-a", "--analysisType",
      help="Select a single analysisType")

processor = analyzeStatsArgs()

processor.processArgs()

args = processor.args
