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

processor = analyzeStatsArgs()

processor.processArgs()

args = processor.args
