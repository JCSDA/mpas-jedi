#!/usr/bin/env python3

import argparse
import jediApplicationArgs
import JediDBArgs
from ProcessArgs import ProcessArgs

class DiagnoseObsStatisticsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()
    self.argProcessors += [jediApplicationArgs, JediDBArgs]

  @staticmethod
  def add_arguments(parser):
    parser.add_argument("-n", "--nprocs", default = 1, type = int,
                        help="Number of tasks/processors for multiprocessing")
    parser.add_argument("-p", "--meanPath", default = JediDBArgs.default_path, type = str,
                        help="Path to deterministic or mean state UFO feedback files, default = "
                             +JediDBArgs.default_path)
    parser.add_argument("-m", "--nMembers", default = 0, type = int,
                        help="number of ensemble members; must be >1 to produce ensemble diagnostic stats")
    parser.add_argument("-e", "--ensemblePath", default = JediDBArgs.default_path+"/mem{:03d}", type = str,
                        help="Path to ensemble member UFO feedback files; must have substitution string for member integer, e.g., '{:03d}' for 001, 002, etc...")

processor = DiagnoseObsStatisticsArgs()

processor.processArgs()

args = processor.args
