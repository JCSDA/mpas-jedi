#!/usr/bin/env python3

import argparse
import jediApplicationArgs
import JediDBArgs
from ProcessArgs import ProcessArgs

class GenerateABEIFactorsArgs(ProcessArgs):
  def __init__(self):
    super().__init__()
    self.argProcessors += [jediApplicationArgs, JediDBArgs]

  @staticmethod
  def add_arguments(parser):
    parser.add_argument("datetime", type = str,
                        help = "Date for output files in YYYYMMDDHH format")
    parser.add_argument("-p", "--dbPath", default = JediDBArgs.default_path, type = str,
                        help="Path to deterministic or mean state UFO feedback files, default = "
                             +JediDBArgs.default_path)
    parser.add_argument("-i", "--IRInstruments", default = 'abi_g16', type = str,
                        help = "Comma-separated list of IR instruments for which to calculate ABEI factors")
    parser.add_argument("-m", "--modelGridFile", default = "./grid.nc", type = str,
                        help="Path to netcdf file containing MPAS grid description")
    parser.add_argument("-l", "--inflationFile", default = "ABEIlambda.nc", type =str,
                        help="Suffix of inflation factor file to be generated")
    parser.add_argument("-plot", "--plotLambda", default = False, type = bool,
                        help="Suffix of inflation factor file to be generated")

processor = GenerateABEIFactorsArgs()

processor.processArgs()

args = processor.args
