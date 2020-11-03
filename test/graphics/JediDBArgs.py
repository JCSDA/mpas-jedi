#!/usr/bin/env python3

# static parameter
default_path='../Data'

# command-line parameters
obsFKey  = 'ObsSpace'
geoFKey  = 'geoval'
diagFKey = 'ObsDiag'
fPrefixes = {
  obsFKey:  'obsout',
  geoFKey:  'geoval',
  diagFKey: 'ydiags',
}

def add_arguments(parser):
  global obsFKey, geoFKey, diagFKey, fPrefixes
  parser.add_argument("-o", "--oPrefix", default = fPrefixes[obsFKey],
                      help="prefix for ObsSpace files")
  parser.add_argument("-g", "--gPrefix", default = fPrefixes[geoFKey],
                      help="prefix for GeoVaLs files")
  parser.add_argument("-d", "--dPrefix", default = fPrefixes[diagFKey],
                      help="prefix for ObsDiagnostics files")

def parse_args(args):
  global fPrefixes
  fPrefixes = {
    obsFKey:  args.oPrefix,
    geoFKey:  args.gPrefix,
    diagFKey: args.dPrefix,
  }
