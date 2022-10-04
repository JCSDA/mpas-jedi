#!/usr/bin/env python3

jediAppName = 'variational'
nOuterIter = 0
ombgGroup  = 'ombg'
omanGroup  = 'oman'

def add_arguments(parser):
  global jediAppName, nOuterIter, ombgGroup, omanGroup
  parser.add_argument('-app', '--jediAppName', type=str, default=jediAppName, choices = ['variational','hofx'],
                  help='Name of jedi application that produced IODA files')
  parser.add_argument('-nout', '--nOuterIter', type=int, default=nOuterIter,
                  help='Number of outer iterations for jedi application')
  parser.add_argument('-bg', '--ombgGroup', type=str, default=ombgGroup,
                  help='ObsGroup for background departures')
  parser.add_argument('-an', '--omanGroup', type=str, default=omanGroup,
                  help='ObsGroup for analysis departures')

def parse_args(args):
  global jediAppName, nOuterIter, ombgGroup, omanGroup
  jediAppName = args.jediAppName
  nOuterIter = args.nOuterIter
  ombgGroup = args.ombgGroup
  omanGroup = args.omanGroup
