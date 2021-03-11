#!/usr/bin/env python3

jediAppName = 'variational'
nOuterIter = 0
depbgGroup  = 'depbg'
depanGroup  = 'depan'

def add_arguments(parser):
  global jediAppName, nOuterIter, depbgGroup, depanGroup
  parser.add_argument('-app', '--jediAppName', type=str, default=jediAppName, choices = ['variational','hofx'],
                  help='Name of jedi application that produced IODA files')
  parser.add_argument('-nout', '--nOuterIter', type=int, default=nOuterIter, 
                  help='Number of outer iterations for jedi application')
  parser.add_argument('-bg', '--depbgGroup', type=str, default=depbgGroup, 
                  help='ObsGroup for background departures')
  parser.add_argument('-an', '--depanGroup', type=str, default=depanGroup, 
                  help='ObsGroup for analysis departures')

def parse_args(args):
  global jediAppName, nOuterIter, depbgGroup, depanGroup
  jediAppName = args.jediAppName
  nOuterIter = args.nOuterIter
  depbgGroup = args.depbgGroup
  depanGroup = args.depanGroup
