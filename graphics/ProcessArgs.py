#!/usr/bin/env python3

import argparse

class ProcessArgs:
  # make the parser a class variable, to limit running processArgs to only one Args object.
  parser = None

  def __init__(self):
    self.argProcessors = [self]
    # guarantee the args data member is defined, even if processArgs doesn't initialize args.
    self.args = None

  def processArgs(self):
    '''
    Only process one Args parser.
    Some programs load more than one args parser,
    due to importing files which in turn import arg parsers.
    This assumes the desired arg parser is imported first.
    The main program should explicitly import the Args class it needs.
    '''
    if ProcessArgs.parser is None:
      ProcessArgs.parser = argparse.ArgumentParser()
      for processor in self.argProcessors:
        processor.add_arguments(self.parser)

      self.args = self.parser.parse_args()

      for processor in self.argProcessors:
        processor.parse_args(self.args)

  @staticmethod
  def add_arguments(parser):
    pass

  @staticmethod
  def parse_args(args):
    pass
