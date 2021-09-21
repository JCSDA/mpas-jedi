#!/usr/bin/env python3

import argparse

class ProcessArgs:
  def __init__(self):
    self.argProcessors = [self]
    self.parser = argparse.ArgumentParser()

  def processArgs(self):
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
