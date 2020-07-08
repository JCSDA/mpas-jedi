import logging
import os
import re
import sys

loglevel = logging.INFO

## setup logging
class _AnsiColorStreamHandler(logging.StreamHandler):
    DEFAULT = '\x1b[0m'
    RED     = '\x1b[31m'
    GREEN   = '\x1b[32m'
    YELLOW  = '\x1b[33m'
    CYAN    = '\x1b[36m'

    CRITICAL = RED
    ERROR    = RED
    WARNING  = YELLOW
    INFO     = GREEN
    DEBUG    = CYAN

    @classmethod
    def _get_color(cls, level):
        if level >= logging.CRITICAL:  return cls.CRITICAL
        elif level >= logging.ERROR:   return cls.ERROR
        elif level >= logging.WARNING: return cls.WARNING
        elif level >= logging.INFO:    return cls.INFO
        elif level >= logging.DEBUG:   return cls.DEBUG
        else:                          return cls.DEFAULT

    def __init__(self, stream=None):
        logging.StreamHandler.__init__(self, stream)

    def format(self, record):
        text = logging.StreamHandler.format(self, record)
        color = self._get_color(record.levelno)
        return color + text + self.DEFAULT


class _ExitOnExceptionHandler(_AnsiColorStreamHandler):
    def emit(self, record):
        super().emit(record)
        if record.levelno in (logging.ERROR, logging.CRITICAL):
            raise SystemExit(-1)

## import this module in order to set the following logging config
logging.basicConfig(
    level=loglevel,
    format='%(asctime)s [%(levelname)s] : %(name)s : %(message)s',
    datefmt='%Y-%m-%d_%H:%M:%S',
    handlers=[
        _ExitOnExceptionHandler(),
    ]
)
