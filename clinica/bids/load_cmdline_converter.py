import clinica.engine.cmdparser
import os
from os.path import realpath, split, join,basename
import re
import sys
import inspect
import importlib

def load_cmdline_converters():
    base_converters = join(split(realpath(__file__))[0], 'converters')
    for file in os.listdir(base_converters):
        if file.endswith(".py"):
            if re.match('__init__',file):
                continue
            if re.match('(.*)\.pyc$',file):
                continue

            f = "clinica.bids.converters.%s" % file.split('.')[0]
            mod = importlib.import_module(f)
            for name, obj in inspect.getmembers(mod):
                if name != 'CmdParser' and name != 'Converter' and inspect.isclass(obj):
                    x = obj()
                    if isinstance(x, clinica.engine.cmdparser.CmdParser):
                        yield x

