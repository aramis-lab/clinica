import abc
from argparse import ArgumentParser

class CmdParser:
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.build()

    def build(self):
        self.reset()
        self.define_name()
        self.define_options()

    def reset(self):
        self._args = ArgumentParser()
        self._name = None

    @property
    def options(self): return self._args

    @property
    def name(self): return self._name

    @abc.abstractmethod
    def define_name(self): pass

    @abc.abstractmethod
    def define_options(self): pass

    @abc.abstractmethod
    def parse_args(self): pass

    @abc.abstractmethod
    def run_pipeline(self): pass

class CmdParserT1(CmdParser):

    def define_name(self):
        self._name = 'T1'

    def define_options(self):
        self._args.add_argument("-s", "--source")

    def parse_args(self):
        print "parse args"

    def run_pipeline(self):
        print "run pipeline"

