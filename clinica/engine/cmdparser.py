# coding: utf8

"""Define the command line parser for each utility (pipeline, iotools, etc.)
"""

import abc
from argparse import ArgumentParser
from os.path import join
from os import getcwd
from os.path import expanduser

from colorama import Fore
PIPELINE_CATEGORIES = {
    'CLINICA_COMPULSORY': '%sClinica mandatory arguments%s' % (Fore.YELLOW, Fore.RESET),
    'OPTIONAL': '%sPipeline options%s' % (Fore.YELLOW, Fore.RESET),
    'CLINICA_OPTIONAL': '%sClinica standard options%s' % (Fore.YELLOW, Fore.RESET),
    'ADVANCED': '%sPipeline advanced options%s' % (Fore.YELLOW, Fore.RESET),
}


class CmdParser:
    """Abstract class to extend in order to create your command line parser

    For pipelines, please use the 'clinica generate template' command.
    For converters, see clinica/iotools/converters/*/*_cli.py
    For iotools, see clinica/iotools/utils/data_handling_cli.py
    """
    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.reset()
        self.build()

    def build(self):
        self.define_name()
        self.define_description()
        self.set_content()
        self.define_options()

    def reset(self):
        self._args = ArgumentParser()
        self._name = None

    def set_content(self):
        self._args._positionals.title = '%sMandatory arguments%s' % (Fore.YELLOW, Fore.RESET)
        self._args._optionals.title = '%sOptional arguments%s' % (Fore.YELLOW, Fore.RESET)
        if self._description is None:
            self._description = self._name
            self._args.description = '%sIf you are not familiar with Clinica, see: http://clinica.run/doc/InteractingWithClinica/%s' % \
                                     (Fore.CYAN, Fore.RESET)
        else:
            self._args.description = '%s%s\n\nIf you are not familiar with Clinica, see: http://clinica.run/doc/InteractingWithClinica/%s' % \
                                     (Fore.CYAN, self._description, Fore.RESET)

    @property
    def options(self): return self._args

    @options.setter
    def options(self, x): self._args = x

    @property
    def name(self): return self._name

    @name.setter
    def name(self, x): self._name = x

    @property
    def description(self): return self._description

    @description.setter
    def description(self, x): self._description = x

    @abc.abstractmethod
    def define_name(self): pass

    def define_description(self):
        self._description = None

    @abc.abstractmethod
    def define_options(self): pass

    @abc.abstractmethod
    def run_command(self, args): pass

    def absolute_path(self, arg):
        if arg is None:
            return None
        elif arg[:1] == '~':
            return expanduser(arg)
        elif arg[:1] == '.':
            return getcwd()
        else:
            return join(getcwd(), arg)


def init_cmdparser_objects(root_parser, parser, objects):
    """
    Init all derived CmdParser instances with specific data.

    Args:
        root_parser: The root parser
        parser: The ArgParser node (e.g. 'run' or 'convert')
        objects: All CmdParser instances of this file
    """

    def silent_help():
        pass

    def error_message(p):
        def error(x):
            p.print_help()
            root_parser.print_help = silent_help
            exit(-1)
        return error

    def init(x):
        import argparse
        x.options = parser.add_parser(x.name,
                                      add_help=False,
                                      help=x.description,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
        x.options.error = error_message(x.options)
        x.options.set_defaults(func=x.run_command)
        x.build()

    for x in objects:
        try:
            init(x)
        except:
            pass


def get_cmdparser_names(objects=None):
    """
    Return the names of all pipelines

    Args:
        objects: All CmdParser instances of this file

    Returns:
        The names of all pipelines
    """
    if objects is None:
        objects = get_cmdparser_objects()
    for x in objects:
        yield x.name
