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
    'CLINICA_COMPULSORY': '%sClinica mandatory arguments%s' % (Fore.BLUE, Fore.RESET),
    'OPTIONAL': '%sPipeline options%s' % (Fore.BLUE, Fore.RESET),
    'CLINICA_OPTIONAL': '%sClinica standard options%s' % (Fore.BLUE, Fore.RESET),
    'ADVANCED': '%sPipelines advanced options%s' % (Fore.BLUE, Fore.RESET),
    'IOTOOLS_OPTIONS': '%sOptional arguments%s' % (Fore.BLUE, Fore.RESET)
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
        from colorama import Fore
        self._args._positionals.title = '%sMandatory arguments%s' % (Fore.BLUE, Fore.RESET)
        self._args._optionals.title = '%sOptional arguments%s' % (Fore.BLUE, Fore.RESET)
        if self._description is None:
            self._description = self._name
            self._args.description = ('%sIf you are not familiar with Clinica, see:\n'
                                      'http://clinica.run/doc/InteractingWithClinica/%s'
                                      % (Fore.GREEN, Fore.RESET))
        else:
            self._args.description = (
                '%s%s\n\nIf you are not familiar with Clinica, see:\n'
                'http://clinica.run/doc/InteractingWithClinica%s'
                % (Fore.GREEN, self._description, Fore.RESET)
            )

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

    def add_clinica_standard_arguments(self,
                                       add_tsv_flag=True,
                                       add_wd_flag=True,
                                       add_nprocs_flag=True,
                                       add_overwrite_flag=False,
                                       ):
        clinica_standard_options = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        if add_tsv_flag:
            clinica_standard_options.add_argument(
                "-tsv", "--subjects_sessions_tsv",
                help='TSV file containing a list of subjects with their sessions.')
        if add_wd_flag:
            clinica_standard_options.add_argument(
                "-wd", "--working_directory",
                help='Temporary directory to store pipelines intermediate results.')
        if add_nprocs_flag:
            clinica_standard_options.add_argument(
                "-np", "--n_procs",
                metavar='N', type=int,
                help='Number of cores used to run in parallel.')
        if add_overwrite_flag:
            clinica_standard_options.add_argument(
                "-overwrite", "--overwrite_outputs",
                action='store_true', default=False,
                help='Force overwrite of output files in CAPS folder.')

        return clinica_standard_options

    @abc.abstractmethod
    def run_command(self, args): pass

    @staticmethod
    def list_to_string(list):
        """Convert list (e.g. [8, 8, 8]) to string (e.g. '8 8 8')"""
        string_without_commas_and_brackets = ' '.join(str(item) for item in list)
        return string_without_commas_and_brackets

    @staticmethod
    def absolute_path(arg):
        if arg is None:
            return None
        elif arg[:1] == '~':
            return expanduser(arg)
        elif arg[:2] == './':
            return join(getcwd(), arg[2:])
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
        except BaseException:
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
