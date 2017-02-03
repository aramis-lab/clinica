"""
The 'clinica' executable command line, installed with the clinica packages,
call this module.

The aim of this module is to execute pipeline from command line,
and give to the user some other utils to works with the pipelines.

"""


from __future__ import print_function
import argcomplete, sys, os, subprocess
from clinica.engine.cmdparser import *

__author__ = "Michael Bacci"
__copyright__ = "Copyright 2016,2017 The Aramis Lab Team"
__credits__ = ["Michael Bacci"]
__license__ = "see LICENSE.txt file"
__version__ = "1.0.0"
__maintainer__ = "Michael Bacci"
__email__ = "michael.bacci@inria.fr"
__status__ = "Development"


def visualize(clinicaWorkflow, ids, rebase=False):
    """Open a specific GUI program to display images made by pipeline

    :param clinicaWorkflow: the main pipeline object
    :param ids: list of id of patients
    :param rebase: path to looking for configuration
    """

    if not clinicaWorkflow.data.has_key('visualize'):
        print("No visualization was defined")
        exit(0)

    class chdir:
        def __init__(self, base):
            self.pwd = os.getcwd()
            os.chdir(base)
        def __del__(self):
            os.chdir(self.pwd)

    change_directory = None
    if rebase is False:
        change_directory = chdir(clinicaWorkflow.base_dir)
    else:
        change_directory = chdir(rebase)

    print(clinicaWorkflow.data['visualize'])
    program, arguments, matches = clinicaWorkflow.data['visualize']

    def run_program(id): subprocess.Popen([program] + arguments.replace("${%s}" % matches, id).strip().split(" "))
    [run_program(id) for id in ids]


def shell(clinicaWorkflow):
    """Open a python/ipython shell and re-init the clinicaWorkflow object

    :param clinicaWorkflow: the main pipeline object
    """

    workflow = clinicaWorkflow
    __banner__ = "workflow variable is instantiated for you!"
    namespace = globals().copy()
    namespace.update(locals())

    def load_python_shell():
        import readline
        import code
        shell = code.InteractiveConsole(namespace)
        shell.interact(banner=__banner__)

    def load_ipython_shell():
        from IPython.terminal.embed import InteractiveShellEmbed
        InteractiveShellEmbed(user_ns=namespace,banner1=__banner__)()

    try:
        load_ipython_shell()
    except:
        try:
            load_python_shell()
        except:
            print("Impossible to load ipython or python shell")

def load_conf(args):
    """Load a pipeline serialization

    :param args: the path where looking for
    :return: ClinicaWorkflow object
    """

    import cPickle

    def load(path):
        file = os.path.join(path, "clinica.pkl")
        if os.path.isfile(file): return cPickle.load(open(file))
        return False

    wk = False

    if len(args) == 0:
        wk = load(os.getcwd())
    elif os.path.isdir(args[0]):
        wk = load(args[0])

    if not wk:
        print("No <clinica.pkl> file found!")
        exit(0)

    return wk


class CmdlineCache():
    def __init__(self):
        self.converters = None
        self.io_options = None

    def setup_converters(self):
        from clinica.bids.load_cmdline_converter import load_cmdline_converters
        self.converters = load_cmdline_converters()

    def setup_io_options(self):
        self.io_options = [CmdParserSubsSess(), CmdParserMergeTsv(), CmdParserMissingModalities()]

    def load(self):
        self.setup_converters()
        self.setup_io_options()

class CmdlineHelper():
    def __init__(self):
        self.cmdline_cache = None

    def load_cache(self):
        import pickle
        from os.path import expanduser
        home_dir = expanduser("~")
        if not os.path.exists(home_dir):
            raise Exception("Home user can't be found!")
            exit(-1)

        clinica_dir = join(home_dir, ".clinica")
        if not os.path.exists(clinica_dir):
            try:
                os.makedirs(clinica_dir)
            except:
                raise Exception("Error: is not possible to create the [%s] dir" % clinica_dir)
                exit(-1)

        cmdline_cache_file = join(clinica_dir, "cmdline_cache.pkl")
        if os.path.isfile(cmdline_cache_file):
            pkl_cmdline_cache = open(cmdline_cache_file, 'rb')
            self.cmdline_cache = pickle.load(pkl_cmdline_cache)
            pkl_cmdline_cache.close()
        else:
            self.cmdline_cache = CmdlineCache()
            self.cmdline_cache.load()

        return self.cmdline_cache

def execute():

    cmdline_helper = CmdlineHelper()
    cmdline_cache = cmdline_helper.load_cache()

    """
    Define and parse the command line argument
    """

    parser = ArgumentParser()
    sub_parser = parser.add_subparsers()
    parser.add_argument("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="don't print status messages to stdout")

    """
    visualize option: open image[s] in a specific GUI program, generated by a pipeline
    """
    vis_parser = sub_parser.add_parser('visualize')
    vis_parser.add_argument("-i", "--id", dest="id",
                      required=True,
                      help="unique identifier")
    vis_parser.add_argument("-r", "--rebase", dest="rebase",
                      default=False,
                      help="unique identifier")
    def vis_parser_fun(args):
        visualize(load_conf(args[1:]), args.id.split(","), args.rebase)
    vis_parser.set_defaults(func=vis_parser_fun)

    """
    shell option: re-open a nipype.Workflow object within python/ipython session
    """
    shell_parser = sub_parser.add_parser('shell')
    def shell_parser_fun(args):
        shell(load_conf(args[1:]))
    shell_parser.set_defaults(func=shell_parser_fun)


    """
    pipelines-list option: show all available pipelines
    """
    pipeline_list_parser = sub_parser.add_parser('pipeline-list')
    def pipeline_list_fun(args):
        #display all available pipelines
        print(*get_cmdparser_names())
    pipeline_list_parser.set_defaults(func=pipeline_list_fun)

    """
    run option: run one of the available pipelines
    """
    run_parser = sub_parser.add_parser('run')
    #adding the independent pipeline ArgumentParser objects
    # init_cmdparser_objects(run_parser.add_subparsers())
    pipelines = [CmdParserT1SPMFullPrep(), CmdParserT1SPMSegment(), CmdParserPETPreprocessing(),
                 CmdParserT1FreeSurfer(), CmdParserT1FSL(),
                 CmdParserDWIPreprocessingFieldmapBased(),
                 CmdParserDWIT1Registration(), CmdParserDWIProcessing(),
                 CmdParserStatisticsSurfStat(), CmdParserMachineLearningVBLinearSVM()]
    init_cmdparser_objects(parser, run_parser.add_subparsers(), pipelines)

    """
    convert option: convert one of the supported dataset to the BIDS specification
    """
    convert_parser = sub_parser.add_parser('convert')
    from clinica.bids.load_cmdline_converter import load_cmdline_converters
    init_cmdparser_objects(parser, convert_parser.add_subparsers(), load_cmdline_converters())

    """
    io option
    """
    io_parser = sub_parser.add_parser('io')
    io_tasks = [CmdParserSubsSess(), CmdParserMergeTsv(), CmdParserMissingModalities()]
    init_cmdparser_objects(parser, io_parser.add_subparsers(), io_tasks)

    def silent_help(): pass

    def single_error_message(p):
        def error(x):
            p.print_help()
            parser.print_help = silent_help
            exit(-1)
        return error
    for p in [vis_parser, shell_parser, pipeline_list_parser, run_parser]: p.error = single_error_message(p)

    #Do not want stderr message
    def silent_msg(x): pass
    parser.error = silent_msg

    args = None
    try:
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(-1)

    if args is None or hasattr(args,'func') is False:
            parser.print_help()
            exit(-1)

    #Run the pipeline!
    args.func(args)


if __name__ == '__main__':
    execute()
