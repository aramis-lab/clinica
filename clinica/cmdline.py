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
__license__ = "See LICENSE.txt file"
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


class ClinicaClassLoader:
    from clinica.pipeline.engine import Pipeline
    """
    Load pipelines from a custom locations (general from $HOME/clinica)
    """
    def __init__(self,env='CLINICAPATH',baseclass=Pipeline,reg=r".*_cli\.py$", extra_dir=""):
        self.env = env
        self.baseclass = baseclass
        self.reg = reg
        self.extra_dir = extra_dir

    def load(self):
        import os
        pipeline_cli_parsers = []

        if not os.environ.has_key(self.env):
            return pipeline_cli_parsers

        clinica_pipelines_path = join(os.environ[self.env],self.extra_dir)
        if not os.path.isdir(clinica_pipelines_path):
            return pipeline_cli_parsers

        src_path = self.discover_path_with_subdir(clinica_pipelines_path)
        self.add_to_python_path(src_path)
        files_match = self.find_files(src_path, self.reg)

        for file in files_match:
            pipeline_cli_parsers.append(self.load_class(self.baseclass, file))

        return pipeline_cli_parsers

    def load_class(self, baseclass, file):
        import imp
        import inspect
        py_module_name, ext = os.path.splitext(os.path.split(file)[-1])
        py_module = imp.load_source(py_module_name, file)
        for class_name, class_obj in inspect.getmembers(py_module, inspect.isclass):
            if inspect.isclass(class_obj) and not inspect.isabstract(class_obj):
                x = class_obj()
                if isinstance(x, baseclass):
                    return x

    def add_to_python_path(self, paths):
        for p in paths:
            if not p in sys.path:
                sys.path.append(p)


    def discover_path_with_subdir(self, path):
        return [os.path.join(path, file) for file in os.listdir(path) if os.path.isdir(os.path.join(path, file))]

    def find_files(self, paths, reg):
        import re
        return [os.path.join(path,file) for path in paths for file in os.listdir(path) if re.match(reg, file) is not None]

def load_modular_pipelines_parser():

    import imp
    import os
    import os.path as op
    import re
    import inspect

    pipeline_cli_parsers = []

    # Clinica path
    if 'CLINICAPATH' in os.environ:
        clinica_path = os.environ['CLINICAPATH']
    else:
        print("WARNING: Variable 'CLINICAPATH' is not defined.")
        return pipeline_cli_parsers

    # Current path
    if os.getcwd() not in clinica_path.split(':'):
        clinica_path = clinica_path + ':' + os.getcwd()

    # List pipeline directories and fetch CLI parser class of each pipeline
    for one_clinica_path in clinica_path.split(':'):
        if op.isdir(one_clinica_path):
            for pipeline_dir in os.listdir(one_clinica_path):
                pipeline_path = op.join(one_clinica_path, pipeline_dir)
                if not pipeline_path in sys.path:
                    sys.path.append(pipeline_path)
                if op.isdir(pipeline_path):
                    for pipeline_file in os.listdir(pipeline_path):
                        if re.match(r".*_cli\.py$", pipeline_file) is not None:
                            py_module_name, ext = op.splitext(op.split(pipeline_file)[-1])
                            py_module = imp.load_source(py_module_name, op.join(pipeline_path, pipeline_file))
                            for class_name, class_obj in inspect.getmembers(py_module, inspect.isclass):
                                if re.match(r".*CLI$", class_name) is not None:
                                    pipeline_cli_parsers.append(class_obj())

    return pipeline_cli_parsers


class CmdlineCache():
    def __init__(self):
        self.converters = None
        self.io_options = None

    def setup_converters(self):pass
        #from clinica.iotools.load_cmdline_converter import load_cmdline_converters
        #self.converters = load_cmdline_converters()

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

    #cmdline_helper = CmdlineHelper()
    #cmdline_cache = cmdline_helper.load_cache()

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
    TODO: complete for future release
    """
    #shell_parser = sub_parser.add_parser('shell')
    #def shell_parser_fun(args):
    #    shell(load_conf(args[1:]))
    #shell_parser.set_defaults(func=shell_parser_fun)


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
    from clinica.engine import CmdParser
    from clinica.pipeline.t1_spm_segmentation.t1_spm_segmentation_cli import T1SPMSegmentationCLI
    from clinica.pipeline.fmri_preprocessing.fmri_preprocessing_cli import fMRIPreprocessingCLI
    run_parser = sub_parser.add_parser('run')
    pipelines = ClinicaClassLoader(baseclass=CmdParser, extra_dir="pipelines").load()
    # pipelines = load_modular_pipelines_parser()
    pipelines = pipelines + [CmdParserT1SPMFullPrep(), CmdParserT1SPMSegment(),
                 CmdParserT1SPMDartelTemplate(), CmdParserPETPreprocessing(),
                 CmdParserT1FreeSurfer(), CmdParserT1FSL(),
                 CmdParserDWIPreprocessingPhaseDifferenceFieldmap(), CmdParserDWIPreprocessingTwoPhaseImagesFieldmap(),
                 CmdParserDWIPreprocessingT1Based(),
                 CmdParserDWIProcessing(), T1SPMSegmentationCLI(), fMRIPreprocessingCLI(),
                 CmdParserStatisticsSurfStat(), CmdParserMachineLearningVBLinearSVM(), CmdParserMachineLearningSVMRB()]
    init_cmdparser_objects(parser, run_parser.add_subparsers(), pipelines)

    """
    convert option: convert one of the supported dataset to the BIDS specification
    """
    converters = ClinicaClassLoader(baseclass=CmdParser, extra_dir="iotools/converters", reg=r".*_bids\.py$").load()
    if (len(converters)):
        convert_parser = sub_parser.add_parser('convert')
        init_cmdparser_objects(parser, convert_parser.add_subparsers(), converters)

    """
    generate option: template
    """
    convert_parser = sub_parser.add_parser('generate')
    from clinica.engine.template import CmdGenerateTemplates
    init_cmdparser_objects(parser, convert_parser.add_subparsers(), [CmdGenerateTemplates()])

    """
    iotools option
    """
    io_parser = sub_parser.add_parser('iotools')
    io_tasks = [CmdParserSubsSess(), CmdParserMergeTsv(), CmdParserMissingModalities()]
    init_cmdparser_objects(parser, io_parser.add_subparsers(), io_tasks)

    def silent_help(): pass

    def single_error_message(p):
        def error(x):
            p.print_help()
            parser.print_help = silent_help
            exit(-1)
        return error
    for p in [vis_parser, pipeline_list_parser, run_parser]: p.error = single_error_message(p)

    #Do not want stderr message
    def silent_msg(x): pass
    parser.error = silent_msg

    args = None
    import argparse
    try:
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
    except SystemExit:
        exit(-1)
    except Exception:
        parser.print_help()
        exit(-1)

    if args is None or hasattr(args,'func') is False:
            parser.print_help()
            exit(-1)

    #Run the pipeline!
    args.func(args)


if __name__ == '__main__':
    execute()
