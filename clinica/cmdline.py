"""
The 'clinica' executable command line, installed with the clinica package,
calls this module.

The aim of this module is to execute pipelines from the command line,
and gives to the user some other utils to work with the pipelines.
"""

from __future__ import print_function

import os
import sys

import argcomplete

from clinica.engine.cmdparser import *
from clinica.utils.stream import cprint

__author__ = "Michael Bacci"
__copyright__ = "Copyright 2016,2017 The Aramis Lab Team"
__credits__ = ["Michael Bacci"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Michael Bacci"
__email__ = "michael.bacci@inria.fr"
__status__ = "Development"


class ClinicaClassLoader:
    """
    Load pipelines from a custom locations (general from $HOME/clinica)
    """
    from clinica.pipelines.engine import Pipeline

    def __init__(self, env='CLINICAPATH', baseclass=Pipeline, reg=r".*_cli\.py$", extra_dir=""):
        self.env = env
        self.baseclass = baseclass
        self.reg = reg
        self.extra_dir = extra_dir

    def load(self):
        import os
        pipeline_cli_parsers = []

        if not os.environ.has_key(self.env):
            return pipeline_cli_parsers

        clinica_pipelines_path = join(os.environ[self.env], self.extra_dir)
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
        return [os.path.join(path, file) for path in paths for file in os.listdir(path) if re.match(reg, file) is not None]


def execute():
    """
    Define and parse the command line argument
    """
    parser = ArgumentParser()
    sub_parser = parser.add_subparsers()
    parser.add_argument("-v", "--verbose",
                        dest='verbose',
                        action='store_true', default=False,
                        help='Verbose: print all messages to the console')
    parser.add_argument("-l", "--logname",
                        dest='logname', default="clinica.log",
                        help='Define the log file name (default: clinica.log)')

    """
    run option: run one of the available pipelines
    """
    from clinica.engine import CmdParser

    from clinica.pipelines.t1_freesurfer.t1_freesurfer_cli import T1FreeSurferCLI  # noqa
    from clinica.pipelines.t1_spm_segmentation.t1_spm_segmentation_cli import T1SPMSegmentationCLI  # noqa
    from clinica.pipelines.t1_spm_dartel.t1_spm_dartel_cli import T1SPMDartelCLI  # noqa
    from clinica.pipelines.t1_spm_dartel2mni.t1_spm_dartel2mni_cli import T1SPMDartel2MNICLI  # noqa
    from clinica.pipelines.t1_spm_full_prep.t1_spm_full_prep_cli import T1SPMFullPrepCLI  # noqa

    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_cli import DWIPreprocessingUsingT1CLI  # noqa
    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_cli import DWIPreprocessingUsingPhaseDiffFieldmapCLI # noqa
    from clinica.pipelines.dwi_processing.dwi_processing_cli import DWIProcessingCLI  # noqa

    from clinica.pipelines.fmri_preprocessing.fmri_preprocessing_cli import fMRIPreprocessingCLI  # noqa

    from clinica.pipelines.statistics_surfstat.statistics_surfstat_cli import StatisticsSurfstatCLI  # noqa

    from clinica.pipelines.pet_preprocess_volume.pet_preprocess_volume_cli import PETPreprocessVolumeCLI  # noqa

    run_parser = sub_parser.add_parser('run')
    pipelines = ClinicaClassLoader(baseclass=CmdParser, extra_dir="pipelines").load()
    pipelines += [
        T1FreeSurferCLI(),
        T1SPMSegmentationCLI(),
        T1SPMDartelCLI(),
        T1SPMDartel2MNICLI(),
        T1SPMFullPrepCLI(),
        DWIPreprocessingUsingT1CLI(),
        DWIPreprocessingUsingPhaseDiffFieldmapCLI(),
        DWIProcessingCLI(),
        fMRIPreprocessingCLI(),
        StatisticsSurfstatCLI(),
        CmdParserMachineLearningVBLinearSVM(),
        CmdParserMachineLearningSVMRB(),
        PETPreprocessVolumeCLI()
    ]
    init_cmdparser_objects(parser, run_parser.add_subparsers(), pipelines)

    """
    pipelines-list option: show all available pipelines
    """
    pipeline_list_parser = sub_parser.add_parser('pipelines-list')

    def pipeline_list_fun(args):
        for pipeline in pipelines:
            cprint(pipeline.name)

    pipeline_list_parser.set_defaults(func=pipeline_list_fun)

    """
    convert option: convert one of the supported datasets into the BIDS
    """
    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids_cli import AiblToBidsCLI  # noqa
    from clinica.iotools.converters.adni_to_bids.adni_to_bids_cli import AdniToBidsCLI  # noqa
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_cli import OasisToBidsCLI  # noqa

    convert_parser = sub_parser.add_parser('convert')
    converters = ClinicaClassLoader(baseclass=CmdParser,
                                    extra_dir="iotools/converters").load()
    converters += [
        AiblToBidsCLI(),
        AdniToBidsCLI(),
        OasisToBidsCLI()
    ]
    init_cmdparser_objects(parser, convert_parser.add_subparsers(), converters)

    """
    generate option: template
    """
    template_parser = sub_parser.add_parser('generate')

    from clinica.engine.template import CmdGenerateTemplates
    init_cmdparser_objects(parser, template_parser.add_subparsers(), [
        CmdGenerateTemplates()
    ])

    """
    iotools option
    """
    io_parser = sub_parser.add_parser('iotools')
    io_tools = [
        CmdParserSubsSess(),
        CmdParserMergeTsv(),
        CmdParserMissingModalities()
    ]
    init_cmdparser_objects(parser, io_parser.add_subparsers(), io_tools)

    """
    Silent all sub-parser errors methods except the one which is called
    otherwise the output console will display useless messages
    """
    def silent_help(): pass

    def single_error_message(p):
        def error(x):
            p.print_help()
            parser.print_help = silent_help
            exit(-1)
        return error
    for p in [pipeline_list_parser, run_parser]:
        p.error = single_error_message(p)

    # Do not want stderr message
    def silent_msg(x):
        pass
    parser.error = silent_msg

    """
    Parse the command and check that everything went fine
    """
    args = None
    unknown_args = None
    try:
        argcomplete.autocomplete(parser)
        args, unknown_args = parser.parse_known_args()
    except SystemExit:
        exit(-1)
    except Exception:
        parser.print_help()
        exit(-1)

    if unknown_args:
        raise ValueError('Unknown flag detected: %s' % unknown_args)

    if args is None or hasattr(args, 'func') is False:
        parser.print_help()
        exit(-1)

    import clinica.utils.stream as var
    var.clinica_verbose = args.verbose

    if args.verbose is False:
        """
        Enable only cprint(msg) --> clinica print(msg)
        - All the print() will be ignored!
        - All the logging will be redirect to the log file.
        """
        from clinica.utils.stream import FilterOut
        sys.stdout = FilterOut(sys.stdout)
        import logging as python_logging
        from logging import Filter, ERROR
        import os
        from nipype import config, logging
        from nipype import logging

        # Configure Nipype logger for our needs
        config.update_config({'logging': {'workflow_level': 'INFO',
                                          'log_directory': os.getcwd(),
                                          'log_to_file': True},
                              'execution': {'stop_on_first_crash': False,
                                            'hash_method': 'content'}
                              })
        logging.update_logging(config)

        # Define the LogFilter for ERROR detection
        class LogFilter(Filter):
            """
            The LogFilter class ables to monitor if an ERROR log signal is sent
            from Clinica/Nipype. If detected, the user will be warned.
            """
            def filter(self, record):
                if record.levelno >= ERROR:
                    cprint("An ERROR was generated: please check the log file for more information")
                return True

        logger = logging.getLogger('workflow')
        logger.addFilter(LogFilter())

        # Remove all handlers associated with the root logger object
        for handler in python_logging.root.handlers[:]:
            python_logging.root.removeHandler(handler)

        logging.disable_file_logging()

        # Enable file logging using a filename
        def enable_file_logging(self, filename):
            """
            Hack to define a filename for the log file! It overloads the
            'enable_file_logging' method in 'nipype/utils/logger.py' file.
            """
            import logging
            from logging.handlers import RotatingFileHandler as RFHandler
            config = self._config
            LOG_FILENAME = os.path.join(config.get('logging', 'log_directory'),
                                        filename)
            hdlr = RFHandler(LOG_FILENAME,
                             maxBytes=int(config.get('logging', 'log_size')),
                             backupCount=int(config.get('logging',
                                                        'log_rotate')))
            formatter = logging.Formatter(fmt=self.fmt, datefmt=self.datefmt)
            hdlr.setFormatter(formatter)
            self._logger.addHandler(hdlr)
            self._fmlogger.addHandler(hdlr)
            self._iflogger.addHandler(hdlr)
            self._hdlr = hdlr
        enable_file_logging(logging, args.logname)

        class Stream:
            def write(self, text):
                print(text)
                sys.stdout.flush()

        python_logging.basicConfig(
            format=logging.fmt, datefmt=logging.datefmt, stream=Stream())

    # Finally, run the pipelines
    args.func(args)


if __name__ == '__main__':
    execute()
