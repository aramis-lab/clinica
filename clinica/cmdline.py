# coding: utf8

"""The 'clinica' executable command line, installed with the clinica package, calls this module.

The aim of this module is to execute pipelines from the command line,
and gives to the user some other utils to work with the pipelines.
"""

import os
import sys
from argparse import ArgumentParser

import argcomplete

from clinica.engine.cmdparser import init_cmdparser_objects
from clinica.utils.stream import cprint


class ClinicaClassLoader:
    """Load pipelines from a custom locations (general from $HOME/clinica)."""

    def __init__(
        self, env="CLINICAPATH", baseclass=None, reg=r".*_cli\.py$", extra_dir=""
    ):
        self.env = env
        if baseclass is None:
            import clinica.pipelines.engine as cpe

            self.baseclass = cpe.Pipeline
        else:
            self.baseclass = baseclass
        self.reg = reg
        self.extra_dir = extra_dir

    def load(self):
        import os

        pipeline_cli_parsers = []

        if self.env not in os.environ.keys():

            return pipeline_cli_parsers

        clinica_pipelines_path = os.path.join(os.environ[self.env], self.extra_dir)
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
            if p not in sys.path:
                sys.path.append(p)

    def discover_path_with_subdir(self, path):
        return [
            os.path.join(path, file)
            for file in os.listdir(path)
            if os.path.isdir(os.path.join(path, file))
        ]

    def find_files(self, paths, reg):
        import re

        return [
            os.path.join(path, file)
            for path in paths
            for file in os.listdir(path)
            if re.match(reg, file) is not None
        ]


# Nice display
def custom_traceback(exc_type, exc_value, exc_traceback):
    import math
    import traceback

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.stream import cprint

    if issubclass(exc_type, ClinicaException):
        cprint(exc_value)
    elif issubclass(exc_type, KeyboardInterrupt):
        cprint(
            f"\n{Fore.RESET}[Error] Program interrupted by the user. Clinica will now exit...{Fore.RESET}"
        )
    else:
        cprint(
            f"{Fore.RED}\n{'*' * 23}\n*** Clinica crashed ***\n{'*' * 23}\n{Fore.RESET}"
        )
        cprint(f"{Fore.YELLOW}Exception type:{Fore.RESET} {exc_type.__name__}")
        cprint(f"{Fore.YELLOW}Exception value:{Fore.RESET} {exc_value}")
        cprint(
            "Below are displayed information that were gathered when Clinica crashed. "
            "This will help to understand what happened if you transfer "
            "those information to the Clinica development team.\n"
        )

        frames = traceback.extract_tb(exc_traceback)
        framewidth = int(math.ceil(math.log(len(frames)) / math.log(10)))
        filewidth = 0
        linewidth = 0
        functionwidth = 0
        for frame in frames:
            filewidth = max(filewidth, len(frame[0]))
            linewidth = max(linewidth, frame[1])
            functionwidth = max(functionwidth, len(frame[2]))
        linewidth = int(math.ceil(math.log(linewidth) / math.log(10)))
        cprint("=" * (filewidth + linewidth + functionwidth + linewidth))
        for i in range(len(frames)):
            t = (
                "{}"
                + " " * (1 + framewidth - len(str(i)))
                + Fore.RED
                + "{}"
                + " " * (1 + filewidth - len(frames[i][0]))
                + Fore.RESET
                + "{}"
                + " " * (1 + linewidth - len(str(frames[i][1])))
                + Fore.GREEN
                + "{}"
                + " " * (1 + functionwidth - len(frames[i][2]))
                + Fore.RESET
                + "{}"
            )
            cprint(t.format(i, frames[i][0], frames[i][1], frames[i][2], frames[i][3]))
        cprint("=" * (filewidth + linewidth + functionwidth + linewidth))

    if not issubclass(exc_type, KeyboardInterrupt):
        cprint(
            f"\nDocumentation can be found here: "
            f"{Fore.BLUE}https://aramislab.paris.inria.fr/clinica/docs/public/latest/{Fore.RESET}\n"
            f"If you need support, do not hesitate to ask: "
            f"{Fore.BLUE}https://groups.google.com/forum/#!forum/clinica-user{Fore.RESET}\n"
            f"Alternatively, you can also open an issue on GitHub: "
            f"{Fore.BLUE}https://github.com/aramis-lab/clinica/issues{Fore.RESET}"
        )


def execute():
    import argparse
    import logging
    import os
    import warnings

    from colorama import Fore

    # Suppress potential warnings
    warnings.filterwarnings("ignore")

    # Remove warnings from duecredit package in order to silent
    # "Assuming non interactive session since isatty found missing" message
    logging.getLogger("duecredit.utils").setLevel(logging.ERROR)

    # Add warning message if PYTHONPATH is not empty
    # cf https://groups.google.com/forum/#!topic/clinica-user/bVgifEdkg20
    python_path = os.environ.get("PYTHONPATH", "")
    if python_path:
        print(
            f"{Fore.YELLOW}[Warning] The PYTHONPATH environment variable is not empty."
            f" Make sure there is no interference with Clinica "
            f"(content of PYTHONPATH: {python_path}).{Fore.RESET}"
        )

    # Nice traceback when clinica crashes
    sys.excepthook = custom_traceback

    OPTIONAL_TITLE = f"{Fore.YELLOW}Optional arguments{Fore.RESET}"
    """
    Define and parse the command line argument
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-h", "--help", action="help", default=argparse.SUPPRESS, help=argparse.SUPPRESS
    )
    parser._positionals.title = (
        f"{Fore.YELLOW}clinica expects one of the following keywords{Fore.RESET}"
    )
    parser._optionals.title = OPTIONAL_TITLE

    sub_parser = parser.add_subparsers(metavar="")
    parser.add_argument(
        "-V",
        "--version",
        dest="version",
        action="store_true",
        default=False,
        help="Clinica's installed version",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Verbose: print all messages to the standard output",
    )
    parser.add_argument(
        "-l",
        "--logname",
        dest="logname",
        default="clinica.log",
        metavar=("file.log"),
        help="Define the log file name (default: clinica.log)",
    )

    """
    run category: run one of the available pipelines
    """
    # fmt: off
    from clinica.engine import CmdParser
    from clinica.pipelines.deeplearning_prepare_data.deeplearning_prepare_data_cli import (
        DeepLearningPrepareDataCLI,
    )
    from clinica.pipelines.dwi_connectome.dwi_connectome_cli import DwiConnectomeCli
    from clinica.pipelines.dwi_dti.dwi_dti_cli import DwiDtiCli
    from clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_cli import (
        DwiPreprocessingUsingPhaseDiffFieldmapCli,
    )
    from clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_cli import (
        DwiPreprocessingUsingT1Cli,
    )
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_cli import (
        SpatialSVMCLI,
    )
    from clinica.pipelines.pet_linear.pet_linear_cli import PETLinearCLI
    from clinica.pipelines.pet_surface.pet_surface_cli import PetSurfaceCLI
    from clinica.pipelines.pet_surface.pet_surface_longitudinal_cli import (
        PetSurfaceLongitudinalCLI,
    )
    from clinica.pipelines.pet_volume.pet_volume_cli import PETVolumeCLI
    from clinica.pipelines.statistics_surface.statistics_surface_cli import (
        StatisticsSurfaceCLI,
    )
    from clinica.pipelines.statistics_volume.statistics_volume_cli import (
        StatisticsVolumeCLI,
    )
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_cli import (
        StatisticsVolumeCorrectionCLI,
    )
    from clinica.pipelines.t1_freesurfer.t1_freesurfer_cli import T1FreeSurferCLI
    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_cli import (
        T1FreeSurferLongitudinalCLI,
    )
    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_cli import (
        T1FreeSurferLongitudinalCorrectionCLI,
    )
    from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_cli import (
        T1FreeSurferTemplateCLI,
    )
    from clinica.pipelines.t1_linear.t1_linear_cli import T1LinearCLI
    from clinica.pipelines.t1_volume.t1_volume_cli import T1VolumeCLI
    from clinica.pipelines.t1_volume_create_dartel.t1_volume_create_dartel_cli import (
        T1VolumeCreateDartelCLI,
    )
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_cli import (
        T1VolumeDartel2MNICLI,
    )
    from clinica.pipelines.t1_volume_existing_template.t1_volume_existing_template_cli import (
        T1VolumeExistingTemplateCLI,
    )
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_cli import (
        T1VolumeParcellationCLI,
    )
    from clinica.pipelines.t1_volume_register_dartel.t1_volume_register_dartel_cli import (
        T1VolumeRegisterDartelCLI,
    )
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_cli import (
        T1VolumeTissueSegmentationCLI,
    )

    # fmt: on

    pipelines = ClinicaClassLoader(baseclass=CmdParser, extra_dir="pipelines").load()
    # The order in pipelines var will the same when typing `clinica run`
    # Pipelines are sorted by main / advanced pipelines then by modality
    pipelines += [
        # Main pipelines:
        T1FreeSurferCLI(),
        T1VolumeCLI(),
        T1FreeSurferLongitudinalCLI(),
        T1LinearCLI(),
        DwiPreprocessingUsingPhaseDiffFieldmapCli(),
        DwiPreprocessingUsingT1Cli(),
        DwiDtiCli(),
        DwiConnectomeCli(),
        PETLinearCLI(),
        PETVolumeCLI(),
        PetSurfaceCLI(),
        PetSurfaceLongitudinalCLI(),
        DeepLearningPrepareDataCLI(),
        SpatialSVMCLI(),
        StatisticsSurfaceCLI(),
        StatisticsVolumeCLI(),
        StatisticsVolumeCorrectionCLI(),
        # Advanced pipelines:
        T1VolumeExistingTemplateCLI(),
        T1VolumeTissueSegmentationCLI(),
        T1VolumeCreateDartelCLI(),
        T1VolumeRegisterDartelCLI(),
        T1VolumeDartel2MNICLI(),
        T1VolumeParcellationCLI(),
        T1FreeSurferTemplateCLI(),
        T1FreeSurferLongitudinalCorrectionCLI(),
    ]

    run_parser = sub_parser.add_parser(
        "run",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
        help="To run pipelines on BIDS/CAPS datasets.",
    )
    run_parser.description = (
        f"{Fore.GREEN}Run pipelines on BIDS/CAPS datasets.{Fore.RESET}"
    )
    run_parser._positionals.title = (
        f"{Fore.GREEN}clinica run expects one of the following pipelines{Fore.RESET}"
    )

    init_cmdparser_objects(
        parser, run_parser.add_subparsers(metavar="", dest="run"), pipelines
    )

    """
    convert category: convert one of the supported datasets into BIDS hierarchy
    """
    # fmt: off
    from clinica.iotools.converters.adni_to_bids.adni_to_bids_cli import AdniToBidsCLI
    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids_cli import AiblToBidsCLI
    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids_cli import NifdToBidsCLI
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_cli import (
        OasisToBidsCLI,
    )

    # fmt: on

    converters = ClinicaClassLoader(
        baseclass=CmdParser, extra_dir="iotools/converters"
    ).load()
    converters += [
        AdniToBidsCLI(),
        AiblToBidsCLI(),
        OasisToBidsCLI(),
        NifdToBidsCLI(),
    ]

    convert_parser = sub_parser.add_parser(
        "convert",
        add_help=False,
        help="To convert unorganized datasets into a BIDS hierarchy.",
    )
    convert_parser.description = f"{Fore.GREEN}Tools to convert unorganized datasets into a BIDS hierarchy.{Fore.RESET}"
    convert_parser._positionals.title = f"{Fore.YELLOW}clinica convert expects one of the following datasets{Fore.RESET}"
    convert_parser._optionals.title = OPTIONAL_TITLE
    init_cmdparser_objects(
        parser, convert_parser.add_subparsers(metavar="", dest="convert"), converters
    )

    """
    iotools category
    """
    from clinica.iotools.utils.data_handling_cli import (
        CmdParserCenterNifti,
        CmdParserMergeTsv,
        CmdParserMissingModalities,
        CmdParserMissingProcessing,
        CmdParserSubjectsSessions,
    )

    io_tools = [
        CmdParserSubjectsSessions(),
        CmdParserMergeTsv(),
        CmdParserMissingProcessing(),
        CmdParserMissingModalities(),
        CmdParserCenterNifti(),
    ]

    HELP_IO_TOOLS = "Tools to handle BIDS/CAPS datasets."
    io_parser = sub_parser.add_parser(
        "iotools",
        add_help=False,
        help=HELP_IO_TOOLS,
    )
    io_parser.description = f"{Fore.GREEN}{HELP_IO_TOOLS}{Fore.RESET}"
    io_parser._positionals.title = f"{Fore.YELLOW}clinica iotools expects one of the following BIDS/CAPS utilities{Fore.RESET}"
    io_parser._optionals.title = OPTIONAL_TITLE

    init_cmdparser_objects(
        parser, io_parser.add_subparsers(metavar="", dest="iotools"), io_tools
    )

    """
    visualize category: run one of the available pipelines
    """
    from clinica.engine import CmdParser
    from clinica.pipelines.t1_freesurfer.t1_freesurfer_visualizer import (
        T1FreeSurferVisualizer,
    )

    visualizers = ClinicaClassLoader(baseclass=CmdParser, extra_dir="pipelines").load()
    visualizers += [
        T1FreeSurferVisualizer(),
    ]

    visualize_parser = sub_parser.add_parser(
        "visualize",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
        help="To visualize outputs of Clinica pipelines.",
    )
    visualize_parser.description = (
        f"{Fore.GREEN}Visualize outputs of Clinica pipelines.{Fore.RESET}"
    )
    visualize_parser._positionals.title = f"{Fore.YELLOW}clinica visualize expects one of the following pipelines{Fore.RESET}"

    init_cmdparser_objects(
        parser,
        visualize_parser.add_subparsers(metavar="", dest="visualize"),
        visualizers,
    )

    """
    generate category: template
    """
    generate_parser = sub_parser.add_parser(
        "generate",
        add_help=False,
        help=(
            "To generate pre-filled files when creating new pipelines (for developers)."
        ),
    )
    generate_parser.description = f"{Fore.GREEN}Generate pre-filled files when creating new pipelines (for  developers).{Fore.RESET}"
    generate_parser._positionals.title = (
        f"{Fore.YELLOW}clinica generate expects one of the following tools{Fore.RESET}"
    )
    generate_parser._optionals.title = OPTIONAL_TITLE

    from clinica.engine.template import CmdGenerateTemplates

    init_cmdparser_objects(
        parser,
        generate_parser.add_subparsers(metavar="", dest="generate"),
        [CmdGenerateTemplates()],
    )

    """
    Silent all sub-parser errors methods except the one which is called
    otherwise the output console will display useless messages
    """

    def silent_help():
        pass

    def single_error_message(p):
        def error(x):
            from colorama import Fore

            print(f"{Fore.RED}Error {x}{Fore.RESET}\n")
            p.print_help()
            parser.print_help = silent_help
            exit(-1)

        return error

    for p in [run_parser, io_parser, convert_parser, generate_parser]:
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
        if args.version:
            import clinica

            print(f"Clinica version is: {clinica.__version__}")
            exit(0)
        if unknown_args:
            if ("--verbose" in unknown_args) or ("-v" in unknown_args):
                cprint("Verbose detected")
                args.verbose = True
            unknown_args = [i for i in unknown_args if i != "-v"]
            unknown_args = [i for i in unknown_args if i != "--verbose"]
            if unknown_args:
                print(
                    f"{Fore.YELLOW}[Warning] Unknown flag(s) detected: {unknown_args}. "
                    f"This will be ignored by Clinica{Fore.RESET}"
                )
    except SystemExit:
        exit(0)
    except Exception:
        print(
            f"{Fore.RED}\n[Error] You wrote wrong arguments on the command line. "
            f"Clinica will now exit.\n{Fore.RESET}"
        )
        parser.print_help()
        exit(-1)

    if "run" in args and hasattr(args, "func") is False:
        # Case when we type `clinica run` on the terminal
        run_parser.print_help()
        exit(0)
    elif "convert" in args and hasattr(args, "func") is False:
        # Case when we type `clinica convert` on the terminal
        convert_parser.print_help()
        exit(0)
    elif "iotools" in args and hasattr(args, "func") is False:
        # Case when we type `clinica iotools` on the terminal
        io_parser.print_help()
        exit(0)
    elif "visualize" in args and hasattr(args, "func") is False:
        # Case when we type `clinica visualize` on the terminal
        visualize_parser.print_help()
        exit(0)
    elif "generate" in args and hasattr(args, "func") is False:
        # Case when we type `clinica generate` on the terminal
        generate_parser.print_help()
        exit(0)
    elif args is None or hasattr(args, "func") is False:
        # Case when we type `clinica` on the terminal
        parser.print_help()
        exit(0)

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
        import os
        from logging import ERROR, Filter

        from nipype import config, logging

        # Configure Nipype logger for our needs
        config.update_config(
            {
                "logging": {
                    "workflow_level": "INFO",
                    "log_directory": os.getcwd(),
                    "log_to_file": True,
                },
                "execution": {"stop_on_first_crash": False, "hash_method": "content"},
            }
        )
        logging.update_logging(config)

        # Define the LogFilter for ERROR detection
        class LogFilter(Filter):
            """Monitor if an ERROR log signal is sent from Clinica/Nipype.

            If detected, the user will be warned.
            """

            def filter(self, record):
                if record.levelno >= ERROR:
                    import datetime

                    now = datetime.datetime.now().strftime("%H:%M:%S")
                    cprint(
                        f"{Fore.RED}[{now}]{Fore.RESET} An error was found: please "
                        f"check the log file ({args.logname}) or crash file "
                        f"(nipypecli crash <pklz_file>) for details. Clinica is still running."
                    )
                return True

        if args.verbose:
            logger = logging.getLogger("nipype.workflow")
            logger.addFilter(LogFilter())

        # Remove all handlers associated with the root logger object
        for handler in python_logging.root.handlers[:]:
            python_logging.root.removeHandler(handler)

        logging.disable_file_logging()

        # Enable file logging using a filename
        def enable_file_logging(self, filename):
            """Hack to define a filename for the log file.

            It overloads the 'enable_file_logging' method in 'nipype/utils/logger.py' file.
            """
            import logging
            from logging.handlers import RotatingFileHandler as RFHandler

            config = self._config
            LOG_FILENAME = os.path.join(
                config.get("logging", "log_directory"), filename
            )
            hdlr = RFHandler(
                LOG_FILENAME,
                maxBytes=int(config.get("logging", "log_size")),
                backupCount=int(config.get("logging", "log_rotate")),
            )
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
            format=logging.fmt, datefmt=logging.datefmt, stream=Stream()
        )

    # Finally, run the command
    args.func(args)


if __name__ == "__main__":
    execute()
