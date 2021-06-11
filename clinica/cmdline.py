# coding: utf8

"""The 'clinica' executable command line, installed with the clinica package, calls this module.

The aim of this module is to execute pipelines from the command line,
and gives to the user some other utils to work with the pipelines.
"""

import os
import sys
from argparse import ArgumentParser
from typing import Optional

import argcomplete

from clinica.engine.cmdparser import init_cmdparser_objects


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


def setup_logging(verbosity: Optional[int] = 0) -> None:
    """
    Setup Clinica's logging facilities.

    Args:
        verbosity (int): The desired level of verbosity for logging.
            (0 (default): WARNING, 1: INFO, 2: DEBUG)
    """
    from logging import DEBUG, INFO, WARNING, getLogger
    from sys import stdout

    from colorlog import ColoredFormatter, StreamHandler

    # Cap max verbosity level to 2.
    verbosity = min(verbosity, 2)

    # Define the module level logger.
    logger = getLogger("clinica")
    logger.setLevel([WARNING, INFO, DEBUG][verbosity])

    # Add console handler with custom formatting.
    console_handler = StreamHandler(stdout)
    console_handler.setFormatter(ColoredFormatter("%(log_color)s%(asctime)s:%(levelname)s:%(message)s"))
    logger.addHandler(console_handler)


# Nice display
def custom_traceback(exc_type, exc_value, exc_traceback):
    import math
    import traceback

    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.stream import cprint

    if issubclass(exc_type, ClinicaException):
        cprint(exc_value, lvl="error")
    elif issubclass(exc_type, KeyboardInterrupt):
        cprint(
            msg=f"Program interrupted by the user. Clinica will now exit...",
            lvl="warning",
        )
    else:
        cprint(
            msg=(
                "{'*' * 23}\n*** Clinica crashed ***\n{'*' * 23}\n\n"
                f"Exception type: {exc_type.__name__}\n"
                f"Exception value: {exc_value}\n\n"
                "Below are displayed information that were gathered when Clinica crashed. "
                "This will help to understand what happened if you transfer "
                "those information to the Clinica development team."
            ),
            lvl="error",
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
        cprint("=" * (filewidth + linewidth + functionwidth + linewidth), lvl="error")
        for i in range(len(frames)):
            t = (
                "{}"
                + " " * (1 + framewidth - len(str(i)))
                + "{}"
                + " " * (1 + filewidth - len(frames[i][0]))
                + "{}"
                + " " * (1 + linewidth - len(str(frames[i][1])))
                + "{}"
                + " " * (1 + functionwidth - len(frames[i][2]))
                + "{}"
            )
            cprint(t.format(i, frames[i][0], frames[i][1], frames[i][2], frames[i][3]), lvl="error")
        cprint("=" * (filewidth + linewidth + functionwidth + linewidth), lvl="error")

    if not issubclass(exc_type, KeyboardInterrupt):
        cprint(
            msg=(
                "Documentation can be found here:\n"
                "https://aramislab.paris.inria.fr/clinica/docs/public/latest/\n"
                "If you need support, do not hesitate to ask:\n"
                "https://groups.google.com/forum/#!forum/clinica-user\n"
                "Alternatively, you can also open an issue on GitHub:\n"
                "https://github.com/aramis-lab/clinica/issues\n"
            ),
            lvl="warning",
        )


def execute():
    import argparse
    import logging
    import os
    import warnings

    from clinica.utils.stream import cprint

    # Suppress potential warnings
    warnings.filterwarnings("ignore")

    # Remove warnings from duecredit package in order to silent
    # "Assuming non interactive session since isatty found missing" message
    logging.getLogger("duecredit.utils").setLevel(logging.ERROR)

    # Nice traceback when clinica crashes
    sys.excepthook = custom_traceback

    OPTIONAL_TITLE = "Optional arguments"
    """
    Define and parse the command line argument
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-h", "--help", action="help", default=argparse.SUPPRESS, help=argparse.SUPPRESS
    )
    parser._positionals.title = (
        f"clinica expects one of the following keywords"
    )
    parser._optionals.title = OPTIONAL_TITLE

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
        action="count",
        default=0,
        help="Increase logging verbosity",
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

    sub_parser = parser.add_subparsers(metavar="")

    run_parser = sub_parser.add_parser(
        "run",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
        help="To run pipelines on BIDS/CAPS datasets.",
    )
    run_parser.description = (
        f"Run pipelines on BIDS/CAPS datasets."
    )
    run_parser._positionals.title = (
        f"clinica run expects one of the following pipelines"
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
    convert_parser.description = f"Tools to convert unorganized datasets into a BIDS hierarchy."
    convert_parser._positionals.title = f"clinica convert expects one of the following datasets"
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
    io_parser.description = f"{HELP_IO_TOOLS}"
    io_parser._positionals.title = f"clinica iotools expects one of the following BIDS/CAPS utilities"
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
        f"Visualize outputs of Clinica pipelines."
    )
    visualize_parser._positionals.title = f"clinica visualize expects one of the following pipelines"

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
    generate_parser.description = f"Generate pre-filled files when creating new pipelines (for  developers)."
    generate_parser._positionals.title = (
        f"clinica generate expects one of the following tools"
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
            print(f"Error: {x}\n")
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
    try:
        argcomplete.autocomplete(parser)
        args = parser.parse_args()
    except SystemExit:
        exit(0)
    except Exception:
        print("Wrong arguments provided, Clinica will now exit!")
        parser.print_help()
        exit(-1)

    if args.version:
        import clinica

        print(f"Clinica version is: {clinica.__version__}")
        exit(0)

    setup_logging(verbosity=args.verbose)

    # Add warning message if PYTHONPATH is not empty
    # cf https://groups.google.com/forum/#!topic/clinica-user/bVgifEdkg20
    python_path = os.environ.get("PYTHONPATH", "")
    if python_path:
        cprint(
            msg=(
                "The PYTHONPATH environment variable is not empty. "
                "Make sure there is no interference with Clinica "
                f"(content of PYTHONPATH: {python_path})."
            ),
            lvl="warning",
        )

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

    # Finally, run the command
    args.func(args)


if __name__ == "__main__":
    execute()
