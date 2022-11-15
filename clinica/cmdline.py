"""The 'clinica' executable command line, installed with the clinica package, calls this module.

The aim of this module is to execute pipelines from the command line,
and gives to the user some other utils to work with the pipelines.
"""
import os
import sys

import click

import clinica
from clinica.engine.template import cli as generate_cli
from clinica.iotools.converters.cli import cli as convert_cli
from clinica.iotools.utils.cli import cli as iotools_cli
from clinica.utils.exceptions import ClinicaException
from clinica.utils.stream import cprint

CONTEXT_SETTINGS = dict(
    # Extend content width to avoid shortening of pipeline help.
    max_content_width=160,
    # Display help string with -h, in addition to --help.
    help_option_names=["-h", "--help"],
)


def setup_logging(verbose: bool = False) -> None:
    """Setup Clinica's logging facilities.

    Args:
        verbosity (int): The desired level of verbosity for logging.
            (0 (default): WARNING, 1: INFO, 2: DEBUG)
    """
    import logging
    import sys
    from logging.handlers import RotatingFileHandler as RFHandler

    import nipype
    from colorlog import ColoredFormatter, StreamHandler

    logging_level = "DEBUG" if verbose else "INFO"

    # Root logger configuration.
    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)

    # Clinica logger configuration.
    clinica_logger = logging.getLogger("clinica")
    clinica_logger.setLevel(logging_level)
    console_handler = StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(
        ColoredFormatter("%(log_color)s%(asctime)s:%(levelname)s:%(message)s")
    )
    clinica_logger.addHandler(console_handler)

    # Nipype logger configuration.
    # Monkey-patch nipype to use Python's RFH logger.
    nipype.utils.logger.RFHandler = RFHandler
    # Setup debug logging to file.
    nipype.config.enable_debug_mode()
    nipype.config.update_config(
        {
            "logging": {
                "log_directory": os.getenv("CLINICA_LOGGING_DIR", os.getcwd()),
                "log_to_file": True,
            },
            "execution": {"check_version": False},
        },
    )
    nipype.logging.update_logging(nipype.config)
    # Disable nipype logging to console.
    nipype_logger = logging.getLogger("nipype")
    nipype_logger.removeHandler(nipype_logger.handlers[0])
    nipype_logger.addHandler(logging.NullHandler())


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=clinica.__version__)
@click.option("-v", "--verbose", is_flag=True, help="Increase logging verbosity.")
def cli(verbose: bool) -> None:
    setup_logging(verbose=verbose)


cli.add_command(convert_cli)
cli.add_command(generate_cli)
cli.add_command(iotools_cli)


def main() -> None:
    try:
        cli()
    except ClinicaException as e:
        cprint(msg=e, lvl="error")
        sys.exit(1)


if __name__ == "__main__":
    main()
