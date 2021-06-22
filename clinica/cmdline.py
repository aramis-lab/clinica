# coding: utf8

"""The 'clinica' executable command line, installed with the clinica package, calls this module.

The aim of this module is to execute pipelines from the command line,
and gives to the user some other utils to work with the pipelines.
"""
import click

from .engine.template import cli as generate_cli
from .iotools.converters.cli import cli as convert_cli
from .iotools.utils.cli import cli as iotools_cli
from .pipelines.cli import cli as run_cli

CONTEXT_SETTINGS = dict(
    # Extend content width to avoid shortening of pipeline help.
    max_content_width=120,
    # Display help string with -h, in addition to --help.
    help_option_names=["-h", "--help"],
)


def setup_logging(verbosity: int = 0) -> None:
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
    console_handler.setFormatter(
        ColoredFormatter("%(log_color)s%(asctime)s:%(levelname)s:%(message)s")
    )
    logger.addHandler(console_handler)


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
@click.option("-v", "--verbose", "verbosity", count=True)
def cli(verbosity: int = 0) -> None:
    setup_logging(verbosity=verbosity)


cli.add_command(convert_cli)
cli.add_command(generate_cli)
cli.add_command(iotools_cli)
cli.add_command(run_cli)


if __name__ == "__main__":
    cli()
