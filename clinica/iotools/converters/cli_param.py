"""Common CLI arguments and options used by converters."""

import click

dataset_directory = click.argument(
    "dataset_directory",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
)

clinical_data_directory = click.argument(
    "clinical_data_directory",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
)

bids_directory = click.argument(
    "bids_directory",
    type=click.Path(writable=True, file_okay=False, resolve_path=True),
)

clinical_data_only = click.option(
    "-c",
    "--clinical-data-only",
    is_flag=True,
    help="Convert clinical data only.",
)

subjects_list = click.option(
    "-sl",
    "--subjects_list",
    type=click.Path(exists=True),
    help=(
        "Path to a text file containing a list of specific subjects to extract."
        "The expected format is one subject ID per line."
    ),
)
