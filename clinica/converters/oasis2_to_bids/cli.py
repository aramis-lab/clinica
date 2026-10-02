from pathlib import Path

import click

from clinica.converters import cli_param


@click.command(name="oasis2-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: Path,
    clinical_data_directory: Path,
    bids_directory: Path,
) -> None:
    """OASIS2 to BIDS converter.

    Convert the imaging and clinical data of OASIS2 (https://sites.wustl.edu/oasisbrains/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY (oasis_longitudinal_demographics.csv) respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from ._converter import convert

    convert(dataset_directory, bids_directory, clinical_data_directory)


if __name__ == "__main__":
    cli()
