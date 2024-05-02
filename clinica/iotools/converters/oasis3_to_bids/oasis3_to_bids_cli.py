from os import PathLike
from pathlib import Path

import click

from clinica.iotools.converters import cli_param


@click.command(name="oasis3-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: PathLike,
    clinical_data_directory: PathLike,
    bids_directory: PathLike,
) -> None:
    """OASIS3 to BIDS converter.

    Convert the imaging and clinical data of OASIS3 (https://sites.wustl.edu/oasisbrains/),
    located in DATASET_DIRECTORY and CLINICAL_DATA_DIRECTORY respectively,
    to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from .oasis3_to_bids import convert

    convert(
        Path(dataset_directory), Path(bids_directory), Path(clinical_data_directory)
    )


if __name__ == "__main__":
    cli()
