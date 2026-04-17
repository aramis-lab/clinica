from os import PathLike
from typing import Optional

import click

from clinica.converters import cli_param


@click.command(name="oasis3-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
@cli_param.subjects_list
def cli(
    dataset_directory: PathLike,
    clinical_data_directory: PathLike,
    bids_directory: PathLike,
    subjects_list: Optional[PathLike] = None,
) -> None:
    """OASIS3 to BIDS converter.

    Convert the imaging and clinical data of OASIS3 (https://sites.wustl.edu/oasisbrains/),
    located in DATASET_DIRECTORY and CLINICAL_DATA_DIRECTORY respectively,
    to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from ._converter import convert

    convert(
        dataset_directory,
        bids_directory,
        clinical_data_directory,
        subjects=subjects_list,
    )


if __name__ == "__main__":
    cli()
