from os import PathLike
from typing import Optional

import click

from clinica.iotools.converters import cli_param


@click.command(name="ixi-to-bids")
@cli_param.dataset_directory
@cli_param.bids_directory
@cli_param.clinical_data_directory
@cli_param.subjects_list
def cli(
    dataset_directory: PathLike,
    bids_directory: PathLike,
    clinical_data_directory: PathLike,
    subjects_list: Optional[PathLike] = None,
) -> None:
    """IXI to BIDS converter."""
    from .ixi_to_bids import convert

    convert(dataset_directory, bids_directory, clinical_data_directory, subjects_list)


if __name__ == "__main__":
    cli()
