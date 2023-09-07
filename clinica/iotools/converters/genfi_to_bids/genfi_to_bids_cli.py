from os import PathLike
from typing import Optional

import click

from clinica.iotools.converters import cli_param

clinical_data_directory = click.option(
    "-cdd",
    "--clinical-data-dir",
    "clinical_data_directory",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Path to the clinical data directory",
)

gif = click.option("-gif", is_flag=True, help="Add values from gif to session.tsv")


@click.command(name="genfi-to-bids")
@cli_param.dataset_directory
@cli_param.bids_directory
@clinical_data_directory
@gif
def cli(
    dataset_directory: PathLike,
    bids_directory: PathLike,
    clinical_data_directory: Optional[PathLike] = None,
    gif: bool = False,
) -> None:
    """GENFI to BIDS converter.

    Convert the imaging and clinical data of GENFI, located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.bids_utils import _write_bidsignore
    from clinica.iotools.converters.genfi_to_bids.genfi_to_bids import (
        GenfiToBidsConverter,
    )

    GenfiToBidsConverter(
        dataset_directory,
        bids_directory,
        clinical_data_directory,
        gif=gif,
    ).convert()
    _write_bidsignore(str(bids_directory))


if __name__ == "__main__":
    cli()
