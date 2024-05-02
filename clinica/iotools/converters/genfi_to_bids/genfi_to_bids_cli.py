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
clinical_data_tsv = click.option(
    "-cdt",
    "--clinical-data-tsv",
    "clinical_data_tsv",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    help="Path to a tsv file containing additional clinical data you want to have in the BIDS output.",
)


@click.command(name="genfi-to-bids")
@cli_param.dataset_directory
@cli_param.bids_directory
@clinical_data_directory
@gif
@clinical_data_tsv
def cli(
    dataset_directory: PathLike,
    bids_directory: PathLike,
    clinical_data_directory: Optional[PathLike] = None,
    clinical_data_tsv: Optional[PathLike] = None,
    gif: Optional[bool] = False,
) -> None:
    """GENFI to BIDS converter.

    Convert the imaging and clinical data of GENFI, located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from pathlib import Path

    from .genfi_to_bids import convert

    convert(
        Path(dataset_directory),
        Path(bids_directory),
        Path(clinical_data_directory) if clinical_data_directory else None,
        gif=gif,
        path_to_clinical_tsv=Path(clinical_data_tsv) if clinical_data_tsv else None,
    )


if __name__ == "__main__":
    cli()
