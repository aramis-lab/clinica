from os import PathLike
from typing import Optional

import click

from clinica.converters import cli_param

clinical_data_directory = click.option(
    "-cdd",
    "--clinical-data-dir",
    "clinical_data_directory",
    type=click.Path(exists=True, file_okay=False, resolve_path=True),
    help="Path to the clinical data directory.",
)

gif = click.option("-gif", is_flag=True, help="Add values from gif to session.tsv.")

full = click.option(
    "-full",
    is_flag=True,
    help="Add all clinical data (mandatory + optional) to sessions.tsv.",
)

clinical_data_txt = click.option(
    "-cdt",
    "--clinical-data-txt",
    "clinical_data_txt",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    help="Path to a txt file containing additional clinical data you want to have in the BIDS output.",
)


@click.command(name="genfi-to-bids")
@cli_param.dataset_directory
@cli_param.bids_directory
@clinical_data_directory
@gif
@full
@clinical_data_txt
def cli(
    dataset_directory: PathLike,
    bids_directory: PathLike,
    clinical_data_directory: Optional[PathLike] = None,
    clinical_data_txt: Optional[PathLike] = None,
    gif: Optional[bool] = False,
    full: Optional[bool] = False,
) -> None:
    """GENFI to BIDS converter.

    Convert the imaging and clinical data of GENFI, located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from ._converter import convert

    convert(
        dataset_directory,
        bids_directory,
        clinical_data_directory,
        gif=gif,
        full=full,
        path_to_clinical_txt=clinical_data_txt,
    )


if __name__ == "__main__":
    cli()
