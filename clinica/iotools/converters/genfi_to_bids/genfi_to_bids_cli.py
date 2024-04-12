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

    from clinica.iotools.bids_utils import _write_bidsignore
    from clinica.utils.check_dependency import ThirdPartySoftware, check_software
    from clinica.utils.stream import cprint

    from .genfi_to_bids import convert_images

    check_software(ThirdPartySoftware.DCM2NIIX)
    bids_directory = Path(bids_directory)
    clinical_data_directory = (
        Path(clinical_data_directory) if clinical_data_directory else None
    )
    clinical_data_tsv = Path(clinical_data_tsv) if clinical_data_tsv else None
    convert_images(
        Path(dataset_directory),
        bids_directory,
        clinical_data_directory,
        gif,
        clinical_data_tsv,
    )
    _write_bidsignore(bids_directory)
    cprint("Conversion to BIDS succeeded.")


if __name__ == "__main__":
    cli()
