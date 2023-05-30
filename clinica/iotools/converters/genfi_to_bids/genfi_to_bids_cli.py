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
    from clinica.iotools.converters.genfi_to_bids.genfi_to_bids import convert_images
    from clinica.utils.check_dependency import check_dcm2niix
    from clinica.utils.stream import cprint

    check_dcm2niix()

    convert_images(dataset_directory, bids_directory, clinical_data_directory, gif)
    _write_bidsignore(str(bids_directory))

    cprint("Conversion to BIDS succeeded.")


if __name__ == "__main__":
    cli()
