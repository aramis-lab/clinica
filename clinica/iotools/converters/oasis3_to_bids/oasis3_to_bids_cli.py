from os import PathLike

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

    Convert the imaging and clinical data of OASIS3 (https://www.oasis-brains.org), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import convert_images
    from clinica.utils.check_dependency import check_dcm2niix
    from clinica.utils.stream import cprint

    convert_images(dataset_directory, bids_directory, clinical_data_directory)

    cprint("Conversion to BIDS succeeded.")


if __name__ == "__main__":
    cli()
