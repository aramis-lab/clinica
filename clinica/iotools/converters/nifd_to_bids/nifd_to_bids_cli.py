from pathlib import PurePath

import click

from clinica.iotools.converters import cli_param


@click.command(name="nifd-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: PurePath,
    clinical_data_directory: PurePath,
    bids_directory: PurePath,
) -> None:
    """NIFD to BIDS converter.

    Convert the imaging and clinical data of NIFD (https://4rtni-ftldni.ini.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images
    from clinica.utils.check_dependency import check_dcm2niix
    from clinica.utils.stream import cprint

    check_dcm2niix()

    convert_images(dataset_directory, bids_directory, clinical_data_directory)

    cprint("Conversion to BIDS succeeded.")


if __name__ == "__main__":
    cli()
