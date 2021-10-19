import click

from clinica.iotools.converters import cli_param


@click.command(name="aibl-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
@cli_param.clinical_data_only
@click.option(
    "--overwrite",
    is_flag=True,
    help="Overwrites previously written nifti and json files.",
)
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    clinical_data_only: bool = False,
    overwrite: bool = False,
) -> None:
    """AIBL to BIDS converter.

    Convert the imaging and clinical data of AIBL (https://aibl.csiro.au/adni/index.html), located in DATASET_DIRECTORY
    and CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from os import makedirs

    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import (
        convert_clinical_data,
        convert_images,
    )
    from clinica.utils.check_dependency import check_dcm2niix

    check_dcm2niix()

    makedirs(bids_directory, exist_ok=True)

    if not clinical_data_only:
        convert_images(
            dataset_directory,
            clinical_data_directory,
            bids_directory,
            overwrite,
        )

    convert_clinical_data(bids_directory, clinical_data_directory)
