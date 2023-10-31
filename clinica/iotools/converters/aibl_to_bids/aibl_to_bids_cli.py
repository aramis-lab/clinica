from typing import Optional

import click

from clinica import option
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
@option.global_option_group
@option.n_procs
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    clinical_data_only: bool = False,
    overwrite: bool = False,
    n_procs: Optional[int] = None,
) -> None:
    """AIBL to BIDS converter.

    Convert the imaging and clinical data of AIBL (https://aibl.csiro.au/adni/index.html), located in DATASET_DIRECTORY
    and CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from pathlib import Path

    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert

    convert(
        Path(dataset_directory),
        Path(clinical_data_directory),
        Path(bids_directory),
        overwrite,
        clinical_data_only,
        n_procs=n_procs,
    )
