from typing import Optional

import click

from clinica import option
from clinica.converters import cli_param


@click.command(name="oasis-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
@cli_param.subjects_list
@option.global_option_group
@option.n_procs
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    subjects_list: Optional[str] = None,
    n_procs: Optional[int] = None,
) -> None:
    """OASIS to BIDS converter.

    Convert the imaging and clinical data of OASIS (https://sites.wustl.edu/oasisbrains/),
    located in DATASET_DIRECTORY and CLINICAL_DATA_DIRECTORY respectively,
    to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from ._converter import convert

    convert(
        dataset_directory,
        bids_directory,
        clinical_data_directory,
        subjects=subjects_list,
        n_procs=n_procs,
    )


if __name__ == "__main__":
    cli()
