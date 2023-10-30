from typing import Optional

import click

from clinica import option
from clinica.iotools.converters import cli_param


@click.command(name="oasis-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
@option.global_option_group
@option.n_procs
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    n_procs: Optional[int] = None,
) -> None:
    """OASIS to BIDS converter.

    Convert the imaging and clinical data of OASIS (http://oasis-brains.org/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

    oasis_to_bids = OasisToBids()
    oasis_to_bids.convert_images(dataset_directory, bids_directory, n_procs=n_procs)
    oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)


if __name__ == "__main__":
    cli()
