from pathlib import Path

import click

from clinica.iotools.converters import cli_param


@click.command(name="nifd-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: Path,
    clinical_data_directory: Path,
    bids_directory: Path,
) -> None:
    """NIFD to BIDS converter.

    Convert the imaging and clinical data of NIFD (https://4rtni-ftldni.ini.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from .nifd_to_bids import convert

    convert(dataset_directory, bids_directory, clinical_data_directory)


if __name__ == "__main__":
    cli()
