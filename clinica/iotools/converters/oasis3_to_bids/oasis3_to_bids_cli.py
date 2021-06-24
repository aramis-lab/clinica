import click

from clinica.iotools.converters import cli_param


@click.command(name="oasis3-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
) -> None:
    """OASIS-3 to BIDS converter.

    Convert the imaging and clinical data of OASIS-3 (http://oasis-brains.org/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import Oasis3ToBids

    oasis3_to_bids = Oasis3ToBids()
    oasis3_to_bids.convert_images(dataset_directory, bids_directory)
    oasis3_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)


if __name__ == "__main__":
    cli()
