import click

from .adni_to_bids import adni_to_bids_cli
from .aibl_to_bids import aibl_to_bids_cli
from .habs_to_bids import habs_to_bids_cli
from .nifd_to_bids import nifd_to_bids_cli
from .oasis3_to_bids import oasis3_to_bids_cli
from .oasis_to_bids import oasis_to_bids_cli
from .ukb_to_bids import ukb_to_bids_cli


@click.group("convert")
def cli() -> None:
    """Convert popular neuroimaging datasets to the BIDS format."""
    pass


cli.add_command(adni_to_bids_cli.cli)
cli.add_command(aibl_to_bids_cli.cli)
cli.add_command(habs_to_bids_cli.cli)
cli.add_command(nifd_to_bids_cli.cli)
cli.add_command(oasis_to_bids_cli.cli)
cli.add_command(oasis3_to_bids_cli.cli)
cli.add_command(ukb_to_bids_cli.cli)

if __name__ == "__main__":
    cli()
