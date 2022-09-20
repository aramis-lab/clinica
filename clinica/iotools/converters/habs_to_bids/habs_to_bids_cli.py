from pathlib import Path

import click

import clinica.iotools.bids_utils as bids


@click.command(name="habs-to-bids")
@click.argument(
    "sourcedata",
    type=click.Path(exists=True, resolve_path=True, path_type=Path),
)
@click.argument(
    "rawdata",
    type=click.Path(writable=True, resolve_path=True, path_type=Path),
)
def cli(sourcedata: str, rawdata: str) -> None:
    """HABS to BIDS converter."""
    from pandas import concat

    from .habs_to_bids import (
        find_clinical_data,
        find_imaging_data,
        parse_imaging_data,
        read_clinical_data,
        write_bids,
    )

    clinical_data = {
        k: read_clinical_data(sourcedata / p, c)
        for k, p, c in find_clinical_data(sourcedata)
    }

    imaging_data = concat(
        [parse_imaging_data(x) for x in find_imaging_data(sourcedata)]
    )

    write_bids(
        sourcedata=sourcedata,
        rawdata=rawdata,
        imaging_data=imaging_data,
        clinical_data=clinical_data,
    )


if __name__ == "__main__":
    cli()
