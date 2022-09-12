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
    readme_dict = {
        "link": "https://habs.mgh.harvard.edu",
        "desc": "The overall goal of the Harvard Aging Brain Study (HABS) is to elucidate the earliest changes in molecular, functional and structural imaging markers that signal the transition from normal cognition to progressive cognitive decline along the trajectory of preclinical Alzheimer’s Disease.",
    }
    bids.write_modality_agnostic_files(
        study_name="HABS", data_dict=readme_dict, bids_dir=rawdata
    )


if __name__ == "__main__":
    cli()
