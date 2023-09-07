from typing import List, Optional

import click

from clinica.iotools.converters import cli_param
from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBidsConverter


@click.command(name="adni-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
@cli_param.subjects_list
@click.option(
    "-c",
    "--clinical-data-only",
    is_flag=True,
    help="Convert clinical data only.",
)
@click.option(
    "-f",
    "--force-new-extraction",
    is_flag=True,
    help="Force new data extraction.",
)
@click.option(
    "-m",
    "--modalities",
    multiple=True,
    type=click.Choice(AdniToBidsConverter.get_modalities_supported()),
    default=AdniToBidsConverter.get_modalities_supported(),
    help="Convert only the selected modality. By default, all available modalities are converted.",
)
@click.option(
    "-xml", "--xml_path", help="Path to the root directory containing the xml metadata."
)
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    xml_path: Optional[str] = None,
    subjects_list: Optional[str] = None,
    clinical_data_only: bool = False,
    force_new_extraction: bool = False,
    modalities: Optional[List[str]] = None,
) -> None:
    """ADNI to BIDS converter.

    Convert the imaging and clinical data of ADNI (https://adni.loni.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    AdniToBidsConverter(
        dataset_directory,
        bids_directory,
        clinical_data_directory,
        xml_path,
        clinical_data_only,
        force_new_extraction,
    ).convert(subjects_list, modalities)


if __name__ == "__main__":
    cli()
