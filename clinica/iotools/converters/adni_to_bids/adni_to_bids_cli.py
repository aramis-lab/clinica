from pathlib import Path
from typing import List, Optional

import click

from clinica import option
from clinica.iotools.converters import cli_param

ALL_MODALITIES = (
    "T1",
    "PET_FDG",
    "PET_AMYLOID",
    "PET_TAU",
    "DWI",
    "FLAIR",
    "fMRI",
    "FMAP",
)


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
    type=click.Choice(ALL_MODALITIES),
    default=ALL_MODALITIES,
    help="Convert only the selected modality. By default, all available modalities are converted.",
)
@click.option(
    "-xml", "--xml_path", help="Path to the root directory containing the xml metadata."
)
@option.global_option_group
@option.n_procs
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    xml_path: Optional[str] = None,
    subjects_list: Optional[str] = None,
    clinical_data_only: bool = False,
    force_new_extraction: bool = False,
    modalities: List[str] = ALL_MODALITIES,
    n_procs: Optional[int] = None,
) -> None:
    """ADNI to BIDS converter.

    Convert the imaging and clinical data of ADNI (https://adni.loni.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import convert
    from clinica.utils.exceptions import ClinicaParserError

    if clinical_data_only and force_new_extraction:
        raise ClinicaParserError(
            "Arguments `clinical_data_only` and `force_new_extraction` are mutually exclusive."
        )
    convert(
        Path(dataset_directory),
        Path(bids_directory),
        Path(clinical_data_directory),
        clinical_data_only=clinical_data_only,
        subjects=subjects_list,
        modalities=modalities,
        xml_path=Path(xml_path) if xml_path else None,
        force_new_extraction=force_new_extraction,
        n_procs=n_procs,
    )


if __name__ == "__main__":
    cli()
