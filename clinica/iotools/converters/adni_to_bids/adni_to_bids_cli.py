from typing import List, Optional

import click

from clinica.iotools.converters import cli_param

ALL_MODALITIES = ("T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI")


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
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
    xml_path: Optional[str] = None,
    subjects_list: Optional[str] = None,
    clinical_data_only: bool = False,
    force_new_extraction: bool = False,
    modalities: List[str] = ALL_MODALITIES,
) -> None:
    """ADNI to BIDS converter.

    Convert the imaging and clinical data of ADNI (https://adni.loni.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
    from clinica.utils.exceptions import ClinicaParserError

    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()

    if clinical_data_only and force_new_extraction:
        raise ClinicaParserError(
            "Arguments `clinical_data_only` and `force_new_extraction` are mutually exclusive."
        )

    if not clinical_data_only:
        adni_to_bids.convert_images(
            dataset_directory,
            clinical_data_directory,
            bids_directory,
            subjects_list,
            modalities,
            force_new_extraction,
        )

    adni_to_bids.convert_clinical_data(
        clinical_data_dir=clinical_data_directory,
        out_path=bids_directory,
        clinical_data_only=clinical_data_only,
        subjects_list_path=subjects_list,
        xml_path=xml_path,
    )


if __name__ == "__main__":
    cli()
