"""Convert the UKB dataset into BIDS."""

from pathlib import Path

__all__ = ["convert"]


def convert(
    path_to_dataset: Path,
    bids_dir: Path,
    path_to_clinical: Path,
    **kwargs,
):
    """Convert the entire dataset to BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    from clinica.iotools.bids_utils import StudyName, write_modality_agnostic_files
    from clinica.utils.check_dependency import ThirdPartySoftware, check_software
    from clinica.utils.stream import cprint

    from .ukb_utils import (
        find_clinical_data,
        merge_imaging_and_clinical_data,
        prepare_dataset_to_bids_format,
        read_imaging_data,
        write_bids,
    )

    check_software(ThirdPartySoftware.DCM2NIIX)
    result = prepare_dataset_to_bids_format(
        merge_imaging_and_clinical_data(
            read_imaging_data(path_to_dataset),
            find_clinical_data(path_to_clinical),
        )
    )
    write_bids(
        to=bids_dir,
        participants=result["participants"],
        sessions=result["sessions"],
        scans=result["scans"],
        dataset_directory=path_to_dataset,
    )
    readme_data = {
        "link": "https://www.ukbiobank.ac.uk/",
        "desc": (
            "UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and "
            "health information from half a million UK participants. The database is regularly augmented with "
            "additional data and is globally accessible to approved researchers undertaking vital research into the "
            "most common and life-threatening diseases. It is a major contributor to the advancement of modern "
            "medicine it and has led to the discovery of several scientific advances and numerous treatments to "
            "improve human health."
        ),
    }
    write_modality_agnostic_files(
        study_name=StudyName.UKB,
        readme_data=readme_data,
        bids_dir=bids_dir,
    )
    cprint("Conversion to BIDS succeeded.")
