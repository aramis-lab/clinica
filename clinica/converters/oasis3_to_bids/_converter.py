"""Convert the OASIS3 dataset into BIDS."""

from typing import Optional

from clinica.utils.filemanip import UserProvidedPath

__all__ = ["convert"]


def convert(
    path_to_dataset: UserProvidedPath,
    bids_dir: UserProvidedPath,
    path_to_clinical: UserProvidedPath,
    subjects: Optional[UserProvidedPath] = None,
    n_procs: Optional[int] = 1,
    **kwargs,
):
    """Convert the entire dataset in BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    from clinica.converters.factory import get_converter_name
    from clinica.converters.study_models import StudyName
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path, write_modality_agnostic_files
    from ._utils import (
        dataset_to_bids,
        intersect_data,
        read_clinical_data,
        read_imaging_data,
        write_bids,
    )

    path_to_dataset = validate_input_path(path_to_dataset)
    bids_dir = validate_input_path(bids_dir, check_exist=False)
    path_to_clinical = validate_input_path(path_to_clinical)
    if subjects:
        cprint(
            (
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.OASIS3)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.OASIS3)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )
    dict_df = read_clinical_data(path_to_clinical)
    imaging_data = read_imaging_data(path_to_dataset)
    imaging_data, df_small = intersect_data(imaging_data, dict_df)
    participants, sessions, scans = dataset_to_bids(imaging_data, df_small)
    write_bids(
        to=bids_dir,
        participants=participants,
        sessions=sessions,
        scans=scans,
        dataset_directory=path_to_dataset,
    )
    readme_data = {
        "link": "https://sites.wustl.edu/oasisbrains/#access",
        "desc": (
            "OASIS-3 is a retrospective compilation of data for 1378 participants that were collected across several "
            "ongoing projects through the WUSTL Knight ADRC over the course of 30years. Participants include 755 "
            "cognitively normal adults and 622 individuals at various stages of cognitive decline ranging in age from "
            "42-95yrs. All participants were assigned a new random identifier and all dates were removed and "
            "normalized to reflect days from entry into study. The dataset contains 2842 MR sessions which include "
            "T1w, T2w, FLAIR, ASL, SWI, time of flight, resting-state BOLD, and DTI sequences. Many of the MR sessions "
            "are accompanied by volumetric segmentation files produced through FreeSurfer processing. PET imaging from "
            "different tracers, PIB, AV45, and FDG, totaling over 2157 raw imaging scans and the accompanying "
            "post-processed files from the Pet Unified Pipeline (PUP) are also available in OASIS-3."
        ),
    }
    write_modality_agnostic_files(
        study_name=StudyName.OASIS3,
        readme_data=readme_data,
        bids_dir=bids_dir,
    )
    cprint("Conversion to BIDS succeeded.", lvl="info")
