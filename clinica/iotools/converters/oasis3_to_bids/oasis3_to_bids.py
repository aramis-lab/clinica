"""Convert the NIFD dataset into BIDS."""

from pathlib import Path

__all__ = ["convert_images"]


def convert_images(
    path_to_dataset: Path,
    bids_dir: Path,
    path_to_clinical: Path,
) -> list[str]:
    """Convert the entire dataset in BIDS.

    Scans available files in the path_to_dataset,
    identifies the patients that have images described by the JSON file,
    converts the image with the highest quality for each category.
    """
    import clinica.iotools.bids_utils as bids

    from .oasis3_utils import (
        dataset_to_bids,
        intersect_data,
        read_clinical_data,
        read_imaging_data,
        write_bids,
    )

    dict_df = read_clinical_data(path_to_clinical)
    imaging_data = read_imaging_data(path_to_dataset)
    imaging_data, df_small = intersect_data(imaging_data, dict_df)
    participants, sessions, scans = dataset_to_bids(imaging_data, df_small)
    written = write_bids(
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
    bids.write_modality_agnostic_files(
        study_name=bids.StudyName.OASIS3,
        readme_data=readme_data,
        bids_dir=bids_dir,
    )

    return written
