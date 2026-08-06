"""Convert OASIS2 dataset (https://sites.wustl.edu/oasisbrains/) to BIDS."""

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

        Workflow
    --------
    1. Read the longitudinal demographics CSV (clinical data).
    2. Discover all raw ANALYZE T1w acquisitions (imaging data).
    3. Merge both data sources and build BIDS metadata tables.
    4. Convert every ANALYZE `.img/.hdr` pair to compressed NIfTI and
       write them into the BIDS hierarchy.
    5. Write `dataset_description.json`, `participants.tsv`,
       and per-subject `*_sessions.tsv` files.
    6. Write modality-agnostic files (README, dataset description).

    Parameters
    ----------
    path_to_dataset:
        Path to the root of the OASIS-2 imaging directory.
        Must contain session folders named ``OAS2_XXXX_MRY/``.
    bids_dir:
        Path to the output BIDS directory (will be created if absent).
    path_to_clinical:
        Path to the directory containing
        ``oasis_longitudinal_demographics.csv``.
    subjects:
        Not yet implemented — all available subjects are converted.
    n_procs:
        Not yet implemented — conversion runs single-threaded.
    """
    from clinica.utils.stream import cprint

    from .._utils import validate_input_path, write_modality_agnostic_files
    from ..factory import get_converter_name
    from ..study_models import StudyName
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
                f"Subject filtering is not yet implemented in {get_converter_name(StudyName.OASIS2)} converter. "
                "All subjects available will be converted."
            ),
            lvl="warning",
        )
    if n_procs != 1:
        cprint(
            f"{get_converter_name(StudyName.OASIS2)} converter does not support multiprocessing yet. n_procs set to 1.",
            lvl="warning",
        )

    cprint("Reading clinical data …", lvl="info")
    df_clinical = read_clinical_data(path_to_clinical)

    cprint("Discovering imaging data …", lvl="info")
    df_imaging = read_imaging_data(path_to_dataset)

    cprint("Merging imaging and clinical data …", lvl="info")
    df_merged, df_subjects = intersect_data(df_imaging, df_clinical)

    cprint("Building BIDS metadata tables …", lvl="info")
    participants, sessions, scans = dataset_to_bids(df_merged, df_subjects)

    cprint(
        f"Converting {len(scans)} T1w acquisitions "
        f"for {len(participants)} subjects across {len(sessions)} sessions …",
        lvl="info",
    )
    write_bids(
        to=bids_dir,
        participants=participants,
        sessions=sessions,
        scans=scans,
        dataset_directory=path_to_dataset,
    )

    readme_data = {
        "link": "https://sites.wustl.edu/oasisbrains/",
        "desc": (
            "OASIS-2: Longitudinal MRI Data in Nondemented and Demented Older Adults. "
            "This dataset consists of a longitudinal collection of 150 subjects aged "
            "60 to 96, scanned on two or more visits separated by at least one year. "
            "For each subject, 3 or 4 individual T1-weighted MRI scans obtained in a "
            "single scan session are included. All subjects are right-handed and include "
            "both men and women. "
            "72 subjects were characterised as nondemented throughout the study. "
            "64 subjects were characterised as demented at their initial visit and "
            "remained so for subsequent scans, including 51 individuals with mild to "
            "moderate Alzheimer's disease. "
            "14 subjects were characterised as nondemented at their initial visit and "
            "subsequently characterised as demented at a later visit."
        ),
    }
    write_modality_agnostic_files(
        study_name=StudyName.OASIS2,
        readme_data=readme_data,
        bids_dir=bids_dir,
    )
    cprint("Conversion to BIDS succeeded.", lvl="info")


# todo : mention Nikhil
