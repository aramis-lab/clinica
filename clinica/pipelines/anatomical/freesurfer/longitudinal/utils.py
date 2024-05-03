"""This module contains functions to deal with longitudinal datasets.

Currently, Clinica Pipeline class and clinica/utils folder can not handle the case
where we need to manipulate longitudinal IDs.

When a new longitudinal pipeline will be developed into Clinica, refactoring will be needed.
"""
from pathlib import Path
from typing import List, Optional, Tuple

__all__ = [
    "grab_image_ids_from_caps_directory",
    "move_subjects_dir_to_source_dir",
    "save_part_sess_long_ids_to_tsv",
]


def save_part_sess_long_ids_to_tsv(
    participant_ids: list[str],
    session_ids: list[str],
    long_ids: list[str],
    out_folder: Path,
    file_name: Optional[str] = None,
):
    """Save participant, session and longitudinal IDs to TSV file.

    TODO: Find a way to merge with utils/save_participants_sessions.py::read_participant_tsv into one util
    """
    import pandas as pd

    from clinica.utils.stream import cprint

    out_folder.mkdir(exist_ok=True)
    tsv_file = out_folder / (file_name or "participants.tsv")

    try:
        data = pd.DataFrame(
            {
                "participant_id": participant_ids,
                "session_id": session_ids,
                "long_id": long_ids,
            }
        )
        data.to_csv(tsv_file, sep="\t", index=False, encoding="utf-8")
    except Exception as e:
        cprint(f"Impossible to save {tsv_file} with pandas")
        raise e


def grab_image_ids_from_caps_directory(
    caps_dir: Path,
) -> Tuple[List[str], List[str], List[str]]:
    """Parse CAPS directory to extract participants, sessions and longitudinal IDs.

    Note:
        This function is a simplified version of create_subs_sess_list for CAPS folders with longitudinal IDs

    TODO: Find a way to merge with iotools/utils/data_handling::create_subs_sess_list into one util

    Use case:
    - CAPS
        - subjects
            - sub-CLNC01
                - long-M000M018
                - ses-M000
                - ses-M018
                - ses-M036
            - sub-CLNC02
                - ses-M000
                - ses-M018

    part_ids = ["sub-CLNC01",  "sub-CLNC01",  "sub-CLNC01" ]
    sess_ids = ["ses-M000",     "ses-M018",     "ses-M036"    ]
    long_ids = ["long-M000M018", "long-M000M018", "long-M000M018"]
    (sub-CLNC02 does not have longitudinal ID so it does not appear on the result)

    Parameters
    ----------
    caps_dir : str
        Path to the CAPS directory.
    """
    participants_paths = sorted([f for f in (caps_dir / "subjects").glob("*sub-*")])
    if len(participants_paths) == 0:
        raise IOError("Dataset empty or not CAPS compliant.")
    participant_ids, session_ids, long_ids = [], [], []
    for part_path in participants_paths:
        for ses_path in part_path.glob("*ses-*"):
            for long_path in part_path.glob("*long-*"):
                participant_ids.append(part_path.name)
                session_ids.append(ses_path.name)
                long_ids.append(long_path.name)

    return participant_ids, session_ids, long_ids


def move_subjects_dir_to_source_dir(
    subjects_dir: Path,
    source_dir: Path,
    subject_id: str,
    image_id: Optional[str] = None,
) -> str:
    """
    Move content of `subjects_dir`/`subject_id` to `source_dir`.

    This function will move content of `subject_id` if recon-all has run in $(TMP). This happens when only
    one time point is used. Content of $(TMP)/`subject_id` is copied to `source_dir` before the deletion of $(TMP).

    Parameters
    ----------
    subjects_dir : Path
        $(TMP), if segmentation was performed on 1 time point,
        <base_dir>/<Pipeline.Name>/ReconAll/`subject_id` otherwise

    source_dir : Path
        <base_dir>/<Pipeline.Name>/ReconAll folder.

    subject_id : str
        Subject ID (e.g. "sub-CLNC01_ses-M000" or "sub-CLNC01_ses-M000M018")

    image_id : str, optional

    Returns
    -------
    str :
        subject_id for node connection with Nipype.
    """
    import shutil

    from clinica.utils.stream import cprint

    if str(source_dir) not in str(subjects_dir):
        shutil.copytree(
            src=subjects_dir / subject_id,
            dst=source_dir / (image_id or subject_id) / subject_id,
            symlinks=True,
        )
        shutil.rmtree(subjects_dir)
        cprint(
            msg=(
                f"Segmentation of {subject_id.replace('_', ' | ')} "
                f"has moved to working directory and $SUBJECTS_DIR folder ({subjects_dir}) was deleted."
            ),
            lvl="warning",
        )

    return subject_id
