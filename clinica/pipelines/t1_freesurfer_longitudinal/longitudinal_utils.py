"""This module contains functions to deal with longitudinal datasets.

Currently, Clinica Pipeline class and clinica/utils folder can not handle the case
where we need to manipulate longitudinal IDs.

When a new longitudinal pipeline will be developed into Clinica, refactoring will be needed.
"""


def extract_subject_session_longitudinal_ids_from_filename(bids_or_caps_files):
    """Extract participant/session/longitudinal IDs from filename.

    Example:
        ['sub-CLNC01', 'sub-CLNC01'], ['ses-M000', 'ses-M018'], ['long-M000M018', 'long-M000M018'])

    TODO: Find a way to merge with utils/filemanip.py::extract_subjects_sessions_from_filename into one util
    """
    import re

    id_bids_or_caps_files = [
        re.search(
            r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_(long-[a-zA-Z0-9]+)", file
        ).group()
        for file in bids_or_caps_files
    ]
    split = [image_id.split("_") for image_id in id_bids_or_caps_files]
    part_ids = [p_id[0] for p_id in split]
    sess_ids = [s_id[1] for s_id in split]
    long_ids = [l_id[2] for l_id in split]
    return part_ids, sess_ids, long_ids


def read_part_sess_long_ids_from_tsv(tsv_file):
    """Extract participant, session and longitudinal from TSV file.

    TODO: Find a way to merge with utils/filemanip.py::read_participant_tsv into one util
    """
    import os

    import pandas

    from clinica.utils.exceptions import ClinicaException

    if not os.path.isfile(tsv_file):
        raise ClinicaException(
            "The TSV file you gave is not a file.\n"
            "Error explanations:\n"
            f" - Clinica expected the following path to be a file: {tsv_file}\n"
            " - If you gave relative path, did you run Clinica on the good folder?"
        )
    df = pandas.read_csv(tsv_file, sep="\t")

    def check_key_in_data_frame(file, data_frame, key):
        if key not in list(data_frame.columns.values):
            raise ClinicaException(
                f"The TSV file does not contain {key} column (path: {file})"
            )

    check_key_in_data_frame(tsv_file, df, "participant_id")
    check_key_in_data_frame(tsv_file, df, "session_id")
    check_key_in_data_frame(tsv_file, df, "long_id")

    participants = list(df.participant_id)
    sessions = list(df.session_id)
    longs = list(df.long_id)

    # Remove potential whitespace in participant, session or longitudinal ID
    return (
        [sub.strip(" ") for sub in participants],
        [ses.strip(" ") for ses in sessions],
        [lng.strip(" ") for lng in longs],
    )


def save_part_sess_long_ids_to_tsv(
    participant_ids, session_ids, long_ids, out_folder, file_name=None
):
    """Save participant, session and longitudinal IDs to TSV file.

    TODO: Find a way to merge with utils/save_participants_sessions.py::read_participant_tsv into one util
    """
    import os

    import pandas

    from clinica.utils.stream import cprint

    os.makedirs(out_folder, exist_ok=True)

    if file_name:
        tsv_file = os.path.join(out_folder, file_name)
    else:
        tsv_file = os.path.join(out_folder, "participants.tsv")

    try:
        data = pandas.DataFrame(
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


def extract_participant_long_ids_from_filename(caps_files):
    """TODO: Find a way to merge with utils/filemanip.py::extract_subjects_sessions_from_filename into one util."""
    import re

    caps_files = [
        re.search(r"(sub-[a-zA-Z0-9]+)_(long-[a-zA-Z0-9]+)", file).group()
        for file in caps_files
    ]
    split = [image_id.split("_") for image_id in caps_files]
    part_ids = [p_id[0] for p_id in split]
    long_ids = [l_id[1] for l_id in split]
    return part_ids, long_ids


def grab_image_ids_from_caps_directory(caps_dir):
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

    Args:
        caps_dir (str): Path to the CAPS directory.
    """
    import os
    from glob import glob

    participants_paths = glob(os.path.join(caps_dir, "subjects", "*sub-*"))

    participants_paths.sort()
    if len(participants_paths) == 0:
        raise IOError("Dataset empty or not CAPS compliant.")

    part_ids = []
    sess_ids = []
    long_ids = []
    for part_path in participants_paths:
        part_id = part_path.split(os.sep)[-1]
        session_paths = glob(os.path.join(part_path, "*ses-*"))
        longitudinal_paths = glob(os.path.join(part_path, "*long-*"))
        if len(session_paths) and len(longitudinal_paths):
            for ses_path in session_paths:
                sess_id = ses_path.split(os.sep)[-1]
                for long_path in longitudinal_paths:
                    long_id = long_path.split(os.sep)[-1]

                    part_ids.append(part_id)
                    sess_ids.append(sess_id)
                    long_ids.append(long_id)

    return part_ids, sess_ids, long_ids
