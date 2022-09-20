"""This module contains utilities for longitudinal pipelines.

See CAPS specifications for details about long ID.
"""

from clinica.utils.filemanip import read_participant_tsv


def get_unique_subjects(in_subject_list, in_session_list):
    """Get unique participant IDs with their sessions.

    This function generates a list of unique participant IDs from `in_subject_list` with their sessions.
    Args:
        in_subject_list (list[str]): list of participant IDs
        in_session_list (list[str]): list of session IDs

    Returns:
        out_unique_subject_list (list[str]): list of participant IDs, where each participant appears only once
        out_per_subject_session_list (list[list[str]]): list of list
            (list2) of session_id associated to any single participant

    Example:
        >>> from clinica.utils.participant import get_unique_subjects
        >>> get_unique_subjects(['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02'], ['ses-M000', 'ses-M018', 'ses-M000'])
        (['sub-CLNC01', 'sub-CLNC02'], [['ses-M000', 'ses-M018'], ['ses-M000']])
    """
    import numpy as np

    subject_array = np.array(in_subject_list)
    session_array = np.array(in_session_list)

    # The second returned element indicates for each participant_id the
    # element they correspond to in the 'unique' list. We will use this
    # to link each session_id in the repeated list of session_id to
    # their corresponding unique participant_id

    unique_subject_array, out_inverse_positions = np.unique(
        subject_array, return_inverse=True
    )
    out_unique_subject_list = unique_subject_array.tolist()

    subject_number = len(out_unique_subject_list)
    out_per_subject_session_list = [
        session_array[out_inverse_positions == subject_index].tolist()
        for subject_index in range(subject_number)
    ]

    assert len(out_unique_subject_list) == len(
        out_per_subject_session_list
    ), "Problem while getting unique subjects and sessions lists"

    return out_unique_subject_list, out_per_subject_session_list


def unique_subjects_sessions_to_subjects_sessions(
    unique_subject_list, per_subject_session_list
):
    """Do reverse operation of get_unique_subjects function.

    Example:
        >>> from clinica.utils.participant import unique_subjects_sessions_to_subjects_sessions
        >>> unique_subjects_sessions_to_subjects_sessions(['sub-01', 'sub-02'], [['ses-M000', 'ses-M018'], ['ses-M000']])
        (['sub-CLNC01', 'sub-01', 'sub-02'], ['ses-M000', 'ses-M018', 'ses-M000'])

    """
    list_participants = []
    list_sessions = []
    for idx, participant_id in enumerate(unique_subject_list):
        for session_id in per_subject_session_list[idx]:
            list_participants.append(participant_id)
            list_sessions.append(session_id)

    return list_participants, list_sessions


def get_subject_session_list(
    input_dir, ss_file=None, is_bids_dir=True, use_session_tsv=False, tsv_dir=None
):
    """Parse a BIDS or CAPS directory to get the subjects and sessions.

    This function lists all the subjects and sessions based on the content of
    the BIDS or CAPS directory or (if specified) on the provided
    subject-sessions TSV file.

    Args:
        input_dir: A BIDS or CAPS directory path.
        ss_file: A subjects-sessions file (.tsv format).
        is_bids_dir: Indicates if input_dir is a BIDS or CAPS directory
        use_session_tsv (boolean): Specify if the list uses the sessions listed in the sessions.tsv files
        tsv_dir (str): if TSV file does not exist, it will be created in output_dir. If
            not specified, output_dir will be in <tmp> folder

    Returns:
        subjects: A subjects list.
        sessions: A sessions list.

    Notes:
        This is a generic method based on folder names. If your <BIDS> dataset contains e.g.:
        - sub-CLNC01/ses-M000/anat/sub-CLNC01_ses-M000_T1w.nii
        - sub-CLNC02/ses-M000/dwi/sub-CLNC02_ses-M000_dwi.{bval|bvec|json|nii}
        - sub-CLNC02/ses-M000/anat/sub-CLNC02_ses-M000_T1w.nii
        get_subject_session_list(<BIDS>, None, True) will return
        ['ses-M000', 'ses-M000'], ['sub-CLNC01', 'sub-CLNC02'].

        However, if your pipeline needs both T1w and DWI files, you will need to check
        with e.g. clinica_file_reader_function.
    """
    import os
    import tempfile
    from time import localtime, strftime, time

    import clinica.iotools.utils.data_handling as cdh

    if not ss_file:
        if tsv_dir:
            output_dir = tsv_dir
        else:
            output_dir = tempfile.mkdtemp()
        timestamp = strftime("%Y%m%d_%H%M%S", localtime(time()))
        tsv_file = f"subjects_sessions_list_{timestamp}.tsv"
        ss_file = os.path.join(output_dir, tsv_file)

        cdh.create_subs_sess_list(
            input_dir=input_dir,
            output_dir=output_dir,
            file_name=tsv_file,
            is_bids_dir=is_bids_dir,
            use_session_tsv=use_session_tsv,
        )

    participant_ids, session_ids = read_participant_tsv(ss_file)
    return session_ids, participant_ids
