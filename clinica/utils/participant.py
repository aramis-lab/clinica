"""This module contains utilities for longitudinal pipelines.

See CAPS specifications for details about long ID.
"""
from os import PathLike
from typing import List, Optional, Tuple

from clinica.utils.filemanip import read_participant_tsv


def get_unique_subjects(
    subjects: List[str], sessions: List[str]
) -> Tuple[List[str], List[List[str]]]:
    """Get unique participant IDs with their sessions.

    This function generates a list of unique participant IDs from `in_subject_list` with their sessions.

    Parameters
    ----------
    subjects : list of str
        List of participant IDs.

    sessions : list of str
        List of session IDs.

    Returns
    -------
    unique_subjects : list of str
        List of participant IDs, where each participant appears only once.

    sessions_per_subject : list of lists of str
        List of session IDs associated to any single participant.

    Examples
    --------
    >>> from clinica.utils.participant import get_unique_subjects
    >>> get_unique_subjects(['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02'], ['ses-M000', 'ses-M018', 'ses-M000'])
    (['sub-CLNC01', 'sub-CLNC02'], [['ses-M000', 'ses-M018'], ['ses-M000']])
    """
    import numpy as np

    if len(subjects) != len(sessions):
        raise ValueError(
            "The number of subjects should match the number of sessions.\n"
            f"You provided the following subjects: {subjects}.\n"
            f"And the following sessions: {sessions}."
        )
    subjects = np.array(subjects)
    sessions = np.array(sessions)

    # The second returned element indicates for each participant_id the
    # element they correspond to in the 'unique' list. We will use this
    # to link each session_id in the repeated list of session_id to
    # their corresponding unique participant_id

    unique_subjects, inverse_positions = np.unique(subjects, return_inverse=True)
    sessions_per_subject = [
        sessions[inverse_positions == subject_index].tolist()
        for subject_index in range(len(unique_subjects))
    ]
    if len(unique_subjects) != len(sessions_per_subject):
        raise ValueError("Problem while getting unique subjects and sessions lists.")

    return unique_subjects.tolist(), sessions_per_subject


def unique_subjects_sessions_to_subjects_sessions(
    subjects: List[str],
    sessions_per_subject: List[List[str]],
) -> Tuple[List[str], List[str]]:
    """Do reverse operation of get_unique_subjects function.

    Parameters
    ----------
    subjects : list of str
        List of unique subject identifiers.

    sessions_per_subject : list of lists of str
        The sessions for each subject in subjects.

    Returns
    -------
    participants : list of str
        The list of subjects of same length as sessions.

    sessions : list of str
        The list of sessions of same length as participants.

    Examples
    --------
    >>> from clinica.utils.participant import unique_subjects_sessions_to_subjects_sessions
    >>> unique_subjects_sessions_to_subjects_sessions(['sub-01', 'sub-02'], [['ses-M000', 'ses-M018'], ['ses-M000']])
    (['sub-01', 'sub-01', 'sub-02'], ['ses-M000', 'ses-M018', 'ses-M000'])
    """
    if len(subjects) != len(sessions_per_subject):
        raise ValueError(
            "The number of unique subjects should match the number of session lists.\n"
            f"You provided the following subjects: {subjects}.\n"
            f"And the following session lists: {sessions_per_subject}."
        )
    participants, sessions = [], []
    for idx, participant_id in enumerate(subjects):
        for session_id in sessions_per_subject[idx]:
            participants.append(participant_id)
            sessions.append(session_id)

    return participants, sessions


def get_subject_session_list(
    input_dir: PathLike,
    subject_session_file: Optional[PathLike] = None,
    is_bids_dir: bool = True,
    use_session_tsv: bool = False,
    tsv_dir: Optional[PathLike] = None,
) -> Tuple[List[str], List[str]]:
    """Parse a BIDS or CAPS directory to get the subjects and sessions.

    This function lists all the subjects and sessions based on the content of
    the BIDS or CAPS directory or (if specified) on the provided
    subject-sessions TSV file.

    Parameters
    ----------
    input_dir : PathLike
        A BIDS or CAPS directory path.

    subject_session_file : PathLike, optional
        A subjects-sessions file in TSV format.

    is_bids_dir : bool, optional
        Indicates if input_dir is a BIDS or CAPS directory.
        Default=True.

    use_session_tsv : bool, optional
        Specify if the list uses the sessions listed in the sessions.tsv files.
        Default=False.

    tsv_dir : PathLike, optional
        If TSV file does not exist, it will be created in output_dir.
        If not specified, output_dir will be in <tmp> folder

    Returns
    -------
    subjects : list of str
        Subjects list.

    sessions : list of str
        Sessions list.

    Notes
    -----
    This is a generic method based on folder names. If your <BIDS> dataset contains e.g.:
        - sub-CLNC01/ses-M000/anat/sub-CLNC01_ses-M000_T1w.nii
        - sub-CLNC02/ses-M000/dwi/sub-CLNC02_ses-M000_dwi.{bval|bvec|json|nii}
        - sub-CLNC02/ses-M000/anat/sub-CLNC02_ses-M000_T1w.nii
        get_subject_session_list(<BIDS>, None, True) will return
        ['ses-M000', 'ses-M000'], ['sub-CLNC01', 'sub-CLNC02'].

    However, if your pipeline needs both T1w and DWI files, you will need to check
    with e.g. clinica_file_reader_function.
    """
    import tempfile
    from pathlib import Path
    from time import localtime, strftime, time

    from clinica.iotools.utils.data_handling import create_subs_sess_list

    if not subject_session_file:
        output_dir = Path(tsv_dir) if tsv_dir else Path(tempfile.mkdtemp())
        timestamp = strftime("%Y%m%d_%H%M%S", localtime(time()))
        tsv_file = f"subjects_sessions_list_{timestamp}.tsv"
        subject_session_file = output_dir / tsv_file
        create_subs_sess_list(
            input_dir=input_dir,
            output_dir=output_dir,
            file_name=tsv_file,
            is_bids_dir=is_bids_dir,
            use_session_tsv=use_session_tsv,
        )

    return read_participant_tsv(subject_session_file)
