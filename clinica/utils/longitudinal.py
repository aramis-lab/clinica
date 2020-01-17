# coding: utf8


"""This module contains utilities for longitudinal pipelines. See CAPS specifications for details about long ID."""


def get_long_id(list_session_id):
    """Extract longitudinal ID from a set of session IDs.

    This will create a unique identifier for a participant and its corresponding sessions. Sessions labels are sorted
    alphabetically before being merged in order to generate the longitudinal ID.

    Args:
        list_session_id (list[str]): List of session IDs
            (e.g. ["ses-M00"] or ["ses-M00", "ses-M18", "ses-M36"])

    Returns:
        Longitudinal ID (str)

    Example:
        >>> from clinica.utils.longitudinal import get_long_id
        >>> get_long_id(['ses-M00'])
        'long-M00'
        >>> get_long_id(['ses-M00', 'ses-M18', 'ses-M36'])
        'long-M00M18M36'
        >>> get_long_id(['ses-M18', 'ses-M36', 'ses-M00']) # Session IDs do not need to be sorted
        'long-M00M18M36'
    """
    sorted_list = sorted(list_session_id)
    list_session_label = [session_id[4:] for session_id in sorted_list]
    long_id = 'long-' + ''.join(list_session_label)

    return long_id


def get_participants_long_id(list_participant_id, list_session_id):
    """Extract list of longitudinal IDs from a set of participant and session IDs."""
    from .participant import get_unique_subjects

    unique_subject_list, per_subject_session_list = get_unique_subjects(list_participant_id, list_session_id)

    list_long_id = []
    for i in range(0, len(unique_subject_list)):
        for s_id in per_subject_session_list[i]:
            list_long_id.append(get_long_id(per_subject_session_list[i]))

    return list_long_id


def save_long_id(list_session_id, output_dir, file_name=None):
    """Save long ID to `file_name`."""
    import os

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    long_id = get_long_id(list_session_id)
    if file_name is None:
        file_name = long_id + '_sessions.tsv'
    sessions_tsv = open(os.path.join(output_dir, file_name), 'w')
    sessions_tsv.write('session_id\n')

    for session_id in sorted(list_session_id):
        sessions_tsv.write(session_id + '\n')

    sessions_tsv.close()


def extract_session_ids(tsv_sessions):
    """Extract sessions IDs from TSV file.

    Raise:
        ClinicaException if participant_id or session_id column is missing from TSV file
    """
    import pandas as pd
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaException

    ss_df = pd.io.parsers.read_csv(tsv_sessions, sep='\t')
    if 'session_id' not in list(ss_df.columns.values):
        raise ClinicaException(
            "\n%s[Error] The TSV file does not contain session_id column (path: %s)%s" %
            (Fore.RED, tsv_sessions, Fore.RESET)
        )

    return list(ss_df.session_id)
