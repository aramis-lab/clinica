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
    long_id = "long-" + "".join(list_session_label)

    return long_id


def get_participants_long_id(list_participant_id, list_session_id):
    """Extract list of longitudinal IDs from a set of participant and session IDs.

    Example:
        >>> from clinica.utils.longitudinal import get_participants_long_id
        >>> get_participants_long_id(['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02'], ['ses-M00', 'ses-M18', 'ses-M00'])
        ['long-M00M18', 'long-M00M18', 'long-M00']
    """
    from .participant import get_unique_subjects

    unique_subject_list, per_subject_session_list = get_unique_subjects(
        list_participant_id, list_session_id
    )

    list_long_id = []
    for i in range(0, len(unique_subject_list)):
        list_long_id = list_long_id + [get_long_id(per_subject_session_list[i])] * len(
            per_subject_session_list[i]
        )

    return list_long_id


def save_long_id(list_session_id, output_dir, file_name=None):
    """Save long ID to `file_name`."""
    import os

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    long_id = get_long_id(list_session_id)
    file_name = file_name or f"{long_id}_sessions.tsv"
    sessions_tsv = open(os.path.join(output_dir, file_name), "w")
    sessions_tsv.write("session_id\n")

    for session_id in sorted(list_session_id):
        sessions_tsv.write(session_id + "\n")

    sessions_tsv.close()


def read_sessions(caps_dir, participant_id, long_id):
    """Extract sessions IDs from `caps_dir`/subjects/`participant_id`/`long_id`/`long_id`_sessions.tsv."""
    import os

    import pandas

    from clinica.utils.exceptions import ClinicaException

    sessions_file = os.path.join(
        os.path.expanduser(caps_dir),
        "subjects",
        participant_id,
        long_id,
        f"{long_id}_sessions.tsv",
    )
    if not os.path.isfile(sessions_file):
        raise ClinicaException(
            "The TSV file with sessions associated "
            f"to {participant_id} for longitudinal ID {long_id} is missing "
            f"(expected path: {sessions_file})."
        )
    ss_df = pandas.read_csv(sessions_file, sep="\t")
    if "session_id" not in list(ss_df.columns.values):
        raise ClinicaException(
            "The TSV file does not contain session_id column "
            f"(path: {sessions_file})."
        )

    return list(ss_df.session_id)
