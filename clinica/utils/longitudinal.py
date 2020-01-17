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


def grab_image_ids_from_caps_directory(caps_dir):
    """
    Parse CAPS directory to extract participants, sessions and longitudinal IDs.

    Use case:

    - CAPS
        - subjects
            - sub-CLNC01
                - long-M00M18
                - ses-M00
                - ses-M18
                - ses-M36
            - sub-CLNC02
                - ses-M00
                - ses-M18

    part_ids = ["sub-CLNC01",  "sub-CLNC01",  "sub-CLNC01" ]
    sess_ids = ["ses-M00",     "ses-M18",     "ses-M36"    ]
    long_ids = ["long-M00M18", "long-M00M18", "long-M00M18"]

    Args:
        caps_dir (str): Path to the CAPS directory.
    """
    from glob import glob
    import os

    participants_paths = glob(os.path.join(caps_dir, 'subjects', '*sub-*'))

    participants_paths.sort()
    if len(participants_paths) == 0:
        raise IOError('Dataset empty or not CAPS compliant.')

    part_ids = []
    sess_ids = []
    long_ids = []
    for part_path in participants_paths:
        part_id = part_path.split(os.sep)[-1]
        session_paths = glob(os.path.join(part_path, '*ses-*'))
        longitudinal_paths = glob(os.path.join(part_path, '*long-*'))
        if len(session_paths) and len(longitudinal_paths):
            for ses_path in session_paths:
                sess_id = ses_path.split(os.sep)[-1]
                for long_path in longitudinal_paths:
                    long_id = long_path.split(os.sep)[-1]

                    part_ids.append(part_id)
                    sess_ids.append(sess_id)
                    long_ids.append(long_id)

    return part_ids, sess_ids, long_ids


def read_participant_tsv(tsv_file):
    """Extract participant IDs and session IDs from TSV file.

    Raise:
        ClinicaException if tsv_file is not a file
        ClinicaException if participant_id or session_id column is missing from TSV file
    """
    import os
    import pandas as pd
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaException

    if not os.path.isfile(tsv_file):
        raise ClinicaException(
            "\n%s[Error] The TSV file you gave is not a file.%s\n"
            "\n%sError explanations:%s\n"
            " - Clinica expected the following path to be a file: %s%s%s\n"
            " - If you gave relative path, did you run Clinica on the good folder?" %
            (Fore.RED, Fore.RESET,
             Fore.YELLOW, Fore.RESET,
             Fore.BLUE, tsv_file, Fore.RESET)
        )
    ss_df = pd.io.parsers.read_csv(tsv_file, sep='\t')
    if 'participant_id' not in list(ss_df.columns.values):
        raise ClinicaException(
            "\n%s[Error] The TSV file does not contain participant_id column (path: %s)%s" %
            (Fore.RED, tsv_file, Fore.RESET)
        )
    if 'session_id' not in list(ss_df.columns.values):
        raise ClinicaException(
            "\n%s[Error] The TSV file does not contain session_id column (path: %s)%s" %
            (Fore.RED, tsv_file, Fore.RESET)
        )
    participants = list(ss_df.participant_id)
    sessions = list(ss_df.session_id)

    # Remove potential whitespace in participant_id or session_id
    return [sub.strip(' ') for sub in participants], [ses.strip(' ') for ses in sessions]
