"""This module contains utilities for longitudinal pipelines. See CAPS specifications for details about long ID."""
from pathlib import Path
from typing import List, Optional


def get_long_id(session_ids: List[str]) -> str:
    """Extract longitudinal ID from a list of session IDs.

    This will create a unique identifier for a participant and
    its corresponding sessions. Sessions labels are sorted alphabetically
    before being merged in order to generate the longitudinal ID.

    Parameters
    ----------
    session_ids : list of str
        List of session IDs (e.g. ["ses-M000"] or ["ses-M000", "ses-M018", "ses-M036"]).

    Returns
    -------
    str :
        Longitudinal ID.

    Examples
    --------
    >>> from clinica.utils.longitudinal import get_long_id
    >>> get_long_id(['ses-M000'])
    'long-M000'
    >>> get_long_id(['ses-M000', 'ses-M018', 'ses-M036'])
    'long-M000M018M036'
    >>> get_long_id(['ses-M018', 'ses-M036', 'ses-M000'])
    'long-M000M018M036'
    """
    if not all([session_id.startswith("ses-") for session_id in session_ids]):
        raise ValueError(
            "Expected a list of session IDs of the form ses-XXX, "
            f"but received {session_ids} instead."
        )
    return "long-" + "".join(
        [session_id.lstrip("ses-") for session_id in sorted(session_ids)]
    )


def get_participants_long_id(
    participant_ids: List[str], session_ids: List[str]
) -> List[str]:
    """Extract list of longitudinal IDs from a set of participant and session IDs.

    Parameters
    ----------
    participant_ids : list of str
        List of participant IDs for which to compute the longitudinal IDs.

    session_ids : list of str
        List of session IDs for which to compute the longitudinal IDs.

    Returns
    -------
    list of str :
        The computed longitudinal IDs.

    Examples
    --------
    >>> from clinica.utils.longitudinal import get_participants_long_id
    >>> get_participants_long_id(['sub-CLNC01', 'sub-CLNC01', 'sub-CLNC02'], ['ses-M000', 'ses-M018', 'ses-M000'])
    ['long-M000M018', 'long-M000M018', 'long-M000']
    """
    from .participant import get_unique_subjects

    _, sessions_for_each_subject = get_unique_subjects(participant_ids, session_ids)

    long_ids = []
    for sessions in sessions_for_each_subject:
        long_ids += [get_long_id(sessions)] * len(sessions)

    return long_ids


def save_long_id(
    session_ids: List[str],
    output_dir: Path,
    file_name: Optional[str] = None,
) -> None:
    """Save the list of session IDs to given `file_name`.

    Parameters
    ----------
    session_ids : list of str
        The list of session IDs to save.

    output_dir : Path
        The path to the output directory in which to save the session IDs.

    file_name : str, optional
        The file name to use. If None, this will be computed as:
        '{longitudinal_id}_sessions.tsv'
    """
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    file_name = file_name or f"{get_long_id(session_ids)}_sessions.tsv"
    content = "\n".join(sorted(session_ids))
    with open(output_dir / file_name, "w") as fp:
        fp.write(f"session_id\n{content}\n")


def read_sessions(
    caps_dir: Path,
    participant_id: str,
    longitudinal_id: str,
) -> List[str]:
    """Extract sessions IDs from `caps_dir`/subjects/`participant_id`/`long_id`/`long_id`_sessions.tsv.

    Parameters
    ----------
    caps_dir : Path
        Path to CAPS folder.

    participant_id : str
        ID of subject for which to extract session IDs.

    longitudinal_id : str
        Longitudinal ID for which to extract session IDs.

    Returns
    -------
    List of str :
        The extracted list of session IDs.

    Raises
    ------
    ClinicaException :
        If expected session TSV file does not exist.
        If 'session_id' is not in the session dataframe.
    """
    import pandas as pd

    from clinica.utils.exceptions import ClinicaException

    sessions_file = (
        caps_dir
        / "subjects"
        / participant_id
        / longitudinal_id
        / f"{longitudinal_id}_sessions.tsv"
    )

    if not sessions_file.is_file():
        raise ClinicaException(
            "The TSV file with sessions associated "
            f"to {participant_id} for longitudinal ID {longitudinal_id} is missing "
            f"(expected path: {sessions_file})."
        )
    df = pd.read_csv(sessions_file, sep="\t")
    if "session_id" not in df.columns:
        raise ClinicaException(
            f"The TSV file does not contain session_id column (path: {sessions_file})."
        )

    return list(df.session_id)
