from pathlib import Path
from typing import Union

__all__ = [
    "get_subjects_from_bids_dataset",
    "get_sessions_for_subject_in_bids_dataset",
    "get_paths_to_subjects_in_bids_dataset",
]


def get_subjects_from_bids_dataset(directory: Union[str, Path]) -> list[str]:
    """Given a BIDS compliant dataset, return the list of all the subjects available.

    Parameters
    ----------
    directory : str or Path
        The path to the BIDS directory.

    Returns
    -------
    list[str] :
        List of subject IDs available in this BIDS dataset.

    See also
    --------
    get_bids_sess_list
    get_bids_subjs_paths
    """
    from ._validation import check_bids_dataset

    return _filter_folder_names(check_bids_dataset(directory), "sub-*")


def _filter_folder_names(folder: Path, pattern: str) -> list[str]:
    return [d.name for d in folder.glob(pattern) if d.is_dir()]


def get_sessions_for_subject_in_bids_dataset(
    directory: Union[str, Path], subject: str
) -> list[str]:
    """Given a BIDS compliant dataset and a specific subject, return the list of sessions available.

    Parameters
    ----------
    directory : str or Path
        The path to the subject folder for which to list the sessions.

    subject : str
        The subject for which to list the sessions.

    Returns
    -------
    List[str] :
        The list of session names for this subject.

    See also
    --------
    get_subjects_from_bids_dataset
    get_paths_to_subjects_in_bids_dataset
    """
    from ._validation import check_bids_dataset

    return _filter_folder_names(check_bids_dataset(directory) / subject, "ses-*")


def get_paths_to_subjects_in_bids_dataset(directory: Union[str, Path]) -> list[Path]:
    """Given a BIDS compliant dataset, returns the list of all paths to the subjects folders.

    Parameters
    ----------
    directory : str or Path
        The path to the BIDS directory.

    Returns
    -------
    List[Path] :
        List of paths to the subjects folders.

    See also
    --------
    get_subjects_from_bids_dataset
    get_sessions_from_bids_dataset
    """
    from ._validation import check_bids_dataset

    return [d for d in check_bids_dataset(directory).glob("sub-*") if d.is_dir()]
