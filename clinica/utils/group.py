"""This module contains utilities to handle groups in Clinica.

See CAPS specifications for details about groups.
"""
from os import PathLike
from typing import List


def check_group_label(group_label: str) -> None:
    """Check that `group_label` is compliant with specifications.

    The group label must be a string composed of alphanumeric
    characters (alphabets or numbers). This means that a group
    label cannot contain spaces nor special characters.

    Parameters
    ----------
    group_label : str
        The group label to verify.

    Raises
    ------
    ValueError :
        If the provided group label is not compliant.
    """
    if not group_label.isalnum():
        raise ValueError(
            "Not valid group_label value: it must be composed only by letters "
            f"and/or numbers (given value: {group_label})."
        )


def _check_group_dir(group_directory: PathLike) -> None:
    from pathlib import Path

    group_directory = Path(group_directory)
    if group_directory.name.startswith("group-"):
        check_group_label(group_directory.name.lstrip("group-"))
    else:
        raise ValueError(f"Group directory {group_directory} is not valid.")


def extract_group_ids(caps_directory: PathLike) -> List[str]:
    """Extract a list of group IDs from a CAPS folder.

    The function searches for sub-folders of `caps_directory/groups`.
    According to the CAPS specifications, these sub-folders should be named
    with their group IDs (e.g. ['group-AD', 'group-HC']).


    Parameters
    ----------
    caps_directory : str
        The CAPS folder to search for group IDs.

    Returns
    -------
    list of str :
        The sorted list of group IDs found inside the provided CAPS folder.

    Raises
    ------
    ValueError :
        If `caps_directory/groups` contains folders which do not
        respect the CAPS naming specifications ("group-{groupID}").
    """
    from pathlib import Path

    caps_directory = Path(caps_directory)

    try:
        groups = [
            group for group in (caps_directory / "groups").iterdir() if group.is_dir()
        ]
    except FileNotFoundError:
        return []
    for group in groups:
        _check_group_dir(group)

    return sorted([group.name for group in groups])
