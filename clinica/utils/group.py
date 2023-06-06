"""This module contains utilities to handle groups in Clinica.

See CAPS specifications for details about groups.
"""


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


def extract_group_ids(caps_directory: str) -> list:
    """Extract a list of group IDs from a CAPS folder.

    The function searches for sub-folders of `caps_directory/groups`.
    According to the CAPS specifications, these sub-folders should be named
    with their group IDs (e.g. ['group-AD', 'group-HC']).

    If the provided CAPS folder does not contain a "groups" subfolder,
    the function returns a list with an empty string.

    Parameters
    ----------
    caps_directory : str
        The CAPS folder to search for group IDs.

    Returns
    -------
    list of str :
        The list of group IDs found inside the provided CAPS folder.
    """
    import os

    try:
        return os.listdir(os.path.join(caps_directory, "groups"))
    except FileNotFoundError:
        return [""]
