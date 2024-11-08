"""This module contains utilities to handle groups in Clinica.

See CAPS specifications for details about groups.
"""
from collections import UserString
from pathlib import Path
from typing import Union

__all__ = [
    "GroupLabel",
    "GroupID",
    "extract_group_ids",
]


class GroupLabel(UserString):
    """This is the type defining a label for groups in Clinica.

    The group label must be a string composed of alphanumeric characters (alphabets or numbers).
    This means that a group label cannot contain spaces nor special characters.
    """

    def __init__(self, value: str):
        super().__init__(self.validate(value))

    @classmethod
    def validate(cls, value: str) -> str:
        if value.isalnum():
            return value
        raise ValueError(
            f"Group label '{value}' is not a valid group label: it must be composed only by letters and/or numbers."
        )


class GroupID(UserString):
    """This is the type defining a group ID in Clinica.

    It is basically a group label with a 'group-' prefix.
    """

    def __init__(self, value: str):
        super().__init__(self.validate(value))

    @classmethod
    def validate(cls, value: str) -> str:
        if not value.startswith("group-"):
            raise ValueError(
                f"Group ID '{value}' is not a valid group ID: it must start with 'group-'."
            )
        GroupLabel.validate(value.removeprefix("group-"))
        return value

    @property
    def label(self) -> GroupLabel:
        return GroupLabel(str(self).removeprefix("group-"))

    @classmethod
    def from_label(cls, label: Union[str, GroupLabel]):
        return cls(f"group-{GroupLabel(label)}")


def extract_group_ids(caps_directory: Path) -> list[GroupID]:
    """Extract a list of group IDs from a CAPS folder.

    The function searches for sub-folders of `caps_directory/groups`.
    According to the CAPS specifications, these sub-folders should be named
    with their group IDs (e.g. ['group-AD', 'group-HC']).

    Parameters
    ----------
    caps_directory : Path
        The path to the CAPS folder to search for group IDs.

    Returns
    -------
    list of GroupID :
        The sorted list of group IDs found inside the provided CAPS folder.

    Raises
    ------
    ValueError :
        If `caps_directory/groups` contains folders which do not
        respect the CAPS naming specifications ("group-{groupID}").
    """
    try:
        groups = [
            group for group in (caps_directory / "groups").iterdir() if group.is_dir()
        ]
    except FileNotFoundError:
        return []
    return sorted([GroupID(group.name) for group in groups])
