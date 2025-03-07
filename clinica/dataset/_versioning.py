from enum import Enum
from typing import Optional, Union

from packaging.version import Version

__all__ = [
    "VersionComparisonPolicy",
    "are_versions_compatible",
]


class VersionComparisonPolicy(str, Enum):
    """Defines the different ways we can compare version numbers in Clinica.

    STRICT: version numbers have to match exactly.
    MINOR : version numbers have to have the same major and minor numbers.
    MAJOR: version numbers only need to share the same major number.
    """

    STRICT = "strict"
    MINOR = "minor"
    MAJOR = "major"


def are_versions_compatible(
    version1: Union[str, Version],
    version2: Union[str, Version],
    policy: Optional[Union[str, VersionComparisonPolicy]] = None,
) -> bool:
    """Returns whether the two provided versions are compatible or not depending on the policy.

    Parameters
    ----------
    version1 : str or Version
        The first version number to compare.

    version2 : str or Version
        The second version number to compare.

    policy : str or VersionComparisonPolicy, optional
        The policy under which to compare version1 with version2.
        By default, a strict policy is used, meaning that version
        numbers have to match exactly.

    Returns
    -------
    bool :
        True if version1 is 'compatible' with version2, False otherwise.
    """
    if isinstance(version1, str):
        version1 = Version(version1)
    if isinstance(version2, str):
        version2 = Version(version2)
    if policy is None:
        policy = VersionComparisonPolicy.STRICT
    else:
        policy = VersionComparisonPolicy(policy)
    if policy == VersionComparisonPolicy.STRICT:
        return version1 == version2
    if policy == VersionComparisonPolicy.MINOR:
        return version1.major == version2.major and version1.minor == version2.minor
    if policy == VersionComparisonPolicy.MAJOR:
        return version1.major == version2.major
