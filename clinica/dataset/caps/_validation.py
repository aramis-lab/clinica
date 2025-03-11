from pathlib import Path
from typing import Union

from clinica.utils.exceptions import ClinicaCAPSError

from .._dataset_type import DatasetType, get_dataset_type

__all__ = [
    "check_caps_dataset",
]


def check_caps_dataset(directory: Union[str, Path]) -> Path:
    """Check if provided `caps_directory`is a CAPS folder.

    Parameters
    ----------
    directory : str or Path
        The path to the CAPS dataset to check.

    Returns
    -------
    Path :
        The path to the validated CAPS dataset.

    Raises
    ------
    ClinicaCAPSError :
        If the provided path does not exist, or is not a directory.
        If the provided path is a BIDS folder (BIDS and CAPS could be
        swapped by user). We simply check that there is not a folder
        whose name starts with 'sub-' in the provided path (that exists
        in BIDS hierarchy).

    Notes
    -----
    Keep in mind that a CAPS folder can be empty.
    """
    directory = Path(directory)
    if get_dataset_type(directory) != DatasetType.DERIVATIVE:
        raise ClinicaCAPSError(
            f"The directory ({directory}) you provided is not a CAPS directory."
        )
    _check_dataset_description_exists_in_caps(directory)

    sub_folders = [f for f in directory.iterdir() if f.name.startswith("sub-")]
    if len(sub_folders) > 0:
        error_string = (
            "Your CAPS directory contains at least one folder whose name "
            "starts with 'sub-'. Check that you did not swap BIDS and CAPS folders.\n"
            "Folder(s) found that match(es) BIDS architecture:\n"
        )
        for directory in sub_folders:
            error_string += f"\t{directory}\n"
        error_string += (
            "A CAPS directory has a folder 'subjects' at its root, in which "
            "are stored the output of the pipeline for each subject."
        )
        raise ClinicaCAPSError(error_string)
    return directory


def _check_dataset_description_exists_in_caps(directory: Path):
    if not (directory / "dataset_description.json").exists():
        raise ClinicaCAPSError(
            f"The CAPS directory ({directory}) you provided is missing a dataset_description.json file."
        )
