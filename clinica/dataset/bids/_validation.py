from pathlib import Path
from typing import Union

from clinica.utils.exceptions import ClinicaBIDSError

from .._dataset_type import DatasetType, get_dataset_type

__all__ = ["check_bids_dataset"]


def check_bids_dataset(directory: Union[str, Path]) -> Path:
    """Check if provided `bids_directory` is a BIDS folder.

    Parameters
    ----------
    directory : str or Path
        The path to the BIDS dataset to check.

    Returns
    -------
    Path :
        The path to the validated BIDS dataset.

    Raises
    ------
    ClinicaBIDSError :
        If the provided path does not exist, or is not a directory.
        If the provided path is a CAPS folder (BIDS and CAPS could
        be swapped by user). We simply check that there is not a folder
        called 'subjects' in the provided path (that exists in CAPS hierarchy).
        If the provided folder is empty.
        If the provided folder does not contain at least one directory whose
        name starts with 'sub-'.
    """
    directory = Path(directory)
    if get_dataset_type(directory) != DatasetType.RAW:
        raise ClinicaBIDSError(
            f"The directory ({directory}) you provided is not a valid BIDS directory."
        )
    _check_bids_is_not_empty(directory)
    _check_bids_has_at_least_one_subject_folder(directory)
    return directory


def _check_bids_is_not_empty(directory: Path):
    if (
        len([f for f in directory.iterdir() if f.name != "dataset_description.json"])
        == 0
    ):
        raise ClinicaBIDSError(
            f"The BIDS directory you provided is empty. ({directory})."
        )


def _check_bids_has_at_least_one_subject_folder(directory: Path):
    if len([f for f in directory.iterdir() if f.name.startswith("sub-")]) == 0:
        raise ClinicaBIDSError(
            "Your BIDS directory does not contains a single folder whose name "
            "starts with 'sub-'. Check that your folder follow BIDS standard."
        )
