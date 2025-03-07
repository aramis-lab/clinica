import json
from enum import Enum
from pathlib import Path
from typing import Union

__all__ = [
    "DatasetType",
    "get_dataset_type",
    "check_dataset",
]


class DatasetType(str, Enum):
    """Defines the possible types of datasets in Clinica."""

    RAW = "raw"
    DERIVATIVE = "derivative"


def get_dataset_type(directory: Union[str, Path]) -> DatasetType:
    """Determine the type of the dataset stored in `input_dir`.

    Parameters
    ----------
    directory : str or Path
        The input folder.

    Returns
    -------
    DatasetType :
        The type of the dataset.

    Raises
    ------
    ClinicaDatasetError :
        If the provided directory is not an existing folder.
    """
    from clinica.utils.exceptions import ClinicaDatasetError

    directory = Path(directory)
    if not directory.is_dir():
        raise ClinicaDatasetError(
            f"The directory you gave is not a folder.\n"
            "Error explanations:\n"
            f"\t- Clinica expected the following path to be a folder: {directory}\n"
            "\t- If you gave relative path, did you run Clinica on the good folder?"
        )
    _check_dataset_description_exists_in_dataset(directory)
    try:
        with open(directory / "dataset_description.json", "r") as fp:
            metadata = json.load(fp)
        dataset_type = DatasetType(metadata["DatasetType"])
    except (json.decoder.JSONDecodeError, KeyError):
        raise ClinicaDatasetError(
            f"The directory ({directory}) you provided has a badly formatted dataset_description.json file."
        )
    return dataset_type


def _check_dataset_description_exists_in_dataset(directory: Path):
    from clinica.utils.exceptions import ClinicaDatasetError

    if not (directory / "dataset_description.json").exists():
        raise ClinicaDatasetError(
            f"The directory ({directory}) you provided is missing a dataset_description.json file."
        )


def _list_subjects_sub_folders(root_dir: Path, groups_dir: Path) -> list[Path]:
    from clinica.utils.stream import cprint

    warning_msg = (
        f"Could not determine if {groups_dir.parent} is a CAPS or BIDS directory. "
        "Clinica will assume this is a CAPS directory."
    )
    folder_content = [f for f in root_dir.iterdir()]
    subjects_sub_folders = [
        sub for sub in folder_content if (sub.name.startswith("sub-") and sub.is_dir())
    ]
    if len(subjects_sub_folders) == 0 and not groups_dir.is_dir():
        cprint(msg=warning_msg, lvl="warning")
    return subjects_sub_folders


def check_dataset(directory: Union[str, Path]) -> None:
    """Check that the provided directory hosts a valid BIDS or CAPS dataset.

    Parameters
    ----------
    directory : str or Path
        The path to the dataset to be checked.
    """
    dataset_type = get_dataset_type(directory)
    if dataset_type == DatasetType.RAW:
        from .bids import check_bids_dataset

        check_bids_dataset(directory)
    if dataset_type == DatasetType.DERIVATIVE:
        from .caps import check_caps_dataset

        check_caps_dataset(directory)
