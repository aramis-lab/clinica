"""This submodule contains logic specific to handling BIDS datasets."""

from ._dataset_description import BIDS_VERSION, BIDSDatasetDescription, get_bids_version
from ._filename import BIDSFileName, BIDSLabel
from ._queries import (
    get_paths_to_subjects_in_bids_dataset,
    get_sessions_for_subject_in_bids_dataset,
    get_subjects_from_bids_dataset,
)
from ._readme import BIDSReadme
from ._validation import check_bids_dataset

__all__ = [
    "BIDSDatasetDescription",
    "BIDSLabel",
    "BIDSFileName",
    "BIDSReadme",
    "BIDS_VERSION",
    "check_bids_dataset",
    "get_bids_version",
    "get_subjects_from_bids_dataset",
    "get_sessions_for_subject_in_bids_dataset",
    "get_paths_to_subjects_in_bids_dataset",
]
