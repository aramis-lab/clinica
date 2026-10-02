"""Contains duplicate utils that are modified due to clinicaio usage. Should replace the utils in the future when all converters are ported."""

import json
from typing import IO, Optional

from attrs import define
from clinicaio import BIDSDataset, BIDSDatasetDescription, BIDSDatasetType

from clinica.converters.study_models import StudyName

BIDS_VALIDATOR_CONFIG = {
    "ignore": [
        # Possibly dcm2nii(x) errors
        "NIFTI_UNIT",
        "INCONSISTENT_PARAMETERS",
        # fMRI-specific errors
        "SLICE_TIMING_NOT_DEFINED",
        "NIFTI_PIXDIM4",
        "BOLD_NOT_4D",
        "REPETITION_TIME_MUST_DEFINE",
        "TASK_NAME_MUST_DEFINE",
        # Won't fix errors
        "MISSING_SESSION",  # Allows subjects to have different sessions
        "INCONSISTENT_SUBJECTS",  # Allows subjects to have different modalities
        "SCANS_FILENAME_NOT_MATCH_DATASET",  # Necessary until PET is added to BIDS standard
        "CUSTOM_COLUMN_WITHOUT_DESCRIPTION",  # We won't create these JSON files as clinical description
        # is already done in TSV files of clinica.
        "NO_AUTHORS",  # Optional field in dataset_description.json
    ],
    "warn": [],
    "error": [],
    "ignoredFiles": [],
}


@define
class BIDSReadme:
    """Model representing the content for a BIDS Readme."""

    name: str
    link: str
    description: str

    def write(self, to: IO[str]):
        to.write(self.to_str())

    def to_str(self) -> str:
        from importlib.metadata import version

        return (
            f"This BIDS directory was generated with Clinica v{version('clinica')}.\n"
            f"More information on https://www.clinica.run\n"
            f"\n"
            f"Study: {self.name}\n"
            f"\n"
            f"{self.description}\n\n"
            f"Find more about it and about the data user agreement: {self.link}"
        )


def _write_bidsignore(bids_dataset: BIDSDataset) -> None:
    """Write `.bidsignore` file at the root of the BIDS directory."""
    with bids_dataset.write_root_file(".bidsignore", write_binary=False) as f:
        # FIXME: outdated comment??
        # pet/ is necessary until PET is added to BIDS standard
        f.write("\n".join(["swi/\n"]))
        f.write("\n".join(["conversion_info/"]))


def _write_bids_validator_config(bids_dataset: BIDSDataset) -> None:
    """Write `.bids-validator-config.json` at the root of the BIDS directory."""
    with bids_dataset.write_root_file(
        ".bids-validator-config.json", write_binary=False
    ) as f:
        json.dump(BIDS_VALIDATOR_CONFIG, f, skipkeys=True, indent=4)


def write_modality_agnostic_files(
    bids_dataset: BIDSDataset,
) -> None:
    # todo : should replace the original function when all converters migrate
    """
    Write the files .bidsignore and .bids-validator-config.json
    at the root of the BIDS directory.
    Parameters
    ----------
    study_name : StudyName
        The name of the study (Ex ADNI).
    bids_dataset : clinicaio.BIDSDataset
        The BIDS dataset.
    """
    _write_bids_validator_config(bids_dataset)
    _write_bidsignore(bids_dataset)


def _get_dataset_description(
    study_name: StudyName,
    bids_version: Optional[str] = None,
) -> BIDSDatasetDescription:
    from clinica.dataset.bids._dataset_description import BIDS_VERSION

    return BIDSDatasetDescription.new(
        BIDSDatasetType.RAW,
        name=study_name.value,
        bids_version=bids_version or str(BIDS_VERSION),
    )
