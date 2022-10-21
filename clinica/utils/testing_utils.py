import json
import os
from pathlib import Path


def build_bids_directory(directory: os.PathLike, subjects_sessions: dict) -> None:
    """Build a fake BIDS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : PathLike
        Path to folder where the fake BIDS dataset should be created.

    subjects_sessions : Dict
        Dictionary containing the subjects and their associated sessions.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    It only adds (for now...) T1W nifti images for all subjects and sessions.
    """
    directory = Path(directory)
    data_types = {"anat"}
    suffixes = {"T1w", "flair"}
    extensions = {"nii.gz"}
    with open(directory / "dataset_description.json", "w") as fp:
        json.dump({"Name": "Example dataset", "BIDSVersion": "1.0.2"}, fp)
    for sub, sessions in subjects_sessions.items():
        (directory / sub).mkdir()
        for ses in sessions:
            (directory / sub / ses).mkdir()
            for data_type in data_types:
                (directory / sub / ses / data_type).mkdir()
                for suffix in suffixes:
                    for extension in extensions:
                        (
                            directory
                            / sub
                            / ses
                            / data_type
                            / f"{sub}_{ses}_{suffix}.{extension}"
                        ).touch()


def build_caps_directory(directory: os.PathLike, configuration: dict) -> None:
    """Build a fake CAPS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : PathLike
        Path to folder where the fake CAPS dataset should be created.

    configuration : Dict
        Dictionary containing the configuration for building the fake CAPS
        directory. It should have the following structure:
            - "groups": ["group_labels"...]
            - "pipelines": ["pipeline_names"...]
            - "subjects": {"subject_labels": ["session_labels"...]}.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    """
    directory = Path(directory)
    _build_groups(directory, configuration)
    _build_subjects(directory, configuration)


def _build_groups(directory: Path, configuration: dict) -> None:
    """Build a fake groups file structure in a CAPS directory if requested."""
    if "groups" in configuration:
        (directory / "groups").mkdir()
        for group_label in configuration["groups"]:
            (directory / "groups" / f"group-{group_label}").mkdir()
            for pipeline in configuration["pipelines"]:
                (directory / "groups" / f"group-{group_label}" / pipeline).mkdir()
                (
                    directory
                    / "groups"
                    / f"group-{group_label}"
                    / pipeline
                    / f"group-{group_label}_template.nii.gz"
                ).touch()


def _build_subjects(directory: Path, configuration: dict) -> None:
    """Build a fake subjects file structure in a CAPS directory if requested."""
    if "subjects" in configuration:
        (directory / "subjects").mkdir()
        for sub, sessions in configuration["subjects"].items():
            (directory / "subjects" / sub).mkdir()
            for ses in sessions:
                (directory / "subjects" / sub / ses).mkdir()
                for pipeline in configuration["pipelines"]:
                    (directory / "subjects" / sub / ses / pipeline).mkdir()
                    if pipeline == "t1_linear":
                        _build_t1_linear(directory, sub, ses)
                    if pipeline == "t1":
                        _build_t1(directory, sub, ses, configuration)


def _build_t1_linear(directory: Path, sub: str, ses: str) -> None:
    """Build a fake t1-linear file structure in a CAPS directory."""
    (
        directory
        / "subjects"
        / sub
        / ses
        / "t1_linear"
        / f"{sub}_{ses}_T1w_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz"
    ).touch()


def _build_t1(directory: Path, sub: str, ses: str, configuration: dict) -> None:
    """Build a fake t1 file structure in a CAPS directory."""
    (directory / "subjects" / sub / ses / "t1" / "spm").mkdir()
    for folder in ["dartel", "segmentation"]:
        (directory / "subjects" / sub / ses / "t1" / "spm" / folder).mkdir()
    for group_label in configuration["groups"]:
        (
            directory
            / "subjects"
            / sub
            / ses
            / "t1"
            / "spm"
            / "dartel"
            / f"group-{group_label}"
        ).mkdir()
    (
        directory
        / "subjects"
        / sub
        / ses
        / "t1"
        / "spm"
        / "dartel"
        / f"group-{group_label}"
        / f"{sub}_{ses}_T1w_target-{group_label}_transformation-forward_deformation.nii.gz"
    ).touch()
    for tissue in ["csf", "graymatter", "whitematter"]:
        (
            directory
            / "subjects"
            / sub
            / ses
            / "t1"
            / "spm"
            / "dartel"
            / f"group-{group_label}"
            / f"{sub}_{ses}_T1w_segm-{tissue}_space-Ixi549Space_modulated-on_probability.nii.gz"
        ).touch()
    for folder in ["dartel_input", "native_space", "normalized_space"]:
        (
            directory / "subjects" / sub / ses / "t1" / "spm" / "segmentation" / folder
        ).mkdir()
    for tissue in ["csf", "graymatter", "whitematter"]:
        (
            directory
            / "subjects"
            / sub
            / ses
            / "t1"
            / "spm"
            / "segmentation"
            / "dartel_input"
            / f"{sub}_{ses}_T1w_segm-{tissue}_dartelinput.nii.gz"
        ).touch()
        (
            directory
            / "subjects"
            / sub
            / ses
            / "t1"
            / "spm"
            / "segmentation"
            / "native_space"
            / f"{sub}_{ses}_T1w_segm-{tissue}_probability.nii.gz"
        ).touch()
        (
            directory
            / "subjects"
            / sub
            / ses
            / "t1"
            / "spm"
            / "segmentation"
            / "normalized_space"
            / f"{sub}_{ses}_T1w_segm-{tissue}_space-Ixi549Space_modulated-off_probability.nii.gz"
        ).touch()


def rmtree(f: Path):
    if f.is_file():
        f.unlink()
    else:
        for child in f.iterdir():
            rmtree(child)
        f.rmdir()
