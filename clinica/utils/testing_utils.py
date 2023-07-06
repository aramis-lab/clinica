import json
import os
from functools import partial
from pathlib import Path
from typing import Callable, Optional

import nibabel as nib
from numpy.testing import assert_array_almost_equal, assert_array_equal


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


def _assert_nifti_relation(
    img1: Path,
    img2: Path,
    assertion_func_data: Callable,
    assertion_func_affine: Optional[Callable] = None,
) -> None:
    """Assert that two nifti images satisfy some relationship.

    Parameters
    ----------
    img1 : Path
        Path to the first image.

    img2 : Path
        Path to the second image.

    assertion_func_data : Callable
        Assertion function for dataobj comparison.
        The function should take two numpy arrays as input.

    assertion_func_affine : Callable, optional
        Assertion function for affine comparison.
        The function should take two numpy arrays as input.
        If not specified, `assertion_func_data` will be used to
        check affine matrices.

    See Also
    --------
    assert_nifti_equal : Check strict equality for two images.
    assert_nifti_almost_equal : Check loose equality for two images.
    """
    assertion_func_affine = assertion_func_affine or assertion_func_data
    img1 = nib.load(img1)
    img2 = nib.load(img2)

    assert img1.shape == img2.shape  # Fail fast
    assertion_func_affine(img1, img2)
    assertion_func_data(img1, img2)


def _assert_affine_equal(img1: nib.Nifti1Image, img2: nib.Nifti1Image) -> None:
    assert_array_equal(img1.affine, img2.affine)


def _assert_dataobj_equal(img1: nib.Nifti1Image, img2: nib.Nifti1Image) -> None:
    assert_array_equal(img1.get_fdata(), img2.get_fdata())


def _assert_affine_almost_equal(
    img1: nib.Nifti1Image, img2: nib.Nifti1Image, decimal: int = 6
) -> None:
    assert_array_almost_equal(img1.affine, img2.affine, decimal=decimal)


def _assert_dataobj_almost_equal(
    img1: nib.Nifti1Image, img2: nib.Nifti1Image, decimal: int = 6
) -> None:
    assert_array_almost_equal(img1.get_fdata(), img2.get_fdata(), decimal=decimal)


def _assert_large_image_dataobj_almost_equal(
    img1: nib.Nifti1Image,
    img2: nib.Nifti1Image,
    decimal: int = 6,
    n_samples: Optional[int] = None,
    verbose: bool = False,
) -> None:
    """This alternative implementation can be used when dealing with large nifti images.

    Memory issue:

    `_assert_dataobj_almost_equal` loads both images data arrays in memory before running
    the equality assertion test which can be problematic when these arrays are large.
    This function will work over the last dimension (usually time, dwi directions...).
    It loops over this dimension and loads in memory only the two data chunks to be compared.

    CPU time issue:

    Even without the previously described memory issue, the comparison can take a lot
    of time for very large images. This function enables to sample volumes for equality
    assertion through the `n_samples` argument.

    TODO: Implement in parallel?

    Parameters
    ----------
    img1 : Nifti1Image
        The first image to be compared.

    img2 : Nifti1Image
        The second image to be compared.

    decimal : int, optional
        The number of decimal for which the loose equality should work.
        Default=6.

    n_samples : int, optional
        If specified, the function will only compare a subset of the data.
        If None, it will compare all data.
        The number of samples specifies the number of "volumes" which will be compared.
        Default=None.

    verbose : bool, optional
        If True, the function prints a message for each compared volume.
        Default=False.

    See Also
    --------
    assert_large_nifti_almost_equal : Check loose equality for two large images.
    """
    import random

    import numpy as np

    volumes = range(0, img1.shape[-1])
    if n_samples:
        volumes = random.sample(volumes, n_samples)
    for volume in volumes:
        if verbose:
            print(f"--> Processing volume {volume}...")
        assert_array_almost_equal(
            np.asarray(img1.dataobj[..., volume]),
            np.asarray(img2.dataobj[..., volume]),
            decimal=decimal,
        )


assert_nifti_equal = partial(
    _assert_nifti_relation,
    assertion_func_data=_assert_dataobj_equal,
    assertion_func_affine=_assert_affine_equal,
)


assert_nifti_almost_equal = partial(
    _assert_nifti_relation,
    assertion_func_data=_assert_dataobj_almost_equal,
    assertion_func_affine=_assert_affine_almost_equal,
)


assert_large_nifti_almost_equal = partial(
    _assert_nifti_relation,
    assertion_func_data=_assert_large_image_dataobj_almost_equal,
    assertion_func_affine=_assert_affine_almost_equal,
)
