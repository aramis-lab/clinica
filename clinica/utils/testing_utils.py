import json
import random
from functools import partial
from pathlib import Path
from typing import Callable, Dict, Optional

import nibabel as nib
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from clinica.pipelines.dwi.utils import DWIDataset

__all__ = [
    "build_test_image_cubic_object",
    "build_bids_directory",
    "build_caps_directory",
    "build_dwi_dataset",
    "rmtree",
    "assert_nifti_equal",
    "assert_nifti_almost_equal",
    "assert_large_nifti_almost_equal",
]


def build_test_image_cubic_object(
    shape: tuple[int, int, int],
    background_value: float,
    object_value: float,
    object_size: int,
    affine: Optional[np.ndarray] = None,
) -> nib.Nifti1Image:
    """Returns a 3D nifti image of a centered cube.

    Parameters
    ----------
    shape : (int, int, int)
        The shape of the image to generate.

    background_value : float
        The value that should be put in the background.

    object_value : float
        The value that should be put in the cubic object.

    object_size : int
        The size of the cube (in number of voxels).

    affine : np.ndarray, optional
        The affine to use for the generated image.
        If None, np.eye(4) is used.

    Returns
    -------
    nib.Nifti1Image :
        The test image.
    """
    from clinica.utils.image import Bbox3D

    if object_size > min(shape):
        raise ValueError(
            f"Cannot generate an image of dimension {shape} with an object of size {object_size} in it..."
        )
    bbox_coordinates = []
    for c in (dim / 2 for dim in shape):
        bbox_coordinates.append(int(c - object_size / 2))
        bbox_coordinates.append(int(c + object_size / 2))
    bbox = Bbox3D.from_coordinates(*bbox_coordinates)
    data = np.ones(shape, dtype=np.float32) * background_value
    x, y, z = bbox.get_slices()
    data[x, y, z] = object_value
    return nib.Nifti1Image(data, affine=(affine or np.eye(4)))


def build_bids_directory(
    directory: Path,
    subjects_sessions: dict,
    modalities: Optional[dict[str, tuple[str, ...]]] = None,
    write_tsv_files: bool = False,
    random_seed: int = 42,
) -> Path:
    """Build a fake BIDS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : Path
        The path to folder where the fake BIDS dataset should be created.

    subjects_sessions : Dict
        Dictionary containing the subjects and their associated sessions.

    modalities : dict
        Modalities to be created.

    Returns
    -------
    Path :
        The path to the BIDS dataset.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    It only adds (for now...) T1W nifti images for all subjects and sessions.
    """
    from clinica.dataset import BIDS_VERSION

    random.seed(random_seed)
    directory.mkdir(exist_ok=True, parents=True)
    modalities = modalities or {"anat": ("T1w", "flair")}
    extensions = ("nii.gz", "json")
    with open(directory / "dataset_description.json", "w") as fp:
        json.dump({"Name": "Example dataset", "BIDSVersion": str(BIDS_VERSION)}, fp)
    for sub, sessions in subjects_sessions.items():
        (directory / sub).mkdir()
        if write_tsv_files:
            (directory / sub / f"{sub}_sessions.tsv").write_text(
                "\n".join(["session_id"] + sessions)
            )
        for ses in sessions:
            session_path = directory / sub / ses
            session_path.mkdir()
            scans_tsv_text = "filename\tscan_id\n"
            for modality, suffixes in modalities.items():
                (session_path / modality).mkdir()
                for suffix in suffixes:
                    for extension in extensions:
                        scan_path = (
                            session_path
                            / modality
                            / f"{sub}_{ses}_{suffix}.{extension}"
                        )
                        scan_path.touch()
                        if write_tsv_files and "nii" in extension:
                            scans_tsv_text += f"{scan_path.relative_to(session_path)}\t{random.randint(1, 1000000)}\n"
            if write_tsv_files and len(scans_tsv_text) > 0:
                (session_path / f"{sub}_{ses}_scans.tsv").write_text(scans_tsv_text)
    return directory


def build_caps_directory(directory: Path, configuration: dict) -> Path:
    """Build a fake CAPS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : Path
        The path to folder where the fake CAPS dataset should be created.

    configuration : Dict
        Dictionary containing the configuration for building the fake CAPS
        directory. It should have the following structure:
            - "groups": ["group_labels"...]
            - "pipelines": {"pipeline_names": config}, where config is a dictionary
              specifying more details for the files that should be written.
            - "subjects": {"subject_labels": ["session_labels"...]}.

    Returns
    -------
    Path :
        The path to the CAPS dataset.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    """
    from clinica.dataset import BIDS_VERSION, CAPS_VERSION

    directory.mkdir(exist_ok=True, parents=True)
    with open(directory / "dataset_description.json", "w") as fp:
        json.dump(
            {
                "Name": "Example dataset",
                "BIDSVersion": str(BIDS_VERSION),
                "CAPSVersion": str(CAPS_VERSION),
            },
            fp,
        )
    _build_groups(directory, configuration)
    _build_subjects(directory, configuration)
    return directory


def _build_groups(directory: Path, configuration: dict) -> None:
    """Build a fake groups file structure in a CAPS directory if requested."""
    if "groups" in configuration:
        (directory / "groups").mkdir()
        for group_label in configuration["groups"]:
            (directory / "groups" / f"group-{group_label}").mkdir()
            for pipeline_name, pipeline_config in configuration["pipelines"].items():
                (directory / "groups" / f"group-{group_label}" / pipeline_name).mkdir()
                (
                    directory
                    / "groups"
                    / f"group-{group_label}"
                    / pipeline_name
                    / f"group-{group_label}_template.nii.gz"
                ).touch()


def _build_subjects(directory: Path, configuration: dict) -> None:
    """Build a fake subjects file structure in a CAPS directory if requested."""
    if "subjects" in configuration:
        (directory / "subjects").mkdir(exist_ok=True)
        for sub, sessions in configuration["subjects"].items():
            (directory / "subjects" / sub).mkdir(exist_ok=True)
            for ses in sessions:
                (directory / "subjects" / sub / ses).mkdir(exist_ok=True)
                for pipeline_name, pipeline_config in configuration[
                    "pipelines"
                ].items():
                    (directory / "subjects" / sub / ses / pipeline_name).mkdir(
                        exist_ok=True
                    )
                    if pipeline_name == "t1_linear":
                        _build_t1_linear(directory, sub, ses, pipeline_config)
                    if pipeline_name == "pet_linear":
                        _build_pet_linear(directory, sub, ses, pipeline_config)
                    if pipeline_name == "t1":
                        _build_t1(directory, sub, ses, configuration)


def _build_t1_linear(directory: Path, sub: str, ses: str, config: dict) -> None:
    """Build a fake t1-linear file structure in a CAPS directory."""
    uncropped = config.get("uncropped_image", False)
    for filename in (
        f"{sub}_{ses}_space-MNI152NLin2009cSym{'' if uncropped else '_desc-Crop'}_res-1x1x1_T1w.nii.gz",
        f"{sub}_{ses}_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
    ):
        (directory / "subjects" / sub / ses / "t1_linear" / filename).touch()


def _build_pet_linear(directory: Path, sub: str, ses: str, config: dict) -> None:
    """Build a fake pet-linear file structure in a CAPS directory."""
    from clinica.utils.pet import SUVRReferenceRegion, Tracer

    tracer = Tracer(config["acq_label"])
    suvr = SUVRReferenceRegion(config["suvr_reference_region"])
    if config.get("save_PETinT1w", False):
        (
            directory
            / "subjects"
            / sub
            / ses
            / "pet_linear"
            / f"{sub}_{ses}_trc-{tracer.value}_space-T1w_pet.nii.gz"
        ).touch()
    for filename in (
        f"{sub}_{ses}_trc-{tracer.value}_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-{suvr.value}_pet.nii.gz",
        f"{sub}_{ses}_trc-{tracer.value}_space-T1w_rigid.mat",
    ):
        (directory / "subjects" / sub / ses / "pet_linear" / filename).touch()


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


def build_dwi_dataset(
    folder: Path, nb_dwi_volumes: int, nb_b_values: int, nb_b_vectors: int
) -> DWIDataset:
    dwi_data = 4.0 * np.ones((5, 5, 5, nb_dwi_volumes))
    dwi_img = nib.Nifti1Image(dwi_data, affine=np.eye(4))
    nib.save(dwi_img, folder / "foo.nii.gz")
    np.savetxt(folder / "foo.bval", [1000] * nb_b_values)
    b_vectors_data = np.random.random((3, nb_b_vectors))
    np.savetxt(folder / "foo.bvec", b_vectors_data)

    return DWIDataset(
        dwi=folder / "foo.nii.gz",
        b_values=folder / "foo.bval",
        b_vectors=folder / "foo.bvec",
    )
