import re
from pathlib import Path
from typing import List

import nibabel as nib
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from clinica.utils.atlas import (
    AAL2,
    AICHA,
    LPBA40,
    Hammers,
    JHUDTI811mm,
    JHUTracts01mm,
    JHUTracts251mm,
    JHUTracts501mm,
    Neuromorphometrics,
)


def test_atlas_factory_error():
    from clinica.utils.atlas import atlas_factory

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid AtlasName",
    ):
        atlas_factory("foo")


@pytest.mark.parametrize(
    "atlas_name,atlas",
    [
        ("AAL2", AAL2),
        ("AICHA", AICHA),
        ("Hammers", Hammers),
        ("LPBA40", LPBA40),
        ("Neuromorphometrics", Neuromorphometrics),
        ("JHUDTI81", JHUDTI811mm),
        ("JHUTracts0", JHUTracts01mm),
        ("JHUTracts25", JHUTracts251mm),
        ("JHUTracts50", JHUTracts501mm),
    ],
)
def test_atlas_factory(tmp_path, monkeypatch, atlas_name, atlas):
    from clinica.utils.atlas import atlas_factory

    monkeypatch.setenv("FSLDIR", str(tmp_path))
    assert isinstance(atlas_factory(atlas_name), atlas)


@pytest.mark.parametrize(
    (
        "atlas,expected_name,expected_checksum,expected_atlas_filename,"
        "expected_roi_filename,expected_resolution,expected_size"
    ),
    [
        (
            Neuromorphometrics(),
            "Neuromorphometrics",
            "19a50136cd2f8a14357a19ad8a1dc4a2ecb6beb3fc16cb5441f4f2ebaf64a9a5",
            "atlas-Neuromorphometrics_dseg.nii.gz",
            "atlas-Neuromorphometrics_dseg.tsv",
            "1.5x1.5x1.5",
            141,
        ),
        (
            LPBA40(),
            "LPBA40",
            "20826b572bbbdbcdbf28bbd3801dc0c2fed28d1e54bc4fd5027e64ccc6d50374",
            "atlas-LPBA40_dseg.nii.gz",
            "atlas-LPBA40_dseg.tsv",
            "1.5x1.5x1.5",
            57,
        ),
        (
            Hammers(),
            "Hammers",
            "c034a7bce2dcab390a0b72f4e7d04769eb3fe5b990d0e18d89b0ce73339a5376",
            "atlas-Hammers_dseg.nii.gz",
            "atlas-Hammers_dseg.tsv",
            "1.5x1.5x1.5",
            69,
        ),
        (
            AICHA(),
            "AICHA",
            "cab554d5f546720e60f61f536f82c3d355b31fadb5a4d3ce6a050a606d7ef761",
            "atlas-AICHA_dseg.nii.gz",
            "atlas-AICHA_dseg.tsv",
            "1.5x1.5x1.5",
            385,
        ),
        (
            AAL2(),
            "AAL2",
            "f6bc698f778a4b383abd3ce355bfd4505c4aa14708e4a7848f8ee928c2b56b37",
            "atlas-AAL2_dseg.nii.gz",
            "atlas-AAL2_dseg.tsv",
            "1.5x1.5x1.5",
            121,
        ),
    ],
    ids=(
        "Neuromorphometrics",
        "LPBA40",
        "Hammers",
        "AICHA",
        "AAL2",
    ),
)
def test_atlases(
    tmp_path,
    atlas,
    expected_name,
    expected_checksum,
    expected_atlas_filename,
    expected_roi_filename,
    expected_resolution,
    expected_size,
    mocker,
):
    # The following atlases are supposed to be downloaded from a remote server
    # For unit testing, we mock the get_file_from_server function to return
    # a nifti image with fake data realistic enough to make the tests pass
    if expected_name in ("Hammers", "LPBA40", "Neuromorphometrics"):
        mocker.patch(
            "clinica.utils.inputs.get_file_from_server",
            return_value=get_mocked_atlas(tmp_path, atlas),
        )
    assert atlas.name == expected_name
    assert atlas.expected_checksum == expected_checksum
    assert atlas.atlas_filename == expected_atlas_filename
    assert atlas.roi_filename == expected_roi_filename
    assert atlas.tsv_roi.exists()
    assert atlas.spatial_resolution == expected_resolution
    assert atlas.atlas_folder.exists()
    assert_array_equal(atlas.get_index(), np.arange(expected_size))


def get_mocked_atlas(tmp_path, atlas) -> Path:
    """Return the path to the mocked atlas label image.

    We need the mocked image to have:
        - the right resolution (same for all mocked atlases)
        - the right data shape (same for all mocked atlases)
        - the right label values

    To ensure the last point, this function generate a data array
    composed of the cycled expected labels.
    """
    from itertools import cycle

    data_shape = (121, 145, 121)
    gen = cycle(get_labels(atlas))
    mocked_data = []
    while len(mocked_data) < np.prod(data_shape):
        mocked_data.append(next(gen))
    mocked_data = np.array(mocked_data, dtype="float")
    mocked_data = mocked_data.reshape(data_shape)
    mocked_image = nib.Nifti1Image(mocked_data, np.diag([-1.5, 1.5, 1.5, 1.0]))
    nib.save(mocked_image, tmp_path / "mocked_atlas.nii.gz")

    return tmp_path / "mocked_atlas.nii.gz"


def get_labels(atlas) -> List[int]:
    """Return the expected label values for the mocked atlases."""
    if atlas.name == "Neuromorphometrics":
        labels = list(range(143))
        labels.remove(35)
        labels.remove(23)
    elif atlas.name == "LPBA40":
        labels = [
            0,
            21,
            22,
            23,
            24,
            25,
            26,
            27,
            28,
            29,
            30,
            31,
            32,
            33,
            34,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
            49,
            50,
            61,
            62,
            63,
            64,
            65,
            66,
            67,
            68,
            81,
            82,
            83,
            84,
            85,
            86,
            87,
            88,
            89,
            90,
            91,
            92,
            101,
            102,
            121,
            122,
            161,
            162,
            163,
            164,
            165,
            166,
            181,
            182,
        ]
    elif atlas.name == "Hammers":
        labels = list(range(69))
    else:
        raise ValueError(
            f"Atlas {atlas.name} is not supposed to be mocked in this test."
        )
    return labels


@pytest.fixture
def atlas(expected_name, tmp_path, monkeypatch):
    from clinica.utils.atlas import atlas_factory

    monkeypatch.setenv("FSLDIR", str(tmp_path))
    return atlas_factory(expected_name)


@pytest.mark.parametrize(
    "expected_name,expected_checksum,expected_atlas_filename,expected_roi_filename,expected_resolution,expected_size",
    [
        (
            "JHUTracts50",
            "20ff0216d770686838de26393c0bdac38c8104760631a1a2b5f518bc0bbb470a",
            "JHU-ICBM-tracts-maxprob-thr50-1mm.nii.gz",
            "atlas-JHUTract_dseg.tsv",
            "1x1x1",
            18,
        ),
        (
            "JHUTracts25",
            "7cd85fa2be1918fc83173e9bc0746031fd4c08d70d6c81b7b9224b5d3da6d8a6",
            "JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz",
            "atlas-JHUTract_dseg.tsv",
            "1x1x1",
            21,
        ),
        (
            "JHUTracts0",
            "eb1de9413a46b02d2b5c7b77852097c6f42c8a5d55a5dbdef949c2e63b95354e",
            "JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz",
            "atlas-JHUTract_dseg.tsv",
            "1x1x1",
            21,
        ),
        (
            "JHUDTI81",
            "3c3f5d2f1250a3df60982acff35a75b99fd549a05d5f8124a63f78221aa0ec16",
            "JHU-ICBM-labels-1mm.nii.gz",
            "atlas-JHUDTI81_dseg.tsv",
            "1x1x1",
            51,
        ),
    ],
    ids=("JHUTracts50", "JHUTracts25", "JHUTracts0", "JHUDTI81"),
)
def test_atlases_fsl(
    tmp_path,
    atlas,
    expected_name,
    expected_checksum,
    expected_atlas_filename,
    expected_roi_filename,
    expected_resolution,
    expected_size,
    mocker,
):
    mocker.patch("nipype.interfaces.fsl.Info.version", return_value="6.0.5")
    mocked_fsl_dir = tmp_path / "data" / "atlases" / "JHU"
    mocked_fsl_dir.mkdir(parents=True)
    (mocked_fsl_dir / expected_atlas_filename).touch()

    assert atlas.name == expected_name
    assert atlas.expected_checksum == expected_checksum
    assert atlas.atlas_filename == expected_atlas_filename
    assert atlas.roi_filename == expected_roi_filename
    assert atlas.tsv_roi.exists()
    assert atlas.atlas_folder.exists()


@pytest.mark.parametrize(
    "expected_name",
    [
        "AAL2",
        "JHUDTI81",
        "JHUTracts0",
        "JHUTracts25",
        "JHUTracts50",
        "AICHA",
    ],
    ids=["AAL2", "JHUDTI81", "JHUTracts0", "JHUTracts25", "JHUTracts50", "AICHA"],
)
def test_atlas_checksum_error(atlas, expected_name, mocker):
    mocker.patch("nipype.interfaces.fsl.Info.version", return_value="6.0.5")
    mocker.patch("clinica.utils.inputs.compute_sha256_hash", return_value="123")

    with pytest.raises(
        IOError,
        match=re.escape("has an SHA256 checksum (123) differing from expected"),
    ):
        atlas.labels


@pytest.fixture
def test_image() -> nib.Nifti1Image:
    rng = np.random.RandomState(42)
    affine = np.zeros((4, 4), float)
    np.fill_diagonal(affine, [1, 3, 2.33, 1])
    return nib.Nifti1Image(rng.random((2, 2, 2)), affine=affine)


@pytest.mark.parametrize("axis,expected_resolution", [(0, "1"), (1, "3"), (2, "2.33")])
def test_get_resolution_along_axis(test_image, axis, expected_resolution):
    from clinica.utils.atlas import _get_resolution_along_axis

    assert (
        _get_resolution_along_axis(test_image.header, axis=axis) == expected_resolution
    )


def test_get_resolution_along_axis_error(test_image):
    from clinica.utils.atlas import _get_resolution_along_axis

    with pytest.raises(
        ValueError,
        match=(
            "The label image has dimension 3 and axis 3 is therefore not valid. "
            "Please use a value between 0 and 2."
        ),
    ):
        _get_resolution_along_axis(test_image.header, axis=3)
