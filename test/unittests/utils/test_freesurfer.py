import pytest

from clinica.utils.freesurfer import (
    _get_prefix,
    extract_image_id_from_longitudinal_segmentation,
)


@pytest.mark.parametrize(
    "func", [_get_prefix, extract_image_id_from_longitudinal_segmentation]
)
@pytest.mark.parametrize("sub_id", ["", "fooo"])
def test_subject_id_error(sub_id, func):
    with pytest.raises(
        ValueError,
        match=f"The provided Freesurfer ID {sub_id} could not be parsed.",
    ):
        func(sub_id)


@pytest.mark.parametrize(
    "sub_id,expected",
    [
        ("sub-CLNC01_ses-M00", "sub-CLNC01_ses-M00"),
        ("sub-CLNC01_long-M00M18", "sub-CLNC01_long-M00M18"),
        (
            "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18",
            "sub-CLNC01_ses-M00_long-M00M18",
        ),
    ],
)
def test_get_prefix(sub_id, expected):
    from clinica.utils.freesurfer import _get_prefix

    assert _get_prefix(sub_id) == expected


def test_generate_regional_measures():
    from clinica.utils.freesurfer import generate_regional_measures

    with pytest.raises(
        OSError,
        match="Image sub-CLNC01 | ses-M00 does not contain FreeSurfer segmentation",
    ):
        generate_regional_measures("", "sub-CLNC01_ses-M00", "destrieux")


@pytest.fixture
def expected_image_id(sub_id):
    from collections import namedtuple

    image_id = namedtuple("image_id", ["participant_id", "session_id", "long_id"])
    if sub_id == "sub-CLNC01_ses-M00":
        return image_id("sub-CLNC01", "ses-M00", "")
    if sub_id == "sub-CLNC01_long-M00M18":
        return image_id("sub-CLNC01", "", "long-M00M18")
    if sub_id == "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18":
        return image_id("sub-CLNC01", "ses-M00", "long-M00M18")


@pytest.mark.parametrize(
    "sub_id",
    [
        "sub-CLNC01_ses-M00",
        "sub-CLNC01_long-M00M18",
        "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18",
    ],
)
def test_extract_image_id_from_longitudinal_segmentation(sub_id, expected_image_id):
    from clinica.utils.freesurfer import extract_image_id_from_longitudinal_segmentation

    assert extract_image_id_from_longitudinal_segmentation(sub_id) == expected_image_id


def test_get_secondary_stats(tmp_path):
    from clinica.utils.freesurfer import InfoType, get_secondary_stats

    assert get_secondary_stats(tmp_path, InfoType.MEANCURV) == {}
