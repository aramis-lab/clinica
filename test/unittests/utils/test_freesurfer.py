import pytest

from clinica.utils.freesurfer import _get_prefix  # noqa
from clinica.utils.freesurfer import (
    ColumnType,
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


@pytest.mark.parametrize("column_type", ["foo", "parcellation", "segmentation"])
def test_read_stats_file_columns_errors(tmp_path, column_type):
    from clinica.utils.freesurfer import _read_stats_file

    with pytest.raises(
        ValueError,
        match=f"Unknown column type {column_type}. Use either parcellation or segmentation.",
    ):
        _read_stats_file(tmp_path, column_type)


@pytest.mark.parametrize(
    "column_type", [ColumnType.PARCELLATION, ColumnType.SEGMENTATION]
)
def test_read_stats_file_file_errors(tmp_path, column_type):
    from clinica.utils.freesurfer import _read_stats_file

    with pytest.raises(
        FileNotFoundError,
        match=f"Stats file {tmp_path / 'foo.stats'} could not be found.",
    ):
        _read_stats_file(tmp_path / "foo.stats", column_type)


@pytest.fixture
def fake_stats_dataframe(columns, n_rows):
    """Generate a fake dummy stats Pandas dataframe for testing purposes."""
    import pandas as pd

    return pd.DataFrame({k: ["foo"] * n_rows for k in columns})


@pytest.mark.parametrize(
    "columns,column_type,n_rows",
    [
        (
            [
                "StructName",
                "NumVert",
                "SurfArea",
                "GrayVol",
                "ThickAvg",
                "ThickStd",
                "MeanCurv",
                "GausCurv",
                "FoldInd",
                "CurvInd",
            ],
            ColumnType.PARCELLATION,
            6,
        ),
        (
            [
                "Index",
                "SegId",
                "NVoxels",
                "Volume_mm3",
                "StructName",
                "normMean",
                "normStdDev",
                "normMin",
                "normMax",
                "normRange",
            ],
            ColumnType.SEGMENTATION,
            6,
        ),
    ],
)
def test_read_stats_file(tmp_path, columns, column_type, n_rows, fake_stats_dataframe):
    from clinica.utils.freesurfer import _read_stats_file

    fake_stats_dataframe.to_csv(
        tmp_path / "stats_file.csv", sep=" ", index=False, header=False
    )
    df = _read_stats_file(tmp_path / "stats_file.csv", column_type)
    assert len(df) == n_rows
    assert set(df.columns) == set(columns)


@pytest.mark.parametrize(
    "atlas,expected",
    [
        ("desikan", "aparc.stats"),
        ("destrieux", "aparc.a2009s.stats"),
        ("ba", "BA_exvivo.stats"),
        ("foo", "foo.stats"),
    ],
)
def test_get_stats_filename_for_atlas(tmp_path, atlas, expected):
    from clinica.utils.freesurfer import _get_stats_filename_for_atlas

    df = _get_stats_filename_for_atlas(tmp_path, atlas)
    for hemi in ("lh", "rh"):
        assert df[hemi] == tmp_path / f"{hemi}.{expected}"


def test_generate_tsv_for_parcellation_errors(tmp_path):
    from clinica.utils.freesurfer import _generate_tsv_for_parcellation

    not_supported_atlas = "foo"
    with pytest.raises(
        FileNotFoundError,
        match=f"Stats file {tmp_path / 'lh.foo.stats'} could not be found.",
    ):
        _generate_tsv_for_parcellation(
            tmp_path,
            tmp_path,
            "sub-CLNC01_ses-M00",
            [not_supported_atlas],
        )


"""@pytest.mark.parametrize(
    "prefix",
    [
        ("sub-CLNC01_ses-M00"),
        ("sub-CLNC01_long-M00M18"),
        ("sub-CLNC01_ses-M00_long-M00M18"),
    ],
)
def test_generate_tsv_for_parcellation(tmp_path, prefix):
    from clinica.utils.freesurfer import _generate_tsv_for_parcellation

    stats_folder = tmp_path / "stats"
    stats_folder.mkdir()
    _generate_tsv_for_parcellation(stats_folder, tmp_path, prefix, atlases)"""


def test_generate_regional_measures_error():
    from clinica.utils.freesurfer import generate_regional_measures

    with pytest.raises(
        FileNotFoundError,
        match="Image sub-CLNC01 | ses-M00 does not contain FreeSurfer segmentation",
    ):
        generate_regional_measures("", "sub-CLNC01_ses-M00", ["destrieux"])

    with pytest.raises(
        ValueError,
        match="atlases should be a list of strings. <class 'str'> was provided instead.",
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
