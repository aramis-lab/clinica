import numpy as np
import pandas as pd
import pytest

from clinica.converters.adni_to_bids._utils import ADNIModalityConverter


def assert_frame_equal(
    df1, df2, check_index: bool = True, check_column_type: bool = True
):
    from pandas.testing import assert_frame_equal as _assert_frame_equal

    if check_index:
        _assert_frame_equal(df1, df2, check_dtype=check_column_type)
    else:
        _assert_frame_equal(
            df1.reset_index(drop=True),
            df2.reset_index(drop=True),
            check_dtype=check_column_type,
        )


@pytest.mark.parametrize(
    "step_value,expected",
    [(2, ADNIModalityConverter.PET_FDG), (4, ADNIModalityConverter.PET_FDG_UNIFORM)],
)
def test_get_modality_from_adni_preprocessing_step(step_value, expected):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        ADNIPreprocessingStep,
        _get_modality_from_adni_preprocessing_step,
    )

    assert (
        _get_modality_from_adni_preprocessing_step(
            ADNIPreprocessingStep.from_step_value(step_value)
        )
        == expected
    )


@pytest.mark.parametrize("step_value", [0, 1, 3, 5])
def test_get_modality_from_adni_preprocessing_step_error(step_value):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        ADNIPreprocessingStep,
        _get_modality_from_adni_preprocessing_step,
    )

    with pytest.raises(
        ValueError,
        match="The ADNI preprocessing step",
    ):
        _get_modality_from_adni_preprocessing_step(
            ADNIPreprocessingStep.from_step_value(step_value)
        )


@pytest.fixture
def input_df():
    return pd.DataFrame(
        {
            "foo": ["foo1", "foo2", "foo3"],
            "bar": [1, 2, 3],
            "baz": [True, False, False],
            "foobar": [4, 5, 6],
        }
    )


@pytest.mark.parametrize(
    "required_columns",
    [set(), {"bar"}, {"foo", "baz"}, {"foo", "baz", "bar", "foobar"}],
)
def test_load_df_with_column_check(
    tmp_path, input_df: pd.DataFrame, required_columns: set[str]
):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _load_df_with_column_check,
    )

    input_df.to_csv(tmp_path / "data.csv", index=False)

    assert_frame_equal(
        _load_df_with_column_check(tmp_path, "data", required_columns), input_df
    )


def test_load_df_with_column_check_errors(tmp_path, input_df: pd.DataFrame):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _load_df_with_column_check,
    )

    input_df.to_csv(tmp_path / "data.csv", index=False)

    with pytest.raises(
        ValueError,
        match="Missing",
    ):
        _load_df_with_column_check(tmp_path, "data", {"foo", "foobaz"})


EXPECTED_FDG_DF_COLUMNS = [
    "Phase",
    "Subject_ID",
    "VISCODE",
    "Visit",
    "Sequence",
    "Scan_Date",
    "Study_ID",
    "Series_ID",
    "Image_ID",
    "Original",
]


@pytest.fixture
def expected_fdg_df_columns() -> list[str]:
    return EXPECTED_FDG_DF_COLUMNS


@pytest.fixture
def expected_images_df_columns() -> list[str]:
    return EXPECTED_FDG_DF_COLUMNS + ["Is_Dicom", "Path"]


def test_get_pet_fdg_columns(expected_fdg_df_columns: list[str]):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _get_pet_fdg_columns,
    )

    assert _get_pet_fdg_columns() == expected_fdg_df_columns


def test_compute_fdg_pet_paths_empty(tmp_path, expected_images_df_columns: list[str]):
    """Checks that _compute_fdg_pet_paths returns an empty dataframe with the right
    columns when the provided list of subjects is empty.
    """
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        ADNIPreprocessingStep,
        _compute_fdg_pet_paths,
    )

    images = _compute_fdg_pet_paths(
        tmp_path, tmp_path, [], tmp_path, ADNIPreprocessingStep.from_step_value(2)
    )

    assert len(images) == 0
    assert images.columns.tolist() == expected_images_df_columns


def test_compute_fdg_pet_paths_column_errors(tmp_path):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        ADNIPreprocessingStep,
        _compute_fdg_pet_paths,
    )

    csv_dir = tmp_path / "csv"
    csv_dir.mkdir()
    for f in ("PETQC", "PETC3", "PET_META_LIST"):
        pd.DataFrame(columns=["Subject"]).to_csv(csv_dir / f"{f}.csv")
    with pytest.raises(
        ValueError,
        match="Missing column",
    ):
        _compute_fdg_pet_paths(
            tmp_path,
            csv_dir,
            ["sub-01"],
            tmp_path,
            ADNIPreprocessingStep.from_step_value(2),
        )


@pytest.mark.parametrize(
    "row,expected",
    [
        (pd.Series({"Subject_ID": "031_S_0294", "VISCODE": "bl"}), True),
        (pd.Series({"Subject_ID": "031_S_0294", "VISCODE": "m36"}), False),
        (pd.Series({"Subject_ID": "037_S_1421", "VISCODE": "m36", "foo": "bar"}), True),
        (pd.Series({"Subject_ID": "037_S_1078", "VISCODE": "m36"}), True),
        (pd.Series({"Subject_ID": "941_S_1195", "VISCODE": "m48"}), True),
        (pd.Series({"Subject_ID": "005_S_0223", "VISCODE": "m12"}), True),
        (
            pd.Series({"Subject_ID": "123_S_456", "VISCODE": "m12", "foooo": "bar"}),
            False,
        ),
        (pd.Series({"Subject_ID": "foo", "VISCODE": "bar"}), False),
    ],
)
def test_is_visit_a_conversion_error(row: pd.Series, expected: bool):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _is_visit_a_conversion_error,
    )

    assert _is_visit_a_conversion_error(row) is expected


def test_remove_known_conversion_errors():
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _remove_known_conversion_errors,
    )

    input_df = pd.DataFrame(
        {
            "Subject_ID": ["foo", "031_S_0294", "031_S_0294", "123_S_456", "123_S_456"],
            "VISCODE": ["bar", "bl", "m12", "bl", "m48"],
            "foo": ["bar", "baz", "foobar", "foobaz", "foobarbaz"],
        }
    )
    expected_df = pd.DataFrame(
        {
            "Subject_ID": ["foo", "031_S_0294", "123_S_456", "123_S_456"],
            "VISCODE": ["bar", "m12", "bl", "m48"],
            "foo": ["bar", "foobar", "foobaz", "foobarbaz"],
        }
    )
    assert_frame_equal(
        _remove_known_conversion_errors(input_df), expected_df, check_index=False
    )


def test_convert_subject_to_rid():
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _convert_subject_to_rid,
    )

    assert _convert_subject_to_rid("123_S_4567") == 4567


@pytest.mark.parametrize(
    "subject", ["", "123_S_456a", "123_P_4567", "000S6699", "123", "4567"]
)
def test_convert_subject_to_rid_error(subject: str):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _convert_subject_to_rid,
    )

    with pytest.raises(
        ValueError,
        match="Cannot convert the subject",
    ):
        _convert_subject_to_rid(subject)


def test_build_pet_qc_all_studies_for_subject_empty_dataframes():
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _build_pet_qc_all_studies_for_subject,
    )

    df1 = pd.DataFrame(columns=["PASS", "RID"])
    df2 = pd.DataFrame(columns=["SCANQLTY", "RID", "SCANDATE"])
    expected = pd.DataFrame(columns=["PASS", "RID", "EXAMDATE", "SCANQLTY", "SCANDATE"])

    assert_frame_equal(
        _build_pet_qc_all_studies_for_subject("123_S_1234", df1, df2),
        expected,
        check_index=False,
        check_column_type=False,
    )


def test_build_pet_qc_all_studies_for_subject():
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        _build_pet_qc_all_studies_for_subject,
    )

    df1 = pd.DataFrame(
        {"PASS": [0, 1, 1, 0, 0, 1], "RID": [1234, 1234, 1234, 2345, 2345, 3456]}
    )
    df2 = pd.DataFrame(
        {
            "SCANQLTY": [1, 0, 1, 1],
            "RID": [1234, 1234, 1234, 6789],
            "SCANDATE": ["01/01/2012", "01/06/2012", "01/12/2012", "01/01/2012"],
        }
    )
    expected = pd.DataFrame(
        {
            "PASS": [1, 1, np.nan, np.nan],
            "RID": [1234, 1234, 1234, 1234],
            "EXAMDATE": [np.nan, np.nan, "01/01/2012", "01/12/2012"],
            "SCANQLTY": [np.nan, np.nan, 1, 1],
            "SCANDATE": [np.nan, np.nan, "01/01/2012", "01/12/2012"],
        }
    )

    assert_frame_equal(
        _build_pet_qc_all_studies_for_subject("123_S_1234", df1, df2), expected
    )


def test_compute_fdg_pet_paths(tmp_path, expected_images_df_columns: list[str]):
    from clinica.converters.adni_to_bids.modality_converters._fdg_pet import (
        ADNIPreprocessingStep,
        _compute_fdg_pet_paths,
    )

    csv_dir = tmp_path / "csv"
    csv_dir.mkdir()
    pd.DataFrame(
        {
            "RID": [1, 1, 2, 3],
            "SCANQLTY": [1, 0, 1, 1],
            "SCANDATE": ["01/01/2018", "01/06/2018", "01/01/2018", "06/11/2020"],
            "VISCODE2": ["bl", "m6", "bl", "bl"],
            "LONIUID": ["I1234", "I1234", "I2345", "I3456"],
        }
    ).to_csv(csv_dir / "PETC3.csv")
    pd.DataFrame(
        {
            "RID": [1, 1, 2, 3],
            "PASS": [1, 0, 1, 1],
            "LONIUID": ["I1234", "I1234", "I2345", "I3456"],
        }
    ).to_csv(csv_dir / "PETQC.csv")
    pd.DataFrame(
        {
            "Subject": ["123_S_0001", "123_S_0002", "123_S_0003", "123_S_0004"],
            "Orig/Proc": ["Original"] * 4,
            "Image ID": [10, 11, 12, 13],
            "Scan Date": ["01/01/2018", "01/06/2018", "01/01/2018", "06/11/2020"],
            "Sequence": ["ADNI Brain PET: Raw"] * 4,
        }
    ).to_csv(csv_dir / "PET_META_LIST.csv")

    images = _compute_fdg_pet_paths(
        tmp_path,
        csv_dir,
        ["123_S_0001"],
        tmp_path,
        ADNIPreprocessingStep.from_step_value(2),
    )

    assert len(images) == 0
    assert images.columns.tolist() == expected_images_df_columns
