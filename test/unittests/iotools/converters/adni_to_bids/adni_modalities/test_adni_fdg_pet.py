from typing import List

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


@pytest.mark.parametrize("step_value,expected", [(2, "fdg"), (4, "fdg_uniform")])
def test_get_modality_from_adni_preprocessing_step(step_value, expected):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
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
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
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
def expected_fdg_df_columns() -> List[str]:
    return [
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
        "Is_Dicom",
        "Path",
    ]


def test_compute_fdg_pet_paths_empty(tmp_path, expected_fdg_df_columns):
    """Checks that _compute_fdg_pet_paths returns an empty dataframe with the right
    columns when the provided list of subjects is empty.
    """
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
        ADNIPreprocessingStep,
        _compute_fdg_pet_paths,
    )

    images = _compute_fdg_pet_paths(
        tmp_path, tmp_path, [], tmp_path, ADNIPreprocessingStep.from_step_value(2)
    )
    assert len(images) == 0
    assert images.columns.tolist() == expected_fdg_df_columns


def test_compute_fdg_pet_paths_column_errors(tmp_path, expected_fdg_df_columns):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
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
def test_is_visit_a_conversion_error(row, expected):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
        _is_visit_a_conversion_error,
    )

    assert _is_visit_a_conversion_error(row) is expected


def test_remove_known_conversion_errors():
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
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
        _remove_known_conversion_errors(input_df).reset_index(drop=True),
        expected_df.reset_index(drop=True),
    )


@pytest.mark.skip(reason="Test is not finished...")
def test_compute_fdg_pet_paths(tmp_path, expected_fdg_df_columns):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fdg_pet import (
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
        }
    ).to_csv(csv_dir / "PETC3.csv")
    pd.DataFrame({"RID": [1, 1, 2, 3], "PASS": [1, 0, 1, 1]}).to_csv(
        csv_dir / "PETQC.csv"
    )
    pd.DataFrame({"Subject": ["sub-0001", "sub-0001", "sub-0002", "sub-0003"]}).to_csv(
        csv_dir / "PET_META_LIST.csv"
    )
    images = _compute_fdg_pet_paths(
        tmp_path,
        csv_dir,
        ["sub-0001"],
        tmp_path,
        ADNIPreprocessingStep.from_step_value(2),
    )
    assert len(images) == 0
    assert images.columns.tolist() == expected_fdg_df_columns
