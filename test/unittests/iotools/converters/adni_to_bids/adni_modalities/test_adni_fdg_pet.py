from typing import List

import pandas as pd
import pytest


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
