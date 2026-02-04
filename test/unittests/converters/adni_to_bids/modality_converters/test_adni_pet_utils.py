from pathlib import Path

import pandas as pd
from pandas.testing import assert_frame_equal


def _build_all_images_df(csv_path: Path) -> Path:
    all_images_df = pd.DataFrame(
        {
            "image_id": [10, 11, 12, 13],
            "subject_id": [
                "123_S_0001",
                "123_S_0002",
                "123_S_0003",
                "123_S_0004",
            ],
            "phase": ["ADNI 1", "ADNI 1", "ADNI 1", "ADNI 1"],
            "study_id": [2815] * 4,
            "image_type": ["Original"] * 4,
            "image_date": [
                "2018-01-01",
                "2018-01-06",
                "2018-01-01",
                "2020-06-11",
            ],
            "image_description": [
                "ADNI Brain PET: Raw",
                "Co-registered, Averaged",
                "Co-registered Dynamic",
                "Coreg, Avg, Standardized Image and Voxel Size",
            ],
        }
    )

    file_path = csv_path / "cohort_1_All_Images_27Jan2026.csv"

    all_images_df.to_csv(file_path, index=False)

    return file_path


def _build_manifest_df(csv_path: Path) -> Path:
    manifest_df = pd.DataFrame(
        {
            "image_id": [10, 11, 12, 13],
            "series_id": [10601, 10602, 10603, 10604],
        }
    )

    file_path = csv_path / "cohort_1_Manifest_27Jan2026.csv"

    manifest_df.to_csv(file_path, index=False)

    return file_path


def _build_expected_df() -> pd.DataFrame:
    expected_df = pd.DataFrame(
        {
            "image_id": [10, 11, 12, 13],
            "subject_id": [
                "123_S_0001",
                "123_S_0002",
                "123_S_0003",
                "123_S_0004",
            ],
            "phase": ["ADNI 1", "ADNI 1", "ADNI 1", "ADNI 1"],
            "study_id": [2815] * 4,
            "image_type": ["Original"] * 4,
            "image_date": [
                "2018-01-01",
                "2018-01-06",
                "2018-01-01",
                "2020-06-11",
            ],
            "image_description": [
                "ADNI Brain PET: Raw",
                "Co-registered, Averaged",
                "Co-registered Dynamic",
                "Coreg, Avg, Standardized Image and Voxel Size",
            ],
            "series_id": [10601, 10602, 10603, 10604],
        }
    )

    return expected_df


def test_load_all_images_metadata(tmp_path):
    from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
        load_all_images_metadata,
    )

    _build_all_images_df(tmp_path)
    _build_manifest_df(tmp_path)
    expected_df = _build_expected_df()

    result = load_all_images_metadata(tmp_path)

    assert_frame_equal(result, expected_df, check_dtype=False)
