import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


def test_load_all_images_metadata(tmp_path):
    from clinica.converters.adni_to_bids.modality_converters._pet_utils import (
        load_all_images_metadata,
    )

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
                "Co-registered, Avereged",
                "Co-registered Dynamic",
                "Coreg, Avg, Standardized Image and Voxel Size",
            ],
        }
    )

    manifest_df = pd.DataFrame(
        {
            "image_id": [10, 11, 12, 13],
            "series_id": [10601, 10602, 10603, 10604],
        }
    )

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
                "Co-registered, Avereged",
                "Co-registered Dynamic",
                "Coreg, Avg, Standardized Image and Voxel Size",
            ],
            "series_id": [10601, 10602, 10603, 10604],
        }
    )

    all_images_df.to_csv(tmp_path / "All_Images.csv", index=False)
    manifest_df.to_csv(tmp_path / "Manifest.csv", index=False)

    result = load_all_images_metadata(tmp_path)

    assert_frame_equal(result, expected_df, check_dtype=False)
