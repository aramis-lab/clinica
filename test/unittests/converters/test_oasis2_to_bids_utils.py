from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


def _build_clinical_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Subject ID": ["OAS2_0001", "OAS2_0001", "OAS2_0002", "OAS2_0003"],
            "MRI ID": [
                "OAS2_0001_MR1",
                "OAS2_0001_MR2",
                "OAS2_0002_MR1",
                "OAS2_0003_MR1",
            ],
            "Group": ["Demented", "Demented", "Nondemented", "Nondemented"],
            "Visit": [1, 2, 1, 1],
            "MR Delay": [0, 432, 0, 0],
            "M/F": ["M", "M", "F", "F"],
            "Hand": ["R", "R", "L", "L"],
            "Age": [54, 54, 81, 81],
            "EDUC": [12, 12, 14, 14],
            "SES": [np.nan, np.nan, 2.0, 2.0],
            "MMSE": [30, 32, 2.0, 2.0],
            "CDR": [0.5, 0.5, 0, 0],
            "eTIV": [1244.4832, 1234.2344, 1100.0543, 1100.0543],
            "nWBV": [0.7322, 0.8653, 0.8543, 0.8543],
            "ASF": [1.0543, 1.8576, 1.2345, 1.2345],
        }
    )


def _build_imaging_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_path": [
                Path("OAS2_0001_MR1/RAW/mpr-1.nifti.img"),
                Path("OAS2_0001_MR1/RAW/mpr-2.nifti.img"),
                Path("OAS2_0001_MR2/RAW/mpr-1.nifti.img"),
                Path("OAS2_0002_MR1/RAW/mpr-1.nifti.img"),
            ],
            "participant_id": [
                "sub-OAS20001",
                "sub-OAS20001",
                "sub-OAS20001",
                "sub-OAS20002",
            ],
            "session_id": ["ses-MR1", "ses-MR1", "ses-MR2", "ses-MR1"],
            "run_number": ["01", "02", "01", "01"],
            "MRI ID": [
                "OAS2_0001_MR1",
                "OAS2_0001_MR1",
                "OAS2_0001_MR2",
                "OAS2_0002_MR1",
            ],
        }
    )


def _build_merged_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_path": [
                Path("OAS2_0001_MR1/RAW/mpr-1.nifti.img"),
                Path("OAS2_0001_MR1/RAW/mpr-2.nifti.img"),
                Path("OAS2_0001_MR2/RAW/mpr-1.nifti.img"),
                Path("OAS2_0002_MR1/RAW/mpr-1.nifti.img"),
            ],
            "participant_id": [
                "sub-OAS20001",
                "sub-OAS20001",
                "sub-OAS20001",
                "sub-OAS20002",
            ],
            "session_id": ["ses-MR1", "ses-MR1", "ses-MR2", "ses-MR1"],
            "run_number": ["01", "02", "01", "01"],
            "MRI ID": [
                "OAS2_0001_MR1",
                "OAS2_0001_MR1",
                "OAS2_0001_MR2",
                "OAS2_0002_MR1",
            ],
            "Subject ID": ["OAS2_0001", "OAS2_0001", "OAS2_0001", "OAS2_0002"],
            "Group": ["Demented", "Demented", "Demented", "Nondemented"],
            "Visit": [1, 1, 2, 1],
            "MR Delay": [0, 0, 432, 0],
            "M/F": ["M", "M", "M", "F"],
            "Hand": ["R", "R", "R", "L"],
            "Age": [54, 54, 54, 81],
            "EDUC": [12, 12, 12, 14],
            "SES": [np.nan, np.nan, np.nan, 2.0],
            "MMSE": [30, 30, 32, 2.0],
            "CDR": [0.5, 0.5, 0.5, 0],
            "eTIV": [1244.4832, 1244.4832, 1234.2344, 1100.0543],
            "nWBV": [0.7322, 0.7322, 0.8653, 0.8543],
            "ASF": [1.0543, 1.0543, 1.8576, 1.2345],
        }
    )


def _build_raw_data(tmp_path: Path) -> Path:
    raw_dataset_path = tmp_path / "raw_data"
    for subject in ("OAS2_0001_MR1", "OAS2_0001_MR2", "OAS2_0002_MR1"):
        image_folder = raw_dataset_path / subject / "RAW"
        image_folder.mkdir(parents=True)
        (image_folder / "mpr-1.nifti.img").touch()
    (raw_dataset_path / "OAS2_0001_MR1" / "RAW" / "mpr-2.nifti.img").touch()
    return raw_dataset_path


def test_read_clinical_data(tmp_path):
    from clinica.converters.oasis2_to_bids._utils import read_clinical_data

    expected_clinical_data = _build_clinical_data()
    expected_clinical_data.to_excel(tmp_path / "clinical_data.xlsx", index=False)

    assert_frame_equal(
        expected_clinical_data,
        read_clinical_data(tmp_path),
        check_like=True,
        check_dtype=False,
    )


def test_read_imaging_data(tmp_path):
    from clinica.converters.oasis2_to_bids._utils import read_imaging_data

    assert_frame_equal(
        read_imaging_data(_build_raw_data(tmp_path)),
        _build_imaging_data(),
        check_like=True,
    )


def test_find_imaging_data(tmp_path):
    from clinica.converters.oasis2_to_bids._utils import _find_imaging_data

    assert set(_find_imaging_data(_build_raw_data(tmp_path))) == {
        Path("OAS2_0001_MR1/RAW/mpr-1.nifti.img"),
        Path("OAS2_0001_MR1/RAW/mpr-2.nifti.img"),
        Path("OAS2_0001_MR2/RAW/mpr-1.nifti.img"),
        Path("OAS2_0002_MR1/RAW/mpr-1.nifti.img"),
    }


@pytest.mark.parametrize(
    "image_file_name, expected",
    [("mpr-1.nifti.img", "01"), ("mpr-4.nifti.img", "04"), ("foo.nifti.img", "01")],
)
def test_identify_run(image_file_name, expected):
    from clinica.converters.oasis2_to_bids._utils import _identify_run

    assert _identify_run(image_file_name) == expected


def test_intersect_data():
    from clinica.converters.oasis2_to_bids._utils import intersect_data

    result = intersect_data(
        df_imaging=_build_imaging_data(), df_clinical=_build_clinical_data()
    )

    expected_merged_data = _build_merged_data()
    expected_merged_data["filename"] = expected_merged_data.apply(
        lambda x: f"anat/{x.participant_id}_{x.session_id}_run-{x.run_number}_T1w.nii.gz",
        axis=1,
    )

    assert_frame_equal(result, expected_merged_data, check_like=True)


def test_build_participants_df():
    from clinica.converters.oasis2_to_bids._utils import _build_participants_df

    expected = pd.DataFrame(
        {
            "participant_id": ["sub-OAS20001", "sub-OAS20002"],
            "sex": ["M", "F"],
            "handedness": ["R", "L"],
            "education": [12, 14],
            "socioeconomic_status": [np.nan, 2.0],
        }
    ).set_index("participant_id")

    assert_frame_equal(
        expected, _build_participants_df(_build_merged_data()), check_like=True
    )


def test_build_sessions_df():
    from clinica.converters.oasis2_to_bids._utils import _build_sessions_df

    expected = pd.DataFrame(
        {
            "participant_id": ["sub-OAS20001", "sub-OAS20001", "sub-OAS20002"],
            "session_id": ["ses-MR1", "ses-MR2", "ses-MR1"],
            "visit": [1, 2, 1],
            "group": ["Demented", "Demented", "Nondemented"],
            "days_since_first_visit": [0, 432, 0],
            "age": [54, 54, 81],
            "mmse": [30, 32, 2.0],
            "cdr": [0.5, 0.5, 0],
            "etiv": [1244.4832, 1234.2344, 1100.0543],
            "nwbv": [0.7322, 0.8653, 0.8543],
            "asf": [1.0543, 1.8576, 1.2345],
        }
    ).set_index(["participant_id", "session_id"])

    assert_frame_equal(
        expected, _build_sessions_df(_build_merged_data()), check_like=True
    )
