from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from clinica.converters.oasis3_to_bids._utils import _find_csv_with_filename


def _build_clinical_data(tmp_path: Path) -> Path:
    clinical_data_directory = tmp_path / "clinical_data"
    clinical_data_directory.mkdir()

    udsb1 = pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30002"],
            "days_to_visit": ["0000", "0101"],
            "age at visit": [62.12, 45.67],
            "WEIGHT": [152, 184.2],
        }
    )
    udsb1.to_csv(
        clinical_data_directory / "OASIS3_UDSb1_physical_eval.csv", index=False
    )

    udsb4 = pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30002"],
            "days_to_visit": ["0000", "0101"],
            "age at visit": [62.12, 45.67],
            "memory": [0, 1],
        }
    )
    udsb4.to_csv(clinical_data_directory / "OASIS3_UDSb4_cdr.csv", index=False)

    demo = pd.DataFrame(
        data={
            "OASISID": ["OAS30001", "OAS30002"],
            "AgeAtEntry": [62, 45],
            "GENDER": [2, 1],
        }
    )
    demo.to_csv(clinical_data_directory / "OASIS3_demographics.csv", index=False)

    return clinical_data_directory


def test_find_csv_with_filename_success(tmp_path):
    expected = pd.DataFrame({"test_data": ["foo", "bar", "foobar"]})
    expected.to_csv(tmp_path / "foo.csv", index=False)
    assert_frame_equal(_find_csv_with_filename(tmp_path, filename="foo"), expected)


def test_find_csv_with_filename_no_file(tmp_path):
    with pytest.raises(FileNotFoundError, match=f"No CSV files found"):
        _find_csv_with_filename(tmp_path, "foo")


def test_find_csv_with_filename_more_files(tmp_path):
    (tmp_path / "foo.csv").touch()
    (tmp_path / "dupe_foo.csv").touch()
    with pytest.raises(FileNotFoundError, match="More than one file"):
        _find_csv_with_filename(tmp_path, filename="foo")


def test_read_mri_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _read_mri_data

    pd.DataFrame().to_csv(tmp_path / "OASIS3_MR_json.csv")
    _read_mri_data(tmp_path)


def test_find_imaging_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _find_imaging_data

    (tmp_path / "images").mkdir()
    (tmp_path / "images" / "test.nii.gz").touch()
    assert "images/test.nii.gz" == str(next(_find_imaging_data(tmp_path)))


@pytest.mark.parametrize(
    "input, expected",
    [
        ("sub-01_ses-02_t1.nii.gz", "run-01"),
        ("imaging_data/sub-01_ses-02_run-12_t1.nii.gz", "run-12"),
    ],
)
def test_identify_runs(input, expected):
    from clinica.converters.oasis3_to_bids._utils import _identify_runs

    assert expected == _identify_runs(Path(input))


def test_read_clinical_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _read_clinical_data

    result = _read_clinical_data(_build_clinical_data(tmp_path))
    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30002"],
            "days_to_visit": [0, 101],
            "age at visit": [62.12, 45.67],
            "WEIGHT": [152, 184.2],
            "memory": [0, 1],
            "session_id": ["d0000", "d0101"],
        }
    )
    assert_frame_equal(expected, result, check_exact=False)


def test_read_demo_data(tmp_path):
    from clinica.converters.oasis3_to_bids._utils import _read_demo_data

    expected = pd.DataFrame(
        data={
            "Subject": ["OAS30001", "OAS30002"],
            "AgeAtEntry": [62, 45],
            "sex": ["F", "M"],
        }
    )
    result = _read_demo_data(_build_clinical_data(tmp_path))
    assert_frame_equal(expected, result, check_exact=False)
