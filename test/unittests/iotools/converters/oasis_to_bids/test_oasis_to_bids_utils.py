from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_utils import (
    create_sessions_dict,
    write_sessions_tsv,
)


@pytest.fixture
def clinical_data_path(tmp_path: Path) -> Path:
    clinical_data_path = tmp_path / "clinical"
    _build_clinical_data(clinical_data_path)
    return clinical_data_path


def _build_clinical_data(clinical_data_path: Path) -> None:
    clinical_data_path.mkdir()

    df = pd.DataFrame(
        {
            "ID": ["OAS1_0001_MR1", "OAS1_0002_MR1"],
            "M/F": ["F", "M"],
            "Hand": ["R", "L"],
            "Age": [74, 67],
            "Educ": [2, 2],
            "SES": [3, 3],
            "MMSE": [29, 29],
            "CDR": [0, 0],
            "eTIV": [1344, 1344],
            "nWBV": [0.704, 0.645],
            "ASF": [1.306, 1.100],
            "Delay": [float("nan"), float("nan")],
        }
    )
    df.to_csv(clinical_data_path / "oasis_cross-sectional.csv", index=False)

    # todo : future with excel


@pytest.fixture
def sessions_path_success(tmp_path: Path) -> Path:
    sessions_path_success = tmp_path / "spec"
    _build_spec_sessions_success(sessions_path_success)
    return sessions_path_success


def _build_spec_sessions_success(sessions_path_success: Path) -> None:
    sessions_path_success.mkdir()
    spec = pd.DataFrame(
        {
            "BIDS CLINICA": ["cdr_global", "MMS", "diagnosis", "foo"],
            "ADNI": [np.nan, np.nan, np.nan, "foo"],
            "OASIS": ["CDR", "MMSE", "CDR", np.nan],
            "OASIS location": [
                "oasis_cross-sectional.csv",
                "oasis_cross-sectional.csv",
                "oasis_cross-sectional.csv",
                np.nan,
            ],
        }
    )
    spec.to_csv(sessions_path_success / "sessions.tsv", index=False, sep="\t")


@pytest.fixture
def sessions_path_error(tmp_path: Path) -> Path:
    sessions_path_error = tmp_path / "spec"
    _build_spec_sessions_error(sessions_path_error)
    return sessions_path_error


def _build_spec_sessions_error(sessions_path_error: Path) -> None:
    sessions_path_error.mkdir()
    spec = pd.DataFrame(
        {
            "BIDS CLINICA": ["foo"],
            "OASIS": ["foo"],
            "OASIS location": [
                "foo.csv",
            ],
        }
    )
    spec.to_csv(sessions_path_error / "sessions.tsv", index=False, sep="\t")


@pytest.fixture
def bids_dir(tmp_path: Path) -> Path:
    bids_dir = tmp_path / "BIDS"
    _build_bids_dir(bids_dir)
    return bids_dir


def _build_bids_dir(bids_dir: Path) -> None:
    (bids_dir / "sub-OASIS10001" / "ses-M000").mkdir(parents=True)
    (bids_dir / "sub-OASIS10001" / "ses-M006").mkdir(parents=True)
    (bids_dir / "sub-OASIS10002" / "ses-M000").mkdir(parents=True)


@pytest.fixture
def expected_dict() -> dict:
    expected = {
        "sub-OASIS10001": {
            "M000": {
                "session_id": "ses-M000",
                "cdr_global": 0,
                "MMS": 29,
                "diagnosis": 0,
            },
            "M006": {
                "session_id": "ses-M006",
                "cdr_global": 0,
                "MMS": 29,
                "diagnosis": 0,
            },
        },
        "sub-OASIS10002": {
            "M000": {
                "session_id": "ses-M000",
                "cdr_global": 0,
                "MMS": 29,
                "diagnosis": 0,
            }
        },
    }

    return expected


def test_create_sessions_dict_success(
    tmp_path,
    clinical_data_path,
    bids_dir,
    sessions_path_success,
    expected_dict,
):
    # todo : how does it handle nan inside excel/csv ? verify with excel

    result = create_sessions_dict(
        clinical_data_path,
        bids_dir,
        sessions_path_success,
        ["sub-OASIS10001", "sub-OASIS10002"],
    )

    assert result == expected_dict


def test_create_sessions_dict_error(
    tmp_path,
    clinical_data_path,
    bids_dir,
    sessions_path_error,
    expected_dict,
):
    # todo : how does it handle nan inside excel/csv ? verify with excel

    with pytest.raises(FileNotFoundError):
        create_sessions_dict(
            clinical_data_path,
            bids_dir,
            sessions_path_error,
            ["sub-OASIS10001", "sub-OASIS10002"],
        )


def test_write_sessions_tsv(
    tmp_path,
    clinical_data_path,
    bids_dir,
    sessions_path_success,
    expected_dict,
):
    sessions = create_sessions_dict(
        clinical_data_path,
        bids_dir,
        sessions_path_success,
        ["sub-OASIS10001", "sub-OASIS10002"],
    )
    write_sessions_tsv(tmp_path / "BIDS", sessions)
    sessions_files = list((tmp_path / "BIDS").rglob("*.tsv"))
    assert len(sessions_files) == 2
    for file in sessions_files:
        assert_frame_equal(
            pd.read_csv(file, sep="\t").set_index("session_id", drop=False),
            pd.DataFrame(expected_dict[file.parent.name]).T.set_index(
                "session_id", drop=False
            ),
            check_like=True,
            check_dtype=False,
        )
        assert file.name == f"{file.parent.name}_sessions.tsv"
