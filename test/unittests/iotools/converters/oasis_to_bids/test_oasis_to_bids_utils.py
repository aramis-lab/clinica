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
def build_clinical_data(tmp_path: Path) -> None:
    (tmp_path / "clinical").mkdir()

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
    df.to_csv(tmp_path / "clinical" / "oasis_cross-sectional.csv", index=False)

    # todo : future with excel


@pytest.fixture
def build_spec_sessions(tmp_path: Path) -> None:
    (tmp_path / "spec").mkdir()
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
    spec.to_csv(tmp_path / "spec" / "sessions.tsv", index=False, sep="\t")


@pytest.fixture
def build_bids_dir(tmp_path: Path) -> None:
    (tmp_path / "BIDS" / "sub-OASIS10001" / "ses-M000").mkdir(parents=True)
    (tmp_path / "BIDS" / "sub-OASIS10001" / "ses-M006").mkdir(parents=True)
    (tmp_path / "BIDS" / "sub-OASIS10002" / "ses-M000").mkdir(parents=True)


@pytest.fixture
def get_expected_dict() -> dict:
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


def test_create_sessions_dict(
    tmp_path,
    build_clinical_data,
    build_bids_dir,
    build_spec_sessions,
    get_expected_dict,
):
    # todo : how does it handle nan inside excel/csv ? verify with excel

    result = create_sessions_dict(
        tmp_path / "clinical",
        tmp_path / "BIDS",
        tmp_path / "spec",
        ["sub-OASIS10001", "sub-OASIS10002"],
    )

    assert result == get_expected_dict


def test_write_sessions_tsv(
    tmp_path,
    build_clinical_data,
    build_bids_dir,
    build_spec_sessions,
    get_expected_dict,
):
    sessions = create_sessions_dict(
        tmp_path / "clinical",
        tmp_path / "BIDS",
        tmp_path / "spec",
        ["sub-OASIS10001", "sub-OASIS10002"],
    )
    write_sessions_tsv(tmp_path / "BIDS", sessions)
    sessions_files = list((tmp_path / "BIDS").rglob("*.tsv"))
    assert len(sessions_files) == 2
    for file in sessions_files:
        assert not assert_frame_equal(
            pd.read_csv(file, sep="\t").set_index("session_id", drop=False),
            pd.DataFrame(get_expected_dict[file.parent.name]).T.set_index(
                "session_id", drop=False
            ),
            check_like=True,
            check_dtype=False,
        )
        assert file.name == f"{file.parent.name}_sessions.tsv"
