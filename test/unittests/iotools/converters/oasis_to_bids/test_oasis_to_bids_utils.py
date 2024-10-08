from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_utils import (
    create_sessions_dict,
    write_sessions_tsv,
)


def build_clinical_data(tmp_path: Path) -> None:
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
    df.to_csv(tmp_path / "oasis_cross-sectional.csv", index=False)

    # todo : do excel for future


def build_bids_dir(tmp_path: Path):
    (tmp_path / "sub-OASIS10001" / "ses-M000").mkdir(parents=True)
    (tmp_path / "sub-OASIS10001" / "ses-M006").mkdir(parents=True)
    (tmp_path / "sub-OASIS10002" / "ses-M000").mkdir(parents=True)


def test_create_sessions_dict(tmp_path):
    # todo : how does it handle nan inside excel/csv ? verify with excel
    build_clinical_data(tmp_path)
    build_bids_dir(tmp_path)

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

    result = create_sessions_dict(
        tmp_path,
        tmp_path,
        Path("clinica/iotools/converters/specifications"),
        ["sub-OASIS10001", "sub-OASIS10002"],
    )

    assert result == expected


def test_write_sessions_tsv(tmp_path):
    build_clinical_data(tmp_path)
    build_bids_dir(tmp_path)

    sessions = create_sessions_dict(
        tmp_path,
        tmp_path,
        Path("clinica/iotools/converters/specifications"),
        ["sub-OASIS10001", "sub-OASIS10002"],
    )

    write_sessions_tsv(tmp_path, sessions)

    sessions_files = list(tmp_path.rglob("*.tsv"))

    # todo : finish test
