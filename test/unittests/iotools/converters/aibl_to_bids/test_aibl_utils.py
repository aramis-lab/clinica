from pathlib import Path

import pandas as pd
import pytest


def test_listdir_nohidden(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import _listdir_nohidden

    (tmp_path / "file").touch()
    (tmp_path / ".hidden_file").touch()
    (tmp_path / "dir").mkdir()
    (tmp_path / ".hidden_dir").mkdir()

    assert _listdir_nohidden(tmp_path) == ["dir"]


def test_get_first_file_matching_pattern(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import (
        _get_first_file_matching_pattern,
    )

    (tmp_path / "foo_bar.nii.gz").touch()
    (tmp_path / "foo_bar_baz.nii").touch()

    assert (
        _get_first_file_matching_pattern(tmp_path, "foo*.nii*")
        == tmp_path / "foo_bar.nii.gz"
    )


@pytest.mark.parametrize(
    "pattern,msg",
    [("", "Pattern is not valid."), ("foo*.txt*", "No file matching pattern")],
)
def test_get_first_file_matching_pattern_error(tmp_path, pattern, msg):
    from clinica.iotools.converters.aibl_to_bids.utils.bids import (
        _get_first_file_matching_pattern,
    )

    (tmp_path / "foo_bar.nii.gz").touch()
    (tmp_path / "foo_bar_baz.nii").touch()

    with pytest.raises(ValueError, match=msg):
        _get_first_file_matching_pattern(tmp_path, pattern)


@pytest.mark.parametrize("birth_date, exam_date, age", [("", [""], []), ("", [""], [])])
def test_compute_age(birth_date, exam_date, age):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _compute_ages_at_each_exam,
    )

    assert _compute_ages_at_each_exam(birth_date, exam_date) == age


def build_sessions_spec(tmp_path: Path) -> Path:
    spec = pd.DataFrame(
        {
            "BIDS CLINICA": [
                "examination_date",
                "age",
                "cdr_global",
                "MMS",
                "diagnosis",
            ],
            "AIBL": ["EXAMDATE", "PTDOB", "CDGLOBAL", "MMSCORE", "DXCURREN"],
            "AIBL location": [
                "aibl_neurobat_*.csv",
                "aibl_ptdemog_*.csv",
                "aibl_cdr_*.csv",
                "aibl_mmse_*.csv",
                "aibl_pdxconv_*.csv",
            ],
        }
    )
    spec.to_csv(tmp_path / "sessions.tsv", index=False, sep="\t")
    return tmp_path


def build_bids_dir(tmp_path: Path) -> Path:
    bids_dir = tmp_path / "BIDS"
    bids_dir.mkdir()
    (bids_dir / "sub-AIBL1" / "ses-M000").mkdir(parents=True)
    (bids_dir / "sub-AIBL100" / "ses-M000").mkdir(parents=True)
    (bids_dir / "sub-AIBL100" / "ses-M012").mkdir(parents=True)
    return bids_dir


def build_clinical_data(tmp_path: Path) -> Path:
    data_path = tmp_path / "clinical_data"
    data_path.mkdir()

    neuro = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100],  # %m/%d/%Y
            "EXAMDATE": [
                "01/01/2001",
                "01/01/2002",
                "01/01/2012",
                "01/01/2100",
                "12/01/2100",
            ],
        }
    )
    neuro.to_csv(data_path / "aibl_neurobat_230ct2024.csv", index=False)

    ptdemog = pd.DataFrame(
        {"RID": [1, 2, 12, 101], "PTDOB": ["/1901", "/1902", "/1912", "/2001"]}
    )
    ptdemog.to_csv(data_path / "aibl_ptdemog_230ct2024.csv", index=False)

    cdr = pd.DataFrame(
        {"RID": [1, 2, 12, 100, 100], "CDGLOBAL": [-4, 1, 0.5, 0, 0]}
    )  # rq:float
    cdr.to_csv(data_path / "aibl_cdr_230ct2024.csv", index=False)

    mmse = pd.DataFrame(
        {"RID": [1, 2, 12, 100, 100], "MMSCORE": [-4, 10, 10, 30, 29]}
    )  # rq:int
    mmse.to_csv(data_path / "aibl_mmse_230ct2024.csv", index=False)

    pdx = pd.DataFrame(
        {"RID": [1, 2, 12, 100, 100], "DXCURREN": [-4, 0, 0, 1, 3]}
    )  # rq : int
    pdx.to_csv(data_path / "aibl_pdxconv_230ct2024.csv", index=False)

    return data_path


def build_expected_sessions() -> pd.DataFrame:
    expectedsub100 = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M012"],
            "months": ["000", "12"],
            "age": [],
            "MMS": [30, 29],
            "cdr_global": [0, 0],
            "diagnosis": ["CN", "AD"],
            "examination_date": ["01/01/2100", "01/12/2100"],
        }
    )

    expectedsub1 = pd.DataFrame(
        {
            "session_id": ["ses-M000"],
            "months": ["ses-M000"],
            "age": [],
            "MMS": [],
            "cdr_global": [],
            "diagnosis": [],
            "examination_date": [],
        }
    )

    return expectedsub100


def test_create_sessions_tsv(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        create_sessions_tsv_file,
    )

    spec_path = build_sessions_spec(tmp_path)
    clinical_data_path = build_clinical_data(tmp_path)
