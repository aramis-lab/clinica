from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.mark.parametrize(
    "visit, visit_list, date_list, expected",
    [
        ("bl", ["bl", "m10"], ["01/01/2000", -4], None),
        ("m10", ["bl", "m10"], ["01/01/2000", -4], "11/01/2000"),
        ("m0006", ["bl"], ["01/01/2000"], "07/01/2000"),
    ],
)
def test_compute_exam_date_from_baseline_success(
    visit, date_list, visit_list, expected
):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _compute_exam_date_from_baseline,
    )

    assert _compute_exam_date_from_baseline(visit, date_list, visit_list) == expected


def test_compute_exam_date_from_baseline_raiseValue():
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _compute_exam_date_from_baseline,
    )

    with pytest.raises(
        ValueError,
        match=f"Unexpected visit code foo. Should be in format mX :"
        "Ex: m0, m6, m12, m048...",
    ):
        _compute_exam_date_from_baseline("foo", [], [])


def test_get_csv_files():
    pass


def test_find_exam_date_in_other_csv_files():
    pass


def test_clean_exam_dates():
    pass


def test_load_specifications_success(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _load_specifications,
    )

    filename = "foo.tsv"
    file = pd.DataFrame(columns=["foo"])
    file.to_csv(tmp_path / filename, sep="\t", index=False)
    assert _load_specifications(tmp_path, filename).equals(file)


def test_load_specifications_error_tmp_path(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _load_specifications,
    )

    with pytest.raises(
        FileNotFoundError,
        match=f"The specifications for bar.tsv were not found. "
        f"The should be located in {tmp_path/'bar.tsv'}.",
    ):
        _load_specifications(tmp_path, "bar.tsv")


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


@pytest.mark.parametrize(
    "birth_date, exam_date, age",
    [
        (
            "/2000",
            ["01/02/2000", "02/01/2000", "01/01/2001", "07/06/2003"],
            [0, 0, 1, 3],
        ),
        ("/2001", ["12/30/2003"], [2]),
    ],
)
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
            "VISCODE": ["bl", "bl", "bl", "bl", "m12"],
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
        {
            "RID": [1, 2, 12, 101],
            "VISCODE": ["bl", "bl", "bl", "bl"],
            "PTDOB": ["/1901", "/1902", "/1912", "/2001"],
        }
    )
    ptdemog.to_csv(data_path / "aibl_ptdemog_230ct2024.csv", index=False)

    cdr = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12"],
            "CDGLOBAL": [-4, 1, 0.5, 0, 0],
        }
    )  # rq:float
    cdr.to_csv(data_path / "aibl_cdr_230ct2024.csv", index=False)

    mmse = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12"],
            "MMSCORE": [-4, 10, 10, 30, 29],
        }
    )  # rq:int
    mmse.to_csv(data_path / "aibl_mmse_230ct2024.csv", index=False)

    pdx = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12"],
            "DXCURREN": [-4, 0, 0, 1, 3],
        }
    )  # rq : int
    pdx.to_csv(data_path / "aibl_pdxconv_230ct2024.csv", index=False)

    return data_path


def test_create_sessions_tsv(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        create_sessions_tsv_file,
    )

    bids_path = build_bids_dir(tmp_path)

    create_sessions_tsv_file(
        input_path=bids_path,
        clinical_data_dir=build_clinical_data(tmp_path),
        clinical_specifications_folder=build_sessions_spec(tmp_path),
    )
    result_sub100_list = list(bids_path.rglob("*sub-AIBL100_sessions.tsv"))
    result_sub1_list = list(bids_path.rglob("*sub-AIBL1_sessions.tsv"))

    assert len(result_sub100_list) == 1
    assert len(result_sub1_list) == 1

    result_sub100 = pd.read_csv(result_sub100_list[0], sep="\t")
    result_sub1 = pd.read_csv(result_sub1_list[0], sep="\t")

    expected_sub100 = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M012"],
            "months": [0, 12],
            "age": [np.nan, np.nan],
            "MMS": [30, 29],
            "cdr_global": [0.0, 0.0],
            "diagnosis": ["CN", "AD"],
            "examination_date": ["01/01/2100", "12/01/2100"],
        }
    )

    expected_sub1 = pd.DataFrame(
        {
            "session_id": ["ses-M000"],
            "months": [0],
            "age": [100],
            "MMS": [np.nan],
            "cdr_global": [np.nan],
            "diagnosis": [np.nan],
            "examination_date": ["01/01/2001"],
        }
    )

    assert expected_sub1.equals(result_sub1)
    assert expected_sub100.equals(result_sub100)
