from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


@pytest.mark.parametrize(
    "birth_date, exam_date, expected",
    [(None, "foo", None), ("foo", None, None), ("/2000", "01/01/2012", 12)],
)
def test_compute_age_at_exam(birth_date, exam_date, expected):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _compute_age_at_exam,
    )

    assert _compute_age_at_exam(birth_date, exam_date) == expected


@pytest.mark.parametrize(
    "diagnosis, expected",
    [
        (-4, "n/a"),
        (1, "CN"),
        (2, "MCI"),
        (3, "AD"),
        (0, "n/a"),
        (None, "n/a"),
    ],
)
def test_mapping_diagnosis(diagnosis, expected):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _mapping_diagnosis,
    )

    assert _mapping_diagnosis(diagnosis) == expected


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
    (bids_dir / "sub-AIBL109" / "ses-M000").mkdir(parents=True)
    (bids_dir / "sub-AIBL109" / "ses-M006").mkdir(parents=True)
    return bids_dir


def build_clinical_data(tmp_path: Path) -> Path:
    data_path = tmp_path / "clinical_data"
    data_path.mkdir()

    neuro = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100, 109, 109],  # %m/%d/%Y
            "VISCODE": ["bl", "bl", "bl", "bl", "m12", "bl", "m06"],
            "EXAMDATE": [
                "01/01/2001",
                "01/01/2002",
                "01/01/2012",
                "01/01/2100",
                "12/01/2100",
                "01/01/2109",
                -4,
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
            "RID": [1, 2, 12, 100, 100, 109, 109],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12", "bl", "m06"],
            "CDGLOBAL": [-4, 1, 0.5, 0, 0, 0, 0],
            "EXAMDATE": [
                "01/01/2001",
                "01/01/2002",
                "01/01/2012",
                "01/01/2100",
                "12/01/2100",
                "01/01/2109",
                -4,
            ],
        }
    )  # rq:float
    cdr.to_csv(data_path / "aibl_cdr_230ct2024.csv", index=False)

    mmse = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100, 109, 109],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12", "bl", "m06"],
            "MMSCORE": [-4, 10, 10, 30, 29, 10, 10],
            "EXAMDATE": [
                "01/01/2001",
                "01/01/2002",
                "01/01/2012",
                "01/01/2100",
                "12/01/2100",
                "01/01/2109",
                -4,
            ],
        }
    )
    mmse.to_csv(data_path / "aibl_mmse_230ct2024.csv", index=False)

    pdx = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100, 109, 109],
            "VISCODE": ["bl", "bl", "bl", "bl", "m12", "bl", "m06"],
            "DXCURREN": [-4, 0, 0, 1, 3, 2, 2],
        }
    )
    pdx.to_csv(data_path / "aibl_pdxconv_230ct2024.csv", index=False)

    mri3 = pd.DataFrame(
        {
            "RID": [1, 2, 12, 100, 100, 109, 109],  # %m/%d/%Y
            "VISCODE": ["bl", "bl", "bl", "bl", "m12", "bl", "m06"],
            "EXAMDATE": [
                "01/01/2001",
                "01/01/2002",
                "01/01/2012",
                "01/01/2100",
                "12/01/2100",
                "01/01/2109",
                -4,
            ],
        }
    )
    mri3.to_csv(data_path / "aibl_mri3meta_230ct2024.csv", index=False)
    return data_path


def test_extract_metadata_df(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _extract_metadata_df,
    )

    clinical_dir = build_clinical_data(tmp_path)
    expected = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M006"],
            "examination_date": ["01/01/2109", "-4"],
        }
    ).set_index("session_id", drop=True)
    result = _extract_metadata_df(
        pd.read_csv(clinical_dir / "aibl_neurobat_230ct2024.csv", dtype={"text": str}),
        109,
        bids_metadata="examination_date",
        source_metadata="EXAMDATE",
    )

    assert_frame_equal(expected, result)


@pytest.mark.parametrize(
    "source_id, session_id, expected",
    [
        (109, "ses-M000", "01/01/2109"),
        (109, "ses-M006", None),
        (109, "ses-M014", None),
        (0, "ses-M014", None),
    ],
)
def test_find_exam_date_in_other_csv_files(tmp_path, source_id, session_id, expected):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _find_exam_date_in_other_csv_files,
    )

    clinical_dir = build_clinical_data(tmp_path)
    assert (
        _find_exam_date_in_other_csv_files(source_id, session_id, clinical_dir)
        == expected
    )


def test_get_csv_paths(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import _get_csv_paths

    assert _get_csv_paths(tmp_path) == ()

    clinical_dir = build_clinical_data(tmp_path)
    csv_paths = [Path(path).name for path in _get_csv_paths(clinical_dir)]

    assert set(csv_paths) == {
        "aibl_cdr_230ct2024.csv",
        "aibl_mmse_230ct2024.csv",
        "aibl_mri3meta_230ct2024.csv",
    }


@pytest.mark.parametrize(
    "rid, session_id, exam_date, expected",
    [
        (0, "foo", "01/01/2000", "01/01/2000"),
        (0, None, "01/01/2000", "01/01/2000"),
        (0, None, None, None),
        (109, "ses-M000", None, "01/01/2109"),
        (109, "ses-M006", None, None),
    ],
)
def test_complete_examination_dates(tmp_path, rid, session_id, exam_date, expected):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        _complete_examination_dates,
    )

    clinical_dir = build_clinical_data(tmp_path)
    assert (
        _complete_examination_dates(rid, session_id, exam_date, clinical_dir)
        == expected
    )


def test_create_sessions_tsv(tmp_path):
    from clinica.iotools.converters.aibl_to_bids.utils.clinical import (
        create_sessions_tsv_file,
    )

    bids_path = build_bids_dir(tmp_path)

    create_sessions_tsv_file(
        bids_dir=bids_path,
        clinical_data_dir=build_clinical_data(tmp_path),
        clinical_specifications_folder=build_sessions_spec(tmp_path),
    )
    result_sub100_list = list(bids_path.rglob("*sub-AIBL100_sessions.tsv"))
    result_sub1_list = list(bids_path.rglob("*sub-AIBL1_sessions.tsv"))
    result_sub109_list = list(bids_path.rglob("*sub-AIBL109_sessions.tsv"))

    assert len(result_sub109_list) == 1
    assert len(result_sub100_list) == 1
    assert len(result_sub1_list) == 1

    result_sub109 = pd.read_csv(result_sub109_list[0], sep="\t", keep_default_na=False)
    result_sub100 = pd.read_csv(result_sub100_list[0], sep="\t", keep_default_na=False)
    result_sub1 = pd.read_csv(result_sub1_list[0], sep="\t", keep_default_na=False)

    expected_sub100 = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M012"],
            "age": ["n/a", "n/a"],
            "MMS": [30, 29],
            "cdr_global": [0.0, 0.0],
            "diagnosis": ["CN", "AD"],
            "examination_date": ["01/01/2100", "12/01/2100"],
        }
    )

    expected_sub1 = pd.DataFrame(
        {
            "session_id": ["ses-M000"],
            "age": [100],
            "MMS": ["n/a"],
            "cdr_global": ["n/a"],
            "diagnosis": ["n/a"],
            "examination_date": ["01/01/2001"],
        }
    )

    expected_sub109 = pd.DataFrame(
        {
            "session_id": ["ses-M000", "ses-M006"],
            "age": ["n/a", "n/a"],
            "MMS": [10, 10],
            "cdr_global": [0.0, 0.0],
            "diagnosis": ["MCI", "MCI"],
            "examination_date": ["01/01/2109", "n/a"],
        }
    )

    assert_frame_equal(result_sub1, expected_sub1, check_like=True)
    assert_frame_equal(result_sub100, expected_sub100, check_like=True)
    assert_frame_equal(result_sub109, expected_sub109, check_like=True)
