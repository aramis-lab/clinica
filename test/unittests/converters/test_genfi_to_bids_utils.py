import math
import re

import numpy as np
import pandas as pd
import pytest
from fsspec.implementations.local import LocalFileSystem
from pandas.testing import assert_frame_equal


@pytest.mark.parametrize(
    "input,expected",
    [
        ("DWI", "dwi"),
        ("T1", "T1"),
        ("T2", "T2w"),
        ("fieldmap", "fieldmap"),
        ("fmri", "rsfmri"),
    ],
)
def test_identify_modality(input: str, expected: str):
    from clinica.converters._utils import identify_modality

    assert identify_modality(input) == expected


def test_identify_modality_is_nan():
    from clinica.converters._utils import identify_modality

    assert math.isnan(identify_modality("blzflbzv"))


@pytest.fixture
def _generate_clinical_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["GRN001", "GRN001", "C9ORF001", "C9ORF001"],
            "source_ses_id": ["01", "11", "17", "21"],
        }
    )


def test_compute_genfi_version(_generate_clinical_data: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_genfi_version

    expected = pd.DataFrame(
        {
            "source_id": ["GRN001", "GRN001", "C9ORF001", "C9ORF001"],
            "source_ses_id": ["01", "11", "17", "21"],
            "genfi_version": ["GENFI1", "GENFI2", "GENFI2", "GENFI3"],
        }
    )
    assert_frame_equal(_compute_genfi_version(_generate_clinical_data), expected)


def test_compute_session_numbers(_generate_clinical_data: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_session_numbers

    expected = pd.DataFrame(
        {
            "source_id": ["GRN001", "GRN001", "C9ORF001", "C9ORF001"],
            "source_ses_id": [1, 11, 17, 21],
            "session_id": ["ses-01", "ses-11", "ses-17", "ses-21"],
        }
    )
    assert_frame_equal(_compute_session_numbers(_generate_clinical_data), expected)


@pytest.fixture
def input_df_compute_runs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["GRN001", "GRN001", "GRN001", "C9ORF001", "C9ORF001"],
            "source_ses_id": [1, 1, 1, 2, 2],
            "suffix": ["a", "a", "c", "a", "b"],
            "number_of_parts": [2, 2, 1, 1, 1],
            "dir_num": [10, 20, 30, 10, 20],
        }
    )


def test_compute_runs(input_df_compute_runs: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_run_numbers_from_parts

    expected = pd.DataFrame(
        {
            "source_id": ["C9ORF001", "C9ORF001", "GRN001", "GRN001", "GRN001"],
            "source_ses_id": [2, 2, 1, 1, 1],
            "suffix": ["a", "b", "a", "a", "c"],
            "number_of_parts": [1, 1, 2, 2, 1],
            "dir_num": [10, 20, 10, 20, 30],
            "run_01_dir_num": [10, 20, 10, 10, 30],
            "run": [False, False, False, True, False],
            "run_number": [1, 1, 1, 2, 1],
            "run_num": ["run-01", "run-01", "run-01", "run-02", "run-01"],
        }
    )
    assert_frame_equal(_compute_run_numbers_from_parts(input_df_compute_runs), expected)


@pytest.fixture
def input_df_compute_philips_parts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["sub-01"] * 5 + ["sub-02"] * 3,
            "source_ses_id": [1, 1, 1, 2, 2, 1, 2, 2],
            "suffix": ["dwi"] * 8,
            "dir_num": [10, 20, 30, 40, 50, 10, 20, 30],
        }
    )


def test_compute_philips_parts(input_df_compute_philips_parts: pd.DataFrame):
    from clinica.converters.genfi_to_bids._utils import _compute_philips_parts

    expected = pd.DataFrame(
        {
            "source_id": ["sub-01"] * 5 + ["sub-02"] * 3,
            "source_ses_id": [1, 1, 1, 2, 2, 1, 2, 2],
            "suffix": ["dwi"] * 8,
            "dir_num": [10, 20, 30, 40, 50, 10, 20, 30],
            "part_01_dir_num": [10, 10, 10, 40, 40, 10, 20, 20],
            "run": [False, True, True, False, True, False, False, True],
            "part_number": [1, 2, 3, 1, 2, 1, 1, 2],
            "number_of_parts": [3, 3, 3, 2, 2, 1, 2, 2],
        }
    )
    assert_frame_equal(_compute_philips_parts(input_df_compute_philips_parts), expected)


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            [False, True, False, True, False, False, False, True, False, True, True],
            [1, 2, 1, 2, 1, 1, 1, 2, 1, 2, 3],
        ),
        ([True, True, True], [1, 2, 3]),
        ([False, False, False], [1, 1, 1]),
    ],
)
def test_compute_scan_sequence_numbers(input: list[bool], expected: list[int]):
    from clinica.converters.genfi_to_bids._utils import _compute_scan_sequence_numbers

    assert _compute_scan_sequence_numbers(input) == expected


def test_compute_scan_sequence_numbers_error():
    from clinica.converters.genfi_to_bids._utils import _compute_scan_sequence_numbers

    with pytest.raises(ValueError):
        _compute_scan_sequence_numbers([])


def test_drop_duplicate_line_with_nans_error():
    from clinica.converters.genfi_to_bids._utils import _drop_duplicate_line_with_nans
    from clinica.utils.exceptions import ClinicaBIDSError

    with pytest.raises(
        ClinicaBIDSError,
        match="Column participant_id was not found in the participants tsv while it is required by BIDS specifications.",
    ):
        _drop_duplicate_line_with_nans(pd.DataFrame())


def test_drop_duplicate_line_with_nans():
    from clinica.converters.genfi_to_bids._utils import _drop_duplicate_line_with_nans

    assert_frame_equal(
        _drop_duplicate_line_with_nans(
            pd.DataFrame(
                {
                    "participant_id": ["sub-1", "sub-1", "sub-2", "sub-3"],
                    "metadata": ["1", np.nan, np.nan, "3"],
                }
            )
        ).reset_index(drop=True),
        pd.DataFrame(
            {
                "participant_id": ["sub-1", "sub-2", "sub-3"],
                "metadata": ["1", np.nan, "3"],
            }
        ),
    )


def test_write_description_and_participants(tmp_path):
    from clinica.converters.genfi_to_bids._utils import (
        _write_description_and_participants,
    )

    participants = pd.DataFrame(
        {
            "participant_id": [
                "sub-C9ORF002",
                "sub-C9ORF002",
                "sub-C9ORF002",
                "sub-MAPT003",
                "sub-MAPT003",
                "sub-MAPT003",
            ],
            "session_id": [
                "ses-01",
                "ses-11",
                np.nan,
                "ses-01",
                "ses-11",
                np.nan,
            ],
            "modality": ["T1", "dwi", np.nan, "T2w", "rsfmri", np.nan],
            "run_num": ["run-01", "run-02", np.nan, "run-03", "run-04", np.nan],
            "bids_filename": [
                "sub-C9ORF002_ses-01_run-01_T1w",
                "sub-C9ORF002_ses-11_run-02_dwi",
                np.nan,
                "sub-MAPT003_ses-02_run-03_T2w",
                "sub-MAPT003_ses-11_task-rest_run-04_bold",
                np.nan,
            ],
            "source": [
                "path/to/file/1",
                "path/to/file/2",
                np.nan,
                "path/to/file/3",
                "path/to/file/4",
                np.nan,
            ],
        }
    )

    fs = LocalFileSystem(auto_mkdir=True)

    _write_description_and_participants(to=tmp_path, fs=fs, participants=participants)

    # Check files existing
    desc_file = tmp_path / "dataset_description.json"
    tsv_file = tmp_path / "participants.tsv"
    assert desc_file.exists()
    assert tsv_file.exists()

    # Check participants.tsv content
    df_out = pd.read_csv(tsv_file, sep="\t")
    assert "participant_id" in df_out.columns
    assert set(df_out["participant_id"]) == {"sub-C9ORF002", "sub-MAPT003"}


def test_write_sessions(tmp_path):
    from clinica.converters.genfi_to_bids._utils import _write_sessions

    participant_id = "sub-C9ORF004"

    (tmp_path / participant_id / "ses-01").mkdir(parents=True)
    (tmp_path / participant_id / "ses-11").mkdir(parents=True)

    sessions = pd.DataFrame(
        {
            "participant_id": [participant_id, participant_id],
            "modality": ["T1w", "dwi"],
            "bids_filename": [
                "sub-C9ORF004_ses-01_run-01_T1w",
                "sub-C9ORF004_ses-11_run-02_dwi",
            ],
            "run_num": ["run-01", "run-02"],
            "session_id": ["ses-01", "ses-11"],
        }
    ).set_index(
        ["participant_id", "modality", "bids_filename", "run_num", "session_id"]
    )

    fs = LocalFileSystem(auto_mkdir=True)

    _write_sessions(to=tmp_path, fs=fs, sessions=sessions)

    # Check files existing
    tsv_path = tmp_path / participant_id / f"{participant_id}_sessions.tsv"
    assert tsv_path.exists()

    # Check sessions.tsv content
    df_out = pd.read_csv(tsv_path, sep="\t")
    assert set(df_out["session_id"]) == {"ses-01", "ses-11"}


def test_write_scans_and_niftis(tmp_path, mocker):
    from clinica.converters.genfi_to_bids._utils import _write_scans_and_niftis

    mock_run = mocker.patch(
        "clinica.converters._utils.run_dcm2niix",
        return_value=True,
    )

    participant_id = "sub-C9ORF006"
    session_id = "ses-01"
    source = tmp_path / "source"
    source.mkdir()

    scans = pd.DataFrame(
        [
            {
                "participant_id": participant_id,
                "session_id": session_id,
                "bids_full_path": f"{participant_id}/{session_id}/anat/{participant_id}_{session_id}_run-01_T1w.nii.gz",
                "bids_filename": f"{participant_id}_{session_id}_run-01_T1w",
                "modality": "T1",
                "run_num": "run-01",
                "source_path": source / "file.dcm",
                "manufacturer": "SIEMENS",
                "number_of_parts": 1,
            }
        ]
    )

    _write_scans_and_niftis(to=tmp_path, source=source, scans=scans)

    # Check run_dcm2niix was called
    assert mock_run.called

    # Check scans.tsv existing
    tsv_path = (
        tmp_path
        / participant_id
        / session_id
        / f"{participant_id}_{session_id}_scans.tsv"
    )
    assert tsv_path.exists()

    # Check scans.tsv content
    df_out = pd.read_csv(tsv_path, sep="\t")

    expected_cols = [
        "index",
        "filename",
        "modality",
        "run_num",
        "source_path",
        "manufacturer",
        "number_of_parts",
    ]

    assert set(df_out.columns) == set(expected_cols)
    assert len(df_out) == 1

    row = df_out.iloc[0]

    assert row["modality"] == "T1"
    assert row["run_num"] == "run-01"
    assert row["filename"] == f"{participant_id}_{session_id}_run-01_T1w"
    assert not str(row["source_path"]).startswith(str(source))
    assert row["manufacturer"] == "SIEMENS"
    assert row["number_of_parts"] == 1


@pytest.mark.parametrize(
    "filename",
    [
        "FINAL_IMAGING.xlsx",
        "FINAL_DEMOGRAPHICS.xlsx",
        "FINAL_CLINICAL.xslx",
        "FINAL_BIOSAMPLES.xslx",
        "FINAL_NEUROPSYCH.xslx",
        "FINAL_GENETICS.xslx",
    ],
)
def test_find_clinical_data(monkeypatch, tmp_path, filename):
    import clinica.converters.genfi_to_bids._utils as genfi_utils

    file_path = tmp_path / filename

    def _check_file(clinical_data_directory, pattern):
        return file_path

    def _read_file(data_file):
        return pd.DataFrame({"source": [data_file.name]})

    monkeypatch.setattr(genfi_utils, "_check_file", _check_file)
    monkeypatch.setattr(genfi_utils, "_read_file", _read_file)

    out = genfi_utils._find_clinical_data(tmp_path, filename)

    assert len(out) == 1
    assert out.iloc[0, 0].endswith(filename)


def test_merge_and_coalesce():
    from clinica.converters.genfi_to_bids._utils import _merge_and_coalesce

    left_df = pd.DataFrame(
        {
            "blinded_code": ["C9ORF001", "C9ORF002"],
            "genetic_status_1": ["P", pd.NA],
            "genetic_status_2": [0, 1],
        }
    )

    right_df = pd.DataFrame(
        {
            "blinded_code": ["C9ORF002", "C9ORF003"],
            "genetic_status_1": ["A", "P"],
            "diagnosis": ["bvFTD", "ALS"],
        }
    )

    on = ["blinded_code"]

    expected_df = pd.DataFrame(
        {
            "blinded_code": ["C9ORF001", "C9ORF002", "C9ORF003"],
            "genetic_status_1": ["P", "A", "P"],
            "genetic_status_2": [0, 1, pd.NA],
            "diagnosis": [pd.NA, "bvFTD", "ALS"],
        }
    )

    result = _merge_and_coalesce(left_df, right_df, on=on)

    assert_frame_equal(result, expected_df, check_like=True, check_dtype=False)


def test_complete_clinical_data():
    from clinica.converters.genfi_to_bids._utils import _complete_clinical_data

    df_imaging = pd.DataFrame(
        {
            "blinded_code": ["C9ORF001", "C9ORF002"],
            "blinded_site": ["GENFI_AA", "GENFI_AA"],
            "visit": [1, 1],
            "scan_for_qc": [1, 1],
            "diagnosis": ["bvFTD", pd.NA],
            "plasma_nfl": [pd.NA, pd.NA],
        }
    )

    df_clinical_list = [
        pd.DataFrame(
            {
                "blinded_code": ["C9ORF001", "C9ORF002"],
                "blinded_site": ["GENFI_AA", "GENFI_AA"],
                "visit": [1, 1],
                "gender": [0, 1],
            }
        ),
        pd.DataFrame(
            {
                "blinded_code": ["C9ORF001", "C9ORF002"],
                "blinded_site": ["GENFI_AA", "GENFI_AA"],
                "visit": [1, 1],
                "diagnosis": ["bvFTD", "ALS"],
            }
        ),
        pd.DataFrame(
            {
                "blinded_code": ["C9ORF001", "C9ORF002"],
                "blinded_site": ["GENFI_AA", "GENFI_AA"],
                "visit": [1, 1],
                "plasma_nfl": [6, 5],
                "diagnosis": [pd.NA, "ALS"],
            }
        ),
        pd.DataFrame(
            {
                "blinded_code": ["C9ORF001", "C9ORF002"],
                "blinded_site": ["GENFI_AA", "GENFI_AA"],
                "visit": [1, 1],
                "ds_f_score": [8, 11],
            }
        ),
        pd.DataFrame(
            {
                "blinded_code": ["C9ORF001", "C9ORF002"],
                "blinded_site": ["GENFI_AA", "GENFI_AA"],
                "visit": [1, 1],
                "genetic_status_1": ["P", "A"],
            }
        ),
    ]

    expected = pd.DataFrame(
        {
            "blinded_code": ["C9ORF001", "C9ORF002"],
            "blinded_site": ["GENFI_AA", "GENFI_AA"],
            "visit": [1, 1],
            "scan_for_qc": [1, 1],
            "gender": [0, 1],
            "diagnosis": ["bvFTD", "ALS"],
            "plasma_nfl": [6, 5],
            "ds_f_score": [8, 11],
            "genetic_status_1": ["P", "A"],
        }
    )

    result = _complete_clinical_data(
        df_imaging=df_imaging, df_clinical_list=df_clinical_list
    )

    assert_frame_equal(result, expected, check_like=True, check_dtype=False)


@pytest.mark.parametrize(
    ("full", "gif", "expected"),
    [
        (False, False, "mandatory_specs"),
        (True, True, "full_specs"),
        (True, False, "full_specs"),
        (False, True, "gif_specs"),
    ],
)
def test_specs_depending_on_option(full, gif, expected):
    from clinica.converters.genfi_to_bids._utils import _specs_depending_on_option

    assert _specs_depending_on_option(full, gif) == expected


TO_COMPLETE_SPECS_DF = pd.DataFrame(
    {
        "participants": ["participant_id", "source", "blinded_code"],
        "sessions": ["participant_id", "session_id", "genfi_version"],
        "scans": ["participant_id", "session_id", "genfi_version"],
    }
)


FULL_SPECS_DF = pd.DataFrame(
    {
        "participants": [
            "participant_id",
            "source",
            "blinded_code",
            "blinded_family",
            "blinded_site",
        ],
        "sessions": ["participant_id", "session_id", "genfi_version", "aad", "aad_1"],
        "scans": [
            "participant_id",
            "session_id",
            "genfi_version",
            "bids_filename",
            "bids_full_path",
        ],
    }
)


@pytest.mark.parametrize(
    ("clinical_data_list", "expected"),
    [
        (
            ["blinded_family\n", "aad_1\n", "bids_filename"],
            ["blinded_family", "aad_1", "bids_filename"],
        ),
        (
            ["  blinded_site  \n", "\n", "aad\n", "\tbids_full_path\n"],
            ["blinded_site", "aad", "bids_full_path"],
        ),  # whitespaces + empty lines
    ],
)
def test_load_clinical_data_list_success(tmp_path, clinical_data_list, expected):
    from clinica.converters.genfi_to_bids._utils import _load_clinical_data_list

    cdt_path = tmp_path / "additional_clinical_data.txt"
    cdt_path.write_text("".join(clinical_data_list), encoding="utf-8")

    out = _load_clinical_data_list(cdt_path, FULL_SPECS_DF)

    assert out == expected


@pytest.mark.parametrize(
    ("clinical_data_list", "expected"),
    [
        ([], "'-clinical_data_txt/cdt' is empty (no valid entries found)."),
        (
            ["\n", "   \n", "\t\n"],
            "'-clinical_data_txt/cdt' is empty (no valid entries found).",
        ),
    ],
)
def test_load_clinical_data_list_empty(tmp_path, clinical_data_list, expected):
    from clinica.converters.genfi_to_bids._utils import _load_clinical_data_list

    cdt_path = tmp_path / "additional_clinical_data.txt"
    cdt_path.write_text("".join(clinical_data_list), encoding="utf-8")

    with pytest.warns(UserWarning, match=re.escape(expected)):
        _load_clinical_data_list(cdt_path, FULL_SPECS_DF)


def test_load_clinical_data_list_unknown_field(tmp_path):
    from clinica.converters.genfi_to_bids._utils import _load_clinical_data_list

    cdt_path = tmp_path / "additional_clinical_data.txt"
    cdt_path.write_text("blinded_family\nfalse_field\naad\n", encoding="utf-8")

    expected = "Line 2: 'false_field' not found in specifications."

    with pytest.warns(UserWarning, match=re.escape(expected)):
        _load_clinical_data_list(cdt_path, FULL_SPECS_DF)


def test_merge_clinical_data_list_into_df_in_matching_columns():
    from clinica.converters.genfi_to_bids._utils import (
        _merge_clinical_data_list_into_df,
    )

    clinical_data_list = ["blinded_family", "aad", "bids_filename", "aad_1"]

    out = _merge_clinical_data_list_into_df(
        clinical_data_list, FULL_SPECS_DF, TO_COMPLETE_SPECS_DF.copy()
    )

    expected = pd.DataFrame(
        {
            "participants": [
                "participant_id",
                "source",
                "blinded_code",
                "blinded_family",
                pd.NA,
            ],
            "sessions": [
                "participant_id",
                "session_id",
                "genfi_version",
                "aad",
                "aad_1",
            ],
            "scans": [
                "participant_id",
                "session_id",
                "genfi_version",
                "bids_filename",
                pd.NA,
            ],
        }
    )

    assert_frame_equal(out, expected)


def test_merge_clinical_data_list_into_df_no_duplicate():
    from clinica.converters.genfi_to_bids._utils import (
        _merge_clinical_data_list_into_df,
    )

    clinical_data_list = ["participant_id", "session_id"]

    out = _merge_clinical_data_list_into_df(
        clinical_data_list, FULL_SPECS_DF, TO_COMPLETE_SPECS_DF.copy()
    )

    expected = pd.DataFrame(
        {
            "participants": ["participant_id", "source", "blinded_code"],
            "sessions": ["participant_id", "session_id", "genfi_version"],
            "scans": ["participant_id", "session_id", "genfi_version"],
        }
    )

    assert_frame_equal(out, expected)
