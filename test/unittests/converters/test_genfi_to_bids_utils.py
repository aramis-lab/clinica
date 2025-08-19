import math

import numpy as np
import pandas as pd
import pytest
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
    import pandas as pd
    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters.genfi_to_bids._utils import (
        _write_description_and_participants,
    )

    participants = pd.DataFrame(
        {
            "participant_id": ["sub-C9ORF002", "sub-MAPT003"],
            "session_id": [None, None],
            "modality": [None, None],
            "run_num": [None, None],
            "bids_filename": [None, None],
            "source": [None, None],
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
    import pandas as pd
    from fsspec.implementations.local import LocalFileSystem

    from clinica.converters.genfi_to_bids._utils import _write_sessions

    participant_id = "sub-C9ORF004"

    (tmp_path / participant_id / "ses-M000").mkdir(parents=True)
    (tmp_path / participant_id / "ses-M070").mkdir(parents=True)

    sessions = pd.DataFrame(
        {
            "participant_id": [participant_id, participant_id],
            "modality": [None, None],
            "bids_filename": [None, None],
            "run_num": [None, None],
            "session_id": ["ses-M000", "ses-M070"],
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
    session_ids = (
        set(df_out["session_id"])
        if "session_id" in df_out.columns
        else set(df_out.index.astype(str))
    )
    assert session_ids == {"ses-M000", "ses-M070"}


def test_write_scans(tmp_path, monkeypatch):
    import subprocess

    import pandas as pd

    from clinica.converters.genfi_to_bids._utils import _write_scans

    # Make dcm2niix succeed without running anything to bypass it
    class FakeCompleted:
        def __init__(self):
            self.returncode = 0
            self.stdout = b""
            self.stderr = b""

    monkeypatch.setattr(subprocess, "run", lambda *a, **k: FakeCompleted())

    participant_id = "sub-C9ORF006"
    session_id = "ses-M000"
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

    _write_scans(to=tmp_path, source=source, scans=scans)

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

    print(list(df_out.columns))

    assert list(df_out.columns) == expected_cols
    assert len(df_out) == 1

    row = df_out.iloc[0]

    assert row["modality"] == "T1"
    assert row["run_num"] == "run-01"
    assert row["filename"] == f"{participant_id}_{session_id}_run-01_T1w"
    assert not str(row["source_path"]).startswith(str(source))
    assert row["manufacturer"] == "SIEMENS"
    assert row["number_of_parts"] == 1
