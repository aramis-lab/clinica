import pandas as pd
import pytest


@pytest.mark.parametrize(
    "subjects,sessions,expected",
    [
        ([], [], ([], [])),
        (
            ["sub-CLNC01", "sub-CLNC01", "sub-CLNC02"],
            ["ses-M000", "ses-M018", "ses-M000"],
            (
                ["sub-CLNC01", "sub-CLNC02"],
                [["ses-M000", "ses-M018"], ["ses-M000"]],
            ),
        ),
    ],
)
def test_get_unique_subjects(subjects, sessions, expected):
    from clinica.utils.participant import get_unique_subjects

    assert get_unique_subjects(subjects, sessions) == expected


def test_get_unique_subjects_error():
    from clinica.utils.participant import get_unique_subjects

    with pytest.raises(
        ValueError,
        match="The number of subjects should match the number of sessions.",
    ):
        get_unique_subjects(["foo"], [])


@pytest.mark.parametrize(
    "subjects,session_per_subject,expected",
    [
        ([], [], ([], [])),
        (
            ["sub-01", "sub-02"],
            [["ses-M000", "ses-M018"], ["ses-M000"]],
            (
                ["sub-01", "sub-01", "sub-02"],
                ["ses-M000", "ses-M018", "ses-M000"],
            ),
        ),
    ],
)
def test_unique_subjects_sessions_to_subjects_sessions(
    subjects, session_per_subject, expected
):
    from clinica.utils.participant import unique_subjects_sessions_to_subjects_sessions

    assert (
        unique_subjects_sessions_to_subjects_sessions(subjects, session_per_subject)
        == expected
    )


def test_unique_subjects_sessions_to_subjects_sessions_error():
    from clinica.utils.participant import unique_subjects_sessions_to_subjects_sessions

    with pytest.raises(
        ValueError,
        match="The number of unique subjects should match the number of session lists.",
    ):
        unique_subjects_sessions_to_subjects_sessions(["foo"], [])


def test_get_subject_session_list_error(tmp_path):
    from clinica.utils.participant import get_subject_session_list

    with pytest.raises(
        OSError,
        match="Dataset empty or not BIDS/CAPS compliant.",
    ):
        get_subject_session_list(tmp_path)


@pytest.mark.parametrize(
    "config,expected",
    [
        ({"sub-01": []}, ([], [])),
        ({"sub-01": ["ses-M000"]}, (["sub-01"], ["ses-M000"])),
        ({"sub-01": ["ses-M000"], "sub-03": []}, (["sub-01"], ["ses-M000"])),
        (
            {"sub-01": ["ses-M000"], "sub-03": ["ses-M000", "ses-M006"]},
            (["sub-01", "sub-03", "sub-03"], ["ses-M000", "ses-M000", "ses-M006"]),
        ),
    ],
)
def test_get_subject_session_list(tmp_path, config, expected):
    from clinica.utils.participant import get_subject_session_list
    from clinica.utils.testing_utils import build_bids_directory

    (tmp_path / "bids").mkdir()
    build_bids_directory(tmp_path / "bids", config)
    (tmp_path / "tsv").mkdir()
    assert (
        get_subject_session_list(tmp_path / "bids", tsv_dir=tmp_path / "tsv")
        == expected
    )


def test_read_participant_tsv_error(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.participant import _read_participant_tsv

    with pytest.raises(
        ClinicaException,
        match="The TSV file you gave is not a file.",
    ):
        _read_participant_tsv(tmp_path / "foo.tsv")


def test_read_participant_tsv(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.participant import _read_participant_tsv

    df = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M000", "ses-M006", "ses-M000"],
        }
    )
    df.to_csv(tmp_path / "foo.tsv", sep="\t")

    assert _read_participant_tsv(tmp_path / "foo.tsv") == (
        ["sub-01", "sub-01", "sub-02"],
        ["ses-M000", "ses-M006", "ses-M000"],
    )

    for column in ("participant_id", "session_id"):
        df.drop(column, axis=1).to_csv(tmp_path / "foo.tsv", sep="\t")
        with pytest.raises(
            ClinicaException,
            match=f"The TSV file does not contain {column} column",
        ):
            _read_participant_tsv(tmp_path / "foo.tsv")
