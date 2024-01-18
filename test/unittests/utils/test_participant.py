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
            (["sub-01", "sub-03", "sub-03"], ["ses-M000", "ses-M006", "ses-M000"]),
        ),
    ],
)
def test_get_subject_session_list(tmp_path, config, expected):
    from clinica.utils.participant import get_subject_session_list
    from clinica.utils.testing_utils import build_bids_directory

    build_bids_directory(tmp_path, config)

    assert get_subject_session_list(tmp_path) == expected
