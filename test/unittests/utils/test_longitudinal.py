import pytest


@pytest.mark.parametrize(
    "session_ids,expected",
    [
        (["ses-M000"], "long-M000"),
        (["ses-M000", "ses-M018", "ses-M036"], "long-M000+M018+M036"),
        (["ses-M018", "ses-M036", "ses-M000"], "long-M000+M018+M036"),
        (["ses-foo", "ses-bar", "ses-baz", "ses-foobar"], "long-bar+baz+foo+foobar"),
    ],
)
def test_get_long_id(session_ids, expected):
    from clinica.utils.longitudinal import get_long_id

    assert get_long_id(session_ids) == expected


@pytest.mark.parametrize(
    "session_ids",
    [
        ["foo", "bar", "baz"],
        ["1", "2", "3"],
        ["SES-M000", "SES-M001"],
    ],
)
def test_get_long_id_errors(session_ids):
    from clinica.utils.longitudinal import get_long_id

    with pytest.raises(
        ValueError,
        match="Expected a list of session IDs of the form ses-XXX",
    ):
        get_long_id(session_ids)


@pytest.mark.parametrize(
    "participant_ids,session_ids,expected",
    [
        (
            ["sub-CLNC01", "sub-CLNC01", "sub-CLNC02"],
            ["ses-M000", "ses-M018", "ses-M000"],
            ["long-M000+M018", "long-M000+M018", "long-M000"],
        )
    ],
)
def test_get_participants_long_id(participant_ids, session_ids, expected):
    from clinica.utils.longitudinal import get_participants_long_id

    assert get_participants_long_id(participant_ids, session_ids) == expected


def test_save_long_id(tmp_path):
    from clinica.utils.longitudinal import save_long_id

    output_dir = tmp_path / "out"
    save_long_id(
        ["ses-M000", "ses-M018", "ses-M036"], output_dir, "longitudinal_ids.tsv"
    )

    assert (output_dir / "longitudinal_ids.tsv").exists()
    saved_ids = (output_dir / "longitudinal_ids.tsv").read_text()
    assert saved_ids == "session_id\nses-M000\nses-M018\nses-M036\n"

    save_long_id(["ses-M000", "ses-M018", "ses-M036"], output_dir)

    assert (output_dir / "long-M000+M018+M036_sessions.tsv").exists()
    saved_ids = (output_dir / "long-M000+M018+M036_sessions.tsv").read_text()
    assert saved_ids == "session_id\nses-M000\nses-M018\nses-M036\n"


def test_read_sessions_no_file_error(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.longitudinal import read_sessions

    with pytest.raises(
        ClinicaException,
        match="The TSV file with sessions associated",
    ):
        read_sessions(tmp_path, "foo", "bar")


def test_read_sessions_wrong_column_error(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.longitudinal import read_sessions

    input_folder = tmp_path / "subjects" / "sub-01" / "long-M000M018M036"
    input_folder.mkdir(parents=True)
    tsv_file = input_folder / "long-M000M018M036_sessions.tsv"
    tsv_file.write_text("foo_bar\nses-M000\nses-M018\nses-M036\n")

    with pytest.raises(
        ClinicaException,
        match="The TSV file does not contain session_id column ",
    ):
        read_sessions(tmp_path, "sub-01", "long-M000M018M036")


def test_read_sessions(tmp_path):
    from clinica.utils.longitudinal import read_sessions

    input_folder = tmp_path / "subjects" / "sub-01" / "long-M000M018M036"
    input_folder.mkdir(parents=True)
    tsv_file = input_folder / "long-M000M018M036_sessions.tsv"
    tsv_file.write_text("session_id\nses-M000\nses-M018\nses-M036\n")

    assert read_sessions(tmp_path, "sub-01", "long-M000M018M036") == [
        "ses-M000",
        "ses-M018",
        "ses-M036",
    ]
