import pytest


@pytest.mark.parametrize(
    "filename",
    [
        "foo",
        "sub-01.json",
        "ses-M000.nii.gz",
        "sub-01_ses-M000.json",
        "sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
    ],
)
def test_get_subject_id_error(filename):
    from clinica.utils.filemanip import get_subject_id

    with pytest.raises(
        ValueError,
        match=(
            f"Input filename {filename} is not in a BIDS or CAPS compliant format. "
            "It does not contain the subject and session information."
        ),
    ):
        get_subject_id(filename)


@pytest.mark.parametrize(
    "filename,expected",
    [
        (
            "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
            "sub-01_ses-M000",
        ),
        (
            "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
            "sub-foo_ses-bar",
        ),
    ],
)
def test_get_subject_id(filename, expected):
    from clinica.utils.filemanip import get_subject_id

    assert get_subject_id(filename) == expected


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("foo.nii.gz", "foo"),
        ("sub-01/ses-M000/sub-01_ses-M000.tar.gz", "sub-01_ses-M000"),
        ("foo/bar/baz/foo-bar_baz.niml.dset", "foo-bar_baz"),
    ],
)
def test_get_filename_no_ext(filename, expected):
    from clinica.utils.filemanip import get_filename_no_ext

    assert get_filename_no_ext(filename) == expected


def test_extract_image_ids_error():
    from clinica.utils.filemanip import extract_image_ids

    with pytest.raises(
        ValueError,
        match=(
            "Input filename foo.bar is not in a BIDS or CAPS compliant format. "
            "It does not contain the subject and session information."
        ),
    ):
        extract_image_ids(["foo.bar"])


def test_extract_image_ids():
    from clinica.utils.filemanip import extract_image_ids

    assert (
        extract_image_ids(
            [
                "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
                "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
                "sub-01_ses-M000.tar.gz",
            ]
        )
        == ["sub-01_ses-M000"] * 3
    )


def test_extract_subjects_sessions_from_filename():
    from clinica.utils.filemanip import extract_subjects_sessions_from_filename

    assert (
        extract_subjects_sessions_from_filename(
            [
                "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
                "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
                "sub-01_ses-M000.tar.gz",
            ]
        )
    ) == (["sub-01", "sub-01", "sub-01"], ["ses-M000", "ses-M000", "ses-M000"])


def test_extract_crash_files_from_log_file_error():
    from clinica.utils.filemanip import extract_crash_files_from_log_file

    with pytest.raises(
        ValueError,
        match="extract_crash_files_from_log_file",
    ):
        extract_crash_files_from_log_file("foo.log")
