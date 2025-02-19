import pytest


@pytest.mark.parametrize(
    "subjects,sessions,expected_message",
    [
        ((), (), "The pipeline will be run on 0 image(s)."),
        (
            ("sub-01",),
            ("ses-M000",),
            "The pipeline will be run on the following 1 image(s):\n\t- sub-01 | ses-M000",
        ),
        (
            ("sub-01", "sub-01"),
            ("ses-M000", "ses-M006"),
            "The pipeline will be run on the following 2 image(s):\n\t- sub-01 | ses-M000, ses-M006",
        ),
        (
            ("sub-01", "sub-01", "sub-02"),
            ("ses-M000", "ses-M006", "ses-M000"),
            "The pipeline will be run on the following 3 image(s):\n\t- sub-01 | ses-M000, ses-M006\n\t- sub-02 | ses-M000",
        ),
        (
            ("sub-01", "sub-01", "sub-02", "sub-03", "sub-03", "sub-03"),
            ("ses-M000", "ses-M006", "ses-M000", "ses-M006", "ses-M012", "ses-M024"),
            "The pipeline will be run on the following 6 image(s):\n\t- sub-01 | ses-M000, ses-M006\n\t- sub-02 | ses-M000\n\t- sub-03 | ses-M006, ses-M012, ses-M024",
        ),
        (
            ("sub-01", "sub-01", "sub-02", "sub-03", "sub-04", "sub-05", "sub-06"),
            (
                "ses-M000",
                "ses-M006",
                "ses-M000",
                "ses-M000",
                "ses-M000",
                "ses-M000",
                "ses-M000",
            ),
            "The pipeline will be run on the following 7 image(s):\n\t- sub-01 | ses-M000, ses-M006\n\t- sub-02 | ses-M000\n\t- sub-03 | ses-M000\n\t- sub-04 | ses-M000\n\t- sub-05 | ses-M000\n\t- sub-06 | ses-M000",
        ),
    ],
)
def test_get_message_images_to_process(subjects, sessions, expected_message):
    from clinica.utils.ux import _get_message_images_to_process

    assert _get_message_images_to_process(subjects, sessions) == expected_message


@pytest.mark.parametrize("max_number_of_lines", [0, 1, 2])
def test_get_message_images_to_process_custom_max_number_of_lines_too_small(
    max_number_of_lines,
):  # subjects, sessions, expected_message):
    from clinica.utils.ux import _get_message_images_to_process

    assert (
        _get_message_images_to_process(
            ("sub-01", "sub-01", "sub-02"),
            ("ses-M000", "ses-M006", "ses-M000"),
            max_number_of_lines=max_number_of_lines,
        )
        == "The pipeline will be run on 3 image(s)."
    )


@pytest.mark.parametrize(
    "subjects,sessions,max_number_of_lines,expected_message",
    [
        (
            ("sub-01", "sub-01", "sub-02"),
            ("ses-M000", "ses-M006", "ses-M000"),
            4,
            (
                "The pipeline will be run on the following 3 image(s):"
                "\n\t- sub-01 | ses-M000, ses-M006\n\t- sub-02 | ses-M000"
            ),
        ),
        (
            ("sub-01", "sub-01", "sub-02", "sub-03", "sub-04", "sub-05", "sub-06"),
            (
                "ses-M000",
                "ses-M006",
                "ses-M000",
                "ses-M000",
                "ses-M000",
                "ses-M000",
                "ses-M000",
            ),
            4,
            (
                "The pipeline will be run on the following 7 image(s):"
                "\n\t- sub-01 | ses-M000, ses-M006\n\t- sub-02 | ses-M000"
                "\n\t\t...\n\t- sub-06 | ses-M000"
            ),
        ),
    ],
)
def test_get_message_images_to_process_custom_max_number_of_lines(
    subjects, sessions, max_number_of_lines, expected_message
):
    from clinica.utils.ux import _get_message_images_to_process

    assert (
        _get_message_images_to_process(
            subjects, sessions, max_number_of_lines=max_number_of_lines
        )
        == expected_message
    )


def test_get_group_message_empty(tmp_path):
    from clinica.utils.ux import _get_group_message

    assert _get_group_message(tmp_path) == "No group was found in CAPS directory"

    (tmp_path / "foo.nii.gz").touch()

    assert _get_group_message(tmp_path) == "No group was found in CAPS directory"


def test_get_group_message(tmp_path):
    from clinica.utils.ux import _get_group_message

    (tmp_path / "groups").mkdir()
    for group in ("group-AD", "group-foo", "group-foobar66"):
        (tmp_path / "groups" / group).mkdir()

    assert (
        _get_group_message(tmp_path)
        == "Groups that exist in your CAPS directory are AD, foo, foobar66."
    )
