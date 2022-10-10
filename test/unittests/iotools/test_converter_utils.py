import pytest


@pytest.mark.parametrize(
    "input_list,expected",
    [
        (
            ["ses-M000", "ses-M006", "ses-M012", "ses-M024", "ses-M048", "ses-M003"],
            ["ses-M000", "ses-M003", "ses-M006", "ses-M012", "ses-M024", "ses-M048"],
        ),
        (
            ["ses-M00", "ses-M06", "ses-M12", "ses-M24", "ses-M48", "ses-M03"],
            ["ses-M00", "ses-M03", "ses-M06", "ses-M12", "ses-M24", "ses-M48"],
        ),
        (
            ["ses-M0", "ses-M6", "ses-M12", "ses-M24", "ses-M48", "ses-M3"],
            ["ses-M0", "ses-M3", "ses-M6", "ses-M12", "ses-M24", "ses-M48"],
        ),
    ],
)
def test_sort_session_list(input_list, expected):
    """Test function `sort_session_list`."""
    from clinica.iotools.converter_utils import sort_session_list

    assert sort_session_list(input_list) == expected


def test_viscode_to_session():
    """Test function `viscode_to_session`."""

    from clinica.iotools.converter_utils import viscode_to_session

    assert viscode_to_session("bl") == "ses-M000"
    assert viscode_to_session("m0") == "ses-M000"
    assert viscode_to_session("m3") == "ses-M003"
    assert viscode_to_session("m03") == "ses-M003"
    assert viscode_to_session("m003") == "ses-M003"
    assert viscode_to_session("m0003") == "ses-M003"
    assert viscode_to_session("m00") == "ses-M000"
