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


@pytest.mark.parametrize(
    "input,expected",
    [
        ("bl", "ses-M000"),
        ("m0", "ses-M000"),
        ("m3", "ses-M003"),
        ("m03", "ses-M003"),
        ("m003", "ses-M003"),
        ("m0003", "ses-M003"),
        ("m00", "ses-M000"),
        ("m0000000", "ses-M000"),
    ],
)
def test_viscode_to_session(input, expected):
    """Test function `viscode_to_session`."""

    from clinica.iotools.converter_utils import viscode_to_session

    assert viscode_to_session(input) == expected
