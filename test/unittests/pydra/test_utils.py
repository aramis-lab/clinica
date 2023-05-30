import pytest


@pytest.mark.parametrize(
    "input_fwhm,expected",
    [
        (3.0, [[3.0, 3.0, 3.0]]),
        ([3.0], [[3.0, 3.0, 3.0]]),
        ((3,), [[3.0, 3.0, 3.0]]),
        ([3.0, 2.0], [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0]]),
        ((3.0, 2.0), [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0]]),
        ([3.0, 2.0, 1.0], [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0], [1.0, 1.0, 1.0]]),
        ([[3.0, 2.0, 1.0], [2.0, 2.0, 1.0]], [[3.0, 2.0, 1.0], [2.0, 2.0, 1.0]]),
    ],
)
def test_sanitize_fwhm(input_fwhm, expected):
    from clinica.pydra.utils import sanitize_fwhm

    assert sanitize_fwhm(input_fwhm) == expected
