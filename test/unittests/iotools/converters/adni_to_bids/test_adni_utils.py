import pandas as pd
import pytest


@pytest.mark.parametrize(
    "input_value,expected",
    [
        ("sub-ADNI000S0000", "000_S_0000"),
        ("sub-ADNI123S4567", "123_S_4567"),
        ("sub-ADNI12S4567", "12_S_4567"),
        ("sub-ADNI123X4567", "123_S_4567"),
        ("sub-ADNI123XYZ4567", "123_S_4567"),
        ("sub-ADNI123XYZ_TT4567", "123_S_4567"),
        ("sub-ADNI123XYZ12TT4567", None),
        ("", None),
        ("foo", None),
        ("12", None),
        ("123_S_4567", "123_S_4567"),
        ("1_XY_22", "1_S_22"),
    ],
)
def test_bids_id_to_loni(input_value, expected):
    """Test function `bids_id_to_loni`."""
    from clinica.iotools.converters.adni_to_bids.adni_utils import bids_id_to_loni

    assert bids_id_to_loni(input_value) == expected


def test_compute_session_id_visit_code_column_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id

    with pytest.raises(
        ValueError,
        match=(
            "DataFrame does not contain a column named 'VISCODE', "
            "which is supposed to encode the visit code."
        ),
    ):
        _compute_session_id(pd.DataFrame(), "foo.csv")


def test_compute_session_id_visit_code_wrong_format_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id

    df = pd.DataFrame({"VISCODE": ["foo", "bar", "baz", "bar", "foo", "foo"]})

    with pytest.raises(ValueError, match="The viscode foo is not correctly formatted."):
        _compute_session_id(df, "foo.csv")
