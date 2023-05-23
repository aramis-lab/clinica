import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


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


@pytest.mark.parametrize(
    "csv_filename,expected_visit_code",
    [
        ("foo.csv", "VISCODE"),
        ("MOCA.csv", "VISCODE2"),
        ("UWNPSYCHSUM_03_07_19.csv", "VISCODE2"),
        ("BHR_EVERYDAY_COGNITION.csv", "Timepoint"),
        ("BHR_BASELINE_QUESTIONNAIRE.csv", "Timepoint"),
        ("BHR_LONGITUDINAL_QUESTIONNAIRE.csv", "Timepoint"),
    ],
)
def test_compute_session_id_visit_code_column_error(csv_filename, expected_visit_code):
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id

    with pytest.raises(
        ValueError,
        match=(
            f"DataFrame does not contain a column named '{expected_visit_code}', "
            "which is supposed to encode the visit code."
        ),
    ):
        _compute_session_id(pd.DataFrame(), csv_filename)


def test_compute_session_id_visit_code_wrong_format_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id

    df = pd.DataFrame({"VISCODE": ["foo", "bar", "baz", "bar", "foo", "foo"]})

    with pytest.raises(ValueError, match="The viscode foo is not correctly formatted."):
        _compute_session_id(df, "foo.csv")


@pytest.mark.parametrize(
    "csv_filename,visit_code_column_name",
    [
        ("foo.csv", "VISCODE"),
        ("MOCA.csv", "VISCODE2"),
        ("UWNPSYCHSUM_03_07_19.csv", "VISCODE2"),
        ("BHR_EVERYDAY_COGNITION.csv", "Timepoint"),
        ("BHR_BASELINE_QUESTIONNAIRE.csv", "Timepoint"),
        ("BHR_LONGITUDINAL_QUESTIONNAIRE.csv", "Timepoint"),
    ],
)
def test_compute_session_id(csv_filename, visit_code_column_name):
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id

    input_data = {visit_code_column_name: ["f", "M00", "uns1", "sc", "M012", "M2368"]}
    expected_data = {
        **input_data,
        **{"session_id": [None, "ses-M000", None, "sc", "ses-M012", "ses-M2368"]},
    }

    assert_frame_equal(
        _compute_session_id(pd.DataFrame(input_data), csv_filename),
        pd.DataFrame(expected_data),
    )
