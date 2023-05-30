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


def test_adni_study_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import ADNIStudy

    with pytest.raises(
        ValueError,
        match="Invalid study name for ADNI: foo.",
    ):
        ADNIStudy.from_string("foo")


@pytest.mark.parametrize(
    "study_name,visit_code,expected",
    [
        ("ADNI3", "bl", "ADNI Screening"),
        ("ADNI3", "M000", "ADNI3 Year 0.0 Visit"),
        ("ADNI3", "M006", "ADNI3 Year 0.5 Visit"),
        ("ADNI3", "M012", "ADNI3 Year 1.0 Visit"),
        ("ADNI3", "M013", "ADNI3 Year 1.0833333333333333 Visit"),
        ("ADNI2", "bl", "ADNI2 Screening MRI-New Pt"),
        ("ADNI2", "m03", "ADNI2 Month 3 MRI-New Pt"),
        ("ADNI2", "m06", "ADNI2 Month 6-New Pt"),
        ("ADNI2", "M000", "ADNI2 Year 0.0 Visit"),
        ("ADNI2", "M006", "ADNI2 Year 0.5 Visit"),
        ("ADNI2", "M012", "ADNI2 Year 1.0 Visit"),
        ("ADNI2", "M013", "ADNI2 Year 1.0833333333333333 Visit"),
        ("ADNI1", "bl", "ADNI Screening"),
        ("ADNIGO", "bl", "ADNIGO Screening MRI"),
        ("ADNI1", "m03", "ADNIGO Month 3 MRI"),
        ("ADNIGO", "m03", "ADNIGO Month 3 MRI"),
        ("ADNI1", "m00", "ADNI1/GO Month 0"),
        ("ADNIGO", "m00", "ADNI1/GO Month 0"),
        ("ADNI1", "m07", "ADNI1/GO Month 7"),
        ("ADNIGO", "m07", "ADNI1/GO Month 7"),
        ("ADNI1", "m12", "ADNI1/GO Month 12"),
        ("ADNIGO", "m12", "ADNI1/GO Month 12"),
        ("ADNI1", "m53", "ADNI1/GO Month 53"),
        ("ADNIGO", "m53", "ADNI1/GO Month 53"),
        ("ADNI1", "m54", "ADNIGO Month 54"),
        ("ADNIGO", "m54", "ADNIGO Month 54"),
    ],
)
def test_get_preferred_visit_name(study_name, visit_code, expected):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        ADNIStudy,
        _get_preferred_visit_name,
    )

    assert (
        _get_preferred_visit_name(ADNIStudy.from_string(study_name), visit_code)
        == expected
    )


def test_get_preferred_visit_name_visit_code_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        ADNIStudy,
        _get_preferred_visit_name,
    )

    with pytest.raises(
        ValueError,
        match="Cannot extract month from visit code",
    ):
        _get_preferred_visit_name(ADNIStudy.from_string("ADNI3"), "")


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
