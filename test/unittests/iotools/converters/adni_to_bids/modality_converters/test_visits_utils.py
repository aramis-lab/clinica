import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from clinica.iotools.converters.adni_to_bids.adni_utils import ADNIStudy


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
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_preferred_visit_name,  # noqa
    )

    assert _get_preferred_visit_name(ADNIStudy(study_name), visit_code) == expected


def test_get_preferred_visit_name_visit_code_error():
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_preferred_visit_name,  # noqa
    )

    with pytest.raises(
        ValueError,
        match="Cannot extract month from visit code",
    ):
        _get_preferred_visit_name(ADNIStudy("ADNI3"), "")


@pytest.mark.parametrize(
    "visits,expected", [([], 0), ([pd.Series({"EXAMDATE": "2010-06-06"})], 1)]
)
def test_get_closest_and_second_closest_visits_errors(visits, expected):
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_closest_and_second_closest_visits,  # noqa
    )

    with pytest.raises(
        ValueError,
        match=f"Expected at least two visits, received {expected}.",
    ):
        _get_closest_and_second_closest_visits(visits, "2003-01-01")


@pytest.fixture
def visits() -> list[pd.Series]:
    return [
        pd.Series({"EXAMDATE": date})
        for date in (
            "2010-06-06",
            "2015-10-23",
            "2005-02-13",
            "2018-02-16",
            "2018-06-01",
        )
    ]


@pytest.mark.parametrize(
    (
        "image_acquisition_date,"
        "expected_closest_visit_date,"
        "expected_second_closest_visit_date,"
        "expected_difference_with_closest,"
        "expected_difference_with_second_closest,"
    ),
    [
        ("2003-01-01", "2005-02-13", "2010-06-06", 774, 2713),
        ("2008-01-01", "2010-06-06", "2005-02-13", 887, 1052),
        ("2015-10-23", "2015-10-23", "2018-02-16", 0, 847),
    ],
)
def test_get_closest_and_second_closest_visits(
    visits,
    image_acquisition_date,
    expected_closest_visit_date,
    expected_second_closest_visit_date,
    expected_difference_with_closest,
    expected_difference_with_second_closest,
):
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_closest_and_second_closest_visits,  # noqa
    )

    closest, second, diff1, diff2 = _get_closest_and_second_closest_visits(
        visits, image_acquisition_date
    )

    assert closest.EXAMDATE == expected_closest_visit_date
    assert second.EXAMDATE == expected_second_closest_visit_date
    assert diff1 == expected_difference_with_closest
    assert diff2 == expected_difference_with_second_closest


def test_get_smallest_time_difference_too_large_message():
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_smallest_time_difference_too_large_message,  # noqa
    )

    assert _get_smallest_time_difference_too_large_message(
        "sub-01",
        "baseline",
        "2020-01-01",
        pd.Series(
            {"EXAMDATE": "2021-01-01", "VISCODE": "ses-M12", "ORIGPROT": "ADNI1"}
        ),
        pd.Series(
            {"EXAMDATE": "2022-01-01", "VISCODE": "ses-M24", "ORIGPROT": "ADNI1"}
        ),
        365,
        365 * 2,
    ) == (
        "More than 90 days for corresponding timepoint in ADNIMERGE for subject sub-01 "
        "in visit baseline on 2020-01-01.\n"
        "Timepoint 1: ses-M12 - ADNI1 on 2021-01-01 (Distance: 365 days)\n"
        "Timepoint 2: ses-M24 - ADNI1 on 2022-01-01 (Distance: 730 days)\n"
    )


def test_get_closest_visit_empty():
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_closest_visit,
    )

    assert _get_closest_visit("2012-03-04", "control visit", [], "sub-01") is None


@pytest.fixture
def closest_visit_timepoints() -> list[pd.Series]:
    return [
        pd.Series(
            {
                "ORIGPROT": "ADNI1",
                "VISCODE": "bl",
                "COLPROT": "foo",
                "EXAMDATE": "2012-01-01",
            }
        ),
        pd.Series(
            {
                "ORIGPROT": "ADNI1",
                "VISCODE": "m03",
                "COLPROT": "bar",
                "EXAMDATE": "2012-03-01",
            }
        ),
        pd.Series(
            {
                "ORIGPROT": "ADNI2",
                "VISCODE": "m06",
                "COLPROT": "baz",
                "EXAMDATE": "2012-06-01",
            }
        ),
        pd.Series(
            {
                "ORIGPROT": "ADNI2",
                "VISCODE": "mfoo",
                "COLPROT": "foobarbaz",
                "EXAMDATE": "2012-08-11",
            }
        ),
        pd.Series(
            {
                "ORIGPROT": "ADNI2",
                "VISCODE": "m12",
                "COLPROT": "foobar",
                "EXAMDATE": "2013-03-20",
            }
        ),
    ]


@pytest.mark.parametrize(
    "image_acquisition_date", ["1809-03-04", "2012-01-01", "2066-12-31"]
)
def test_get_closest_visit_single_visit(
    closest_visit_timepoints, image_acquisition_date
):
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_closest_visit,
    )

    assert_series_equal(
        _get_closest_visit(
            image_acquisition_date,
            "control visit",
            closest_visit_timepoints[0:1],  # fake a single visit list
            "sub-01",
        ),
        closest_visit_timepoints[0],
    )


@pytest.mark.parametrize(
    "image_acquisition_date,expected",
    [
        ("2012-03-04", 1),
        ("1989-06-16", 0),
        ("2020-11-26", -1),
        (
            "2012-12-12",
            -2,
        ),  # Special case where the second-closest visit is preferred over the closest one
    ],
)
def test_get_closest_visit(closest_visit_timepoints, image_acquisition_date, expected):
    from clinica.iotools.converters.adni_to_bids.modality_converters._visits_utils import (
        _get_closest_visit,
    )

    assert_series_equal(
        _get_closest_visit(
            image_acquisition_date,
            "control visit",
            closest_visit_timepoints,
            "sub-01",
        ),
        closest_visit_timepoints[expected],
    )
