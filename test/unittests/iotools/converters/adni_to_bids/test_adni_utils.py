from pathlib import Path
from typing import Iterable, List

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal


@pytest.fixture
def suffix_directory_builder(
    tmp_path: Path,
):
    filenames = {
        "sub-01_ses-M001_magnitude1.nii.gz",
        "sub-01_ses-M001_fmap_real.nii.gz",
        "sub-01_ses-M001_fmap_imaginary.nii",
        "sub-01_ses-M001_T1w_ADC.nii.gz",
        "sub-01_ses-M001_T1w.nii.gz",
        "sub-01_ses-M001_T1wADC.nii.gz",
        "sub-01_ses-M001_FLAIR_ADC.nii.gz",
        "sub-01_ses-M001_pet_ADC.json",
    }

    for file in filenames:
        (tmp_path / file).touch()


@pytest.fixture
def get_expected_path(tmp_path: Path, expected: Iterable[str]):
    return set([tmp_path / name for name in expected])


@pytest.mark.parametrize(
    "suffixes, expected",
    [
        ((), {}),
        (
            ("FLAIR",),
            {
                "sub-01_ses-M001_FLAIR_ADC.nii.gz",
            },
        ),
        (
            ("real", "imaginary"),
            {
                "sub-01_ses-M001_fmap_real.nii.gz",
                "sub-01_ses-M001_fmap_imaginary.nii",
            },
        ),
    ],
)
def test_get_images_with_suffix(
    tmp_path, suffix_directory_builder, suffixes, get_expected_path, expected
):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _get_images_with_suffix,
    )

    suffix_directory_builder
    assert set(_get_images_with_suffix(tmp_path, suffixes)) == get_expected_path


@pytest.mark.parametrize(
    "update, suffixes, expected",
    [
        (True, ("foo",), True),
        (True, ("T1w",), True),
        (False, ("imaginary", "T1w"), False),
        (False, (), True),
    ],
)
def test_remove_existing_images_if_necessary(
    tmp_path, suffix_directory_builder, suffixes, update, expected
):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _remove_existing_images_if_necessary,
    )

    suffix_directory_builder
    assert _remove_existing_images_if_necessary(tmp_path, suffixes, update) == expected


@pytest.mark.parametrize(
    "suffixes, expected",
    [
        (
            ("ADC", "real", "imaginary"),
            {
                "sub-01_ses-M001_fmap_real.nii.gz",
                "sub-01_ses-M001_fmap_imaginary.nii",
                "sub-01_ses-M001_T1w_ADC.nii.gz",
                "sub-01_ses-M001_T1wADC.nii.gz",
                "sub-01_ses-M001_T1wADC.nii.gz",
                "sub-01_ses-M001_FLAIR_ADC.nii.gz",
                "sub-01_ses-M001_pet_ADC.json",
            },
        ),
        (
            ("foo",),
            {},
        ),
    ],
)
def test_remove_files_with_unsupported_suffixes(
    tmp_path, suffix_directory_builder, get_expected_path, suffixes, expected
):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _remove_files_with_unsupported_suffixes,
    )

    suffix_directory_builder
    assert (
        set(_remove_files_with_unsupported_suffixes(tmp_path, suffixes))
        == get_expected_path
    )


@pytest.mark.parametrize(
    "input, expected",
    [
        (
            {"001_S_0001", "001_S_0002", "001_S_0003"},
            {"001_S_0001", "001_S_0002", "001_S_0003"},
        ),
        ({"001_S_0001", "001_S_00014", ".001_S_0001", "001S0001"}, {"001_S_0001"}),
    ],
)
def test_define_subjects_list_directory(tmp_path, input, expected):
    from clinica.iotools.converters.adni_to_bids.adni_utils import _define_subjects_list

    source_dir = tmp_path / "source_dir"
    source_dir.mkdir()

    for subject in input:
        (source_dir / subject).touch()

    assert set(_define_subjects_list(source_dir)) == expected


def test_define_subjects_list_txt(tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_utils import _define_subjects_list

    source_dir = tmp_path / "source_dir"
    subjs_list_path = tmp_path / "subjects_list.txt"
    input = {"001_S_0001", "001_S_00022", "001S0003"}
    with open(subjs_list_path, "w") as f:
        f.write("\n".join(input))

    assert set(_define_subjects_list(source_dir, subjs_list_path)) == input


@pytest.mark.parametrize(
    "write_all, input, expected",
    [
        (True, ["001_S_0001", "001_S_0002"], {"001_S_0001", "001_S_0002"}),
        (False, ["001_S_0001", "001_S_0002"], {"001_S_0001"}),
    ],
)
def test_check_subjects_list(tmp_path, write_all, input, expected):
    from clinica.iotools.converters.adni_to_bids.adni_utils import _check_subjects_list

    clinical_dir = tmp_path / "clinical_dir"
    clinical_dir.mkdir()

    if not write_all:
        input.pop()

    adni_df = pd.DataFrame(columns=["PTID"], data=input)
    adni_df.to_csv(clinical_dir / "ADNIMERGE.csv")

    assert set(_check_subjects_list(input, clinical_dir)) == expected


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


@pytest.mark.parametrize(
    "visits,expected", [([], 0), ([pd.Series({"EXAMDATE": "2010-06-06"})], 1)]
)
def test_get_closest_and_second_closest_visits_errors(visits, expected):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _get_closest_and_second_closest_visits,
    )

    with pytest.raises(
        ValueError,
        match=f"Expected at least two visits, received {expected}.",
    ):
        _get_closest_and_second_closest_visits(visits, "2003-01-01")


@pytest.fixture
def visits() -> List[pd.Series]:
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
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _get_closest_and_second_closest_visits,
    )

    closest, second, diff1, diff2 = _get_closest_and_second_closest_visits(
        visits, image_acquisition_date
    )

    assert closest.EXAMDATE == expected_closest_visit_date
    assert second.EXAMDATE == expected_second_closest_visit_date
    assert diff1 == expected_difference_with_closest
    assert diff2 == expected_difference_with_second_closest


def test_get_smallest_time_difference_too_large_message():
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _get_smallest_time_difference_too_large_message,
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
    from clinica.iotools.converters.adni_to_bids.adni_utils import _get_closest_visit

    assert _get_closest_visit("2012-03-04", "control visit", [], "sub-01") is None


@pytest.fixture
def closest_visit_timepoints() -> List[pd.Series]:
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
    from clinica.iotools.converters.adni_to_bids.adni_utils import _get_closest_visit

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
    from clinica.iotools.converters.adni_to_bids.adni_utils import _get_closest_visit

    assert_series_equal(
        _get_closest_visit(
            image_acquisition_date,
            "control visit",
            closest_visit_timepoints,
            "sub-01",
        ),
        closest_visit_timepoints[expected],
    )


@pytest.fixture
def input_df():
    return pd.DataFrame(
        {
            "foo": ["foo1", "foo2", "foo3"],
            "bar": [1, 2, 3],
            "baz": [True, False, False],
            "foobar": [4, 5, 6],
        }
    )


@pytest.mark.parametrize(
    "csv_name,csv_to_look_for",
    [("adnimerge.csv", "adnimerge"), ("adnimerge_20Oct2023.csv", "adnimerge")],
)
def test_load_clinical_csv(tmp_path, input_df, csv_name, csv_to_look_for):
    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    input_df.to_csv(tmp_path / csv_name, index=False)
    assert_frame_equal(load_clinical_csv(tmp_path, csv_to_look_for), input_df)


def test_load_clinical_csv_error(
    tmp_path,
):
    import re

    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    pattern = r"(_\d{1,2}[A-Za-z]{3}\d{4})?.csv"
    with pytest.raises(
        IOError,
        match=re.escape(
            f"Expecting to find exactly one file in folder {tmp_path} "
            f"matching pattern adnimerge{pattern}. 0 "
            f"files were found instead : \n[- ]"
        ),
    ):
        load_clinical_csv(tmp_path, "adnimerge")


def test_load_clinical_csv_value_error(tmp_path):
    import re

    from clinica.iotools.converters.adni_to_bids.adni_utils import load_clinical_csv

    with open(tmp_path / "adnimerge.csv", "w") as fp:
        fp.write("col1,col2,col3\n1,2,3\n1,2,3,4")

    with pytest.raises(
        ValueError,
        match=f"File {tmp_path}/adnimerge.csv was found but could not "
        "be loaded as a DataFrame. Please check your data.",
    ):
        load_clinical_csv(tmp_path, "adnimerge")
