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


def get_expected_path(tmp_path: Path, expected: Iterable[str]) -> set[Path]:
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
def test_get_images_with_suffix(tmp_path, suffix_directory_builder, suffixes, expected):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _get_images_with_suffix,
    )

    assert set(_get_images_with_suffix(tmp_path, suffixes)) == get_expected_path(
        tmp_path, expected
    )


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
    tmp_path, suffix_directory_builder, suffixes, expected
):
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        _remove_files_with_unsupported_suffixes,
    )

    assert set(
        _remove_files_with_unsupported_suffixes(tmp_path, suffixes)
    ) == get_expected_path(tmp_path, expected)


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


def test_adni_study_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import ADNIStudy  # noqa

    with pytest.raises(
        ValueError,
        match="'foo' is not a valid ADNIStudy",
    ):
        ADNIStudy("foo")


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
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id  # noqa

    with pytest.raises(
        ValueError,
        match=(
            f"DataFrame does not contain a column named '{expected_visit_code}', "
            "which is supposed to encode the visit code."
        ),
    ):
        _compute_session_id(pd.DataFrame(), csv_filename)


def test_compute_session_id_visit_code_wrong_format_error():
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id  # noqa

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
    from clinica.iotools.converters.adni_to_bids.adni_utils import _compute_session_id  # noqa

    input_data = {visit_code_column_name: ["f", "M00", "uns1", "sc", "M012", "M2368"]}
    expected_data = {
        **input_data,
        **{"session_id": [None, "ses-M000", None, "sc", "ses-M012", "ses-M2368"]},
    }

    assert_frame_equal(
        _compute_session_id(pd.DataFrame(input_data), csv_filename),
        pd.DataFrame(expected_data),
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
