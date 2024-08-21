from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    _define_magnetic_field,
    _get_bids_filename_from_image_data,
    _get_img_data,
    _get_subjects_list_from_data,
    _get_subjects_list_from_file,
    _rename_clinical_data,
    _rename_modalities,
    read_clinical_data,
)


def test_get_subjects_list_from_data(tmp_path):
    for filename in ["IXI1", "IXI123", "IXIaaa", "foo"]:
        Path(f"{tmp_path}/{filename}_T1w.nii.gz").touch()
    assert _get_subjects_list_from_data(tmp_path) == ["IXI123"]


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("T1", "T1w"),
        ("T2", "T2w"),
        ("MRA", "angio"),
        ("PD", "PDw"),
        ("DTI", "dti"),
    ],
)
def test_rename_ixi_modalities_success(input_str, expected):
    assert _rename_modalities(input_str) == expected


@pytest.mark.parametrize("input_str", ["t1", "foo", "T1w"])
def test_rename_ixi_modalities_error(input_str):
    with pytest.raises(
        ValueError,
        match=f"The modality {input_str} is not recognized in the IXI dataset.",
    ):
        _rename_modalities(input_str)


def test_get_subjects_list_from_file(tmp_path):
    with open(tmp_path / "subjects.txt", "w") as f:
        f.write("IXI123\nIXI001")
    assert _get_subjects_list_from_file(Path(tmp_path / "subjects.txt")) == [
        "IXI123",
        "IXI001",
    ]


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("SEX_ID (1=m, 2=f)", "sex"),
        ("ETHNIC_ID", "ethnicity"),
        ("MARITAL_ID", "marital status"),
        ("OCCUPATION_ID", "occupation"),
        ("QUALIFICATION_ID", "qualification"),
        ("IXI_ID", "source_id"),
        ("DOB", "date of birth"),
        ("STUDY_DATE", "acq_time"),
        ("FOO", "foo"),
    ],
)
def test_rename_clinical_data(input_str, expected):
    assert _rename_clinical_data(input_str) == expected


def test_read_clinical_data(tmp_path):
    # todo : need to mock the _get functions ?
    pass


@pytest.mark.parametrize(
    "input_str, expected",
    [
        ("Guys", "1.5"),
        ("IOP", "1.5"),
        ("HH", "3"),
    ],
)
def test_define_magnetic_field_success(input_str, expected):
    assert _define_magnetic_field(input_str) == expected


def test_define_magnetic_field_success():
    with pytest.raises(
        ValueError,
        match=f"The hospital {'foo'} was not recognized.",
    ):
        _define_magnetic_field("foo")


def image_dataframe_builder(tmp_path: Path) -> pd.DataFrame:
    tmp_image = tmp_path / "IXI001-Guys-1234-T1.nii.gz"
    tmp_image.touch()
    (tmp_path / "IXI001-Guys-1234-T1").touch()
    (tmp_path / "IXI001-Guys-T1.nii.gz").touch()
    (tmp_path / "IXI001-Guys-1234-T1-00.nii.gz").touch()

    return pd.DataFrame(
        {
            "img_path": [tmp_image],
            "img_name": ["IXI001-Guys-1234-T1.nii.gz"],
            "img_name_no_ext": ["IXI001-Guys-1234-T1"],
            "subject": ["IXI001"],
            "participant_id": ["sub-IXI001"],
            "hospital": ["Guys"],
            "modality": ["T1w"],
            "field": ["1.5"],
            "session": ["ses-M000"],
        }
    )


def test_get_image_data(tmp_path):
    input = image_dataframe_builder(tmp_path)
    assert _get_img_data(tmp_path).equals(input)


def test_get_bids_filename_from_image_data(tmp_path):
    input = image_dataframe_builder(tmp_path)
    assert _get_bids_filename_from_image_data(input.loc[0]) == "sub-IXI001_ses-M000_T1w"
