from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    _define_magnetic_field,
    _get_bids_filename_from_image_data,
    _get_ethnic_mapping,
    _get_img_data,
    _get_marital_mapping,
    _get_occupation_mapping,
    _get_qualification_mapping,
    _get_subjects_list_from_data,
    _get_subjects_list_from_file,
    _rename_clinical_data,
    _rename_modalities,
    read_clinical_data,
)


def test_get_subjects_list_from_data(tmp_path):
    for filename in ("IXI1", "IXI123", "IXIaaa", "foo"):
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
    assert _get_subjects_list_from_file(tmp_path / "subjects.txt") == [
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


def test_define_magnetic_field_error():
    with pytest.raises(
        ValueError,
        match=f"The hospital foo was not recognized.",
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


def test_get_marital_mapping_success(tmp_path):
    marital = pd.DataFrame(
        {
            "ID": [1, 2, 4, 3, 5],
            "MARITAL": [
                "Single",
                "Married",
                "Divorced or Separated",
                "Cohabiting",
                "Widowed",
            ],
        }
    )
    marital.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Marital Status")
    assert _get_marital_mapping(tmp_path).equals(marital.set_index("ID")["MARITAL"])


def test_get_ethnic_mapping_success(tmp_path):
    ethnic = pd.DataFrame(
        {
            "ID": [1, 4, 3, 5, 6],
            "ETHNIC": [
                "White",
                "Black or BlackBritish",
                "Asian or Asian British",
                "Chinese",
                "Other",
            ],
        }
    )
    ethnic.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Ethnicity")
    assert _get_ethnic_mapping(tmp_path).equals(ethnic.set_index("ID")["ETHNIC"])


# todo : write error for key not in excel i guess


def test_get_occupation_mapping_success(tmp_path):
    occup = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4, 5, 6, 7, 8],
            "OCCUPATION": [
                "Go out to full time employment",
                "Go out to part time employment (<25hrs)",
                "Study at college or university",
                "Full-time housework",
                "Retired",
                "Unemployed",
                "Work for pay at home",
                "Other",
            ],
        }
    )
    occup.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Occupation")
    assert _get_occupation_mapping(tmp_path).equals(occup.set_index("ID")["OCCUPATION"])


def test_get_qualification_mapping_success(tmp_path):
    qualif = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4, 5],
            "QUALIFICATION": [
                "No qualifications",
                "O - levels, GCSEs, or CSEs",
                "A - levels",
                "Further education e.g.City & Guilds / NVQs",
                "University or Polytechnic degree",
            ],
        }
    )
    qualif.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Qualification")
    assert _get_qualification_mapping(tmp_path).equals(
        qualif.set_index("ID")["QUALIFICATION"]
    )
