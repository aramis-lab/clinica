import json
from pathlib import Path
from unittest.mock import patch

import nibabel
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal, assert_series_equal

from clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils import (
    ClinicalDataMapping,
    _define_magnetic_field,
    _find_subject_dti_data,
    _get_bids_filename_from_image_data,
    _get_img_data,
    _get_mapping,
    _get_subjects_list_from_data,
    _get_subjects_list_from_file,
    _identify_expected_modalities,
    _merge_dti,
    _padding_source_id,
    _rename_clinical_data_to_bids,
    _rename_modalities,
    _write_json_image,
    _write_subject_dti_if_exists,
    _write_subject_no_dti,
    check_modalities,
    define_participants,
    read_clinical_data,
    write_dwi_b_values,
    write_participants,
    write_scans,
    write_sessions,
)


def test_write_dwi_b_values(tmp_path):
    write_dwi_b_values(tmp_path)
    bvec_files = list(tmp_path.rglob("*.bvec"))
    bval_files = list(tmp_path.rglob("*.bval"))
    assert len(bvec_files) == 1 and bvec_files[0].name == "dwi.bvec"
    assert len(bval_files) == 1 and bval_files[0].name == "dwi.bval"


def test_get_subjects_list_from_data(tmp_path):
    for filename in ("IXI1", "IXI123", "IXIaaa", "foo"):
        (tmp_path / f"{filename}_T1w.nii.gz").touch()
    assert _get_subjects_list_from_data(tmp_path) == ["IXI123"]


def test_get_subjects_list_from_file(tmp_path):
    with open(tmp_path / "subjects.txt", "w") as f:
        f.write("IXI123\nIXI001")
    assert _get_subjects_list_from_file(tmp_path / "subjects.txt") == [
        "IXI123",
        "IXI001",
    ]


def test_define_participants_filter(tmp_path):
    for filename in ("IXI001", "IXI002", "IXI003", "IXI004"):
        (tmp_path / f"{filename}_T1w.nii.gz").touch()
    with open(tmp_path / "subjects.txt", "w") as f:
        f.write("IXI001\nIXI006")
    assert define_participants(
        data_directory=tmp_path, subjs_list_path=tmp_path / "subjects.txt"
    ) == ["IXI001"]


def test_define_participants_optional(tmp_path):
    for filename in ("IXI001", "IXI002"):
        (tmp_path / f"{filename}_T1w.nii.gz").touch()
    assert define_participants(data_directory=tmp_path) == ["IXI001", "IXI002"]


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
def test_rename_clinical_data_to_bids(input_str, expected):
    assert _rename_clinical_data_to_bids(input_str) == expected


@pytest.mark.parametrize(
    "input, expected",
    [
        ("1", "IXI001"),
        ("12", "IXI012"),
        ("123", "IXI123"),
        (1, "IXI001"),
    ],
)
def test_padding_source_id_success(input, expected):
    assert _padding_source_id(input) == expected


def test_padding_source_id_error():
    with pytest.raises(
        ValueError,
        match=f"The source id 1234 has more than 3 digits while IXI"
        f"source ids are expected to be between 1 and 3 digits.",
    ):
        _padding_source_id("1234")


@pytest.mark.parametrize("input", ["IXI_name.xls", "IXI_format.csv"])
def test_read_clinical_data_error(tmp_path, input):
    (tmp_path / input).touch()
    with pytest.raises(
        FileNotFoundError,
        match=f"Clinical data stored in the folder {tmp_path} is expected to be an excel file named 'IXI.xls'. "
        f"In case the file downloaded from the IXI website changed format, please do not hesitate to report to us !",
    ):
        read_clinical_data(tmp_path)


def clinical_data_builder(tmp_path: Path) -> None:
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

    subject = pd.DataFrame(
        {
            "IXI_ID": ["1"],
            "STUDY_DATE": ["2024-08-23"],
            "SEX_ID (1=m, 2=f)": [2],
            "ETHNIC_ID": [6],
            "MARITAL_ID": [1],
            "OCCUPATION_ID": [8],
            "QUALIFICATION_ID": [3],
            "DOB": ["2000-01-01"],
            "WEIGHT": [80],
        }
    )
    with pd.ExcelWriter(tmp_path / "IXI.xls") as writer:
        subject.to_excel(writer, index=False)
        qualif.to_excel(writer, sheet_name="Qualification", index=False)
        occup.to_excel(writer, sheet_name="Occupation", index=False)
        ethnic.to_excel(writer, sheet_name="Ethnicity", index=False)
        marital.to_excel(writer, sheet_name="Marital Status", index=False)


def test_read_clinical_data_success(tmp_path):
    clinical_data_builder(tmp_path)
    assert (
        read_clinical_data(tmp_path)
        .eq(formatted_clinical_data_builder().drop("participant_id", axis=1))
        .all()
        .all()
    )


def test_merge_dti(tmp_path):
    im1 = nibabel.Nifti1Image(
        np.empty(shape=(256, 156, 256), dtype=np.float64), np.eye(4)
    )
    im1.to_filename(tmp_path / "im1.nii.gz")
    im2 = nibabel.Nifti1Image(
        np.empty(shape=(256, 156, 256), dtype=np.float64), np.eye(4)
    )
    im2.to_filename(tmp_path / "im2.nii.gz")
    merged = _merge_dti([tmp_path / "im1.nii.gz", tmp_path / "im2.nii.gz"])
    assert type(merged) == nibabel.Nifti1Image
    assert merged.shape[-1] == 2


def test_write_dti_success(tmp_path):
    im1 = nibabel.Nifti1Image(
        np.empty(shape=(256, 156, 256), dtype=np.float64), np.eye(4)
    )
    im1.to_filename(tmp_path / "IXI001-Guys-1234-DTI-00.nii.gz")
    im2 = nibabel.Nifti1Image(
        np.empty(shape=(256, 156, 256), dtype=np.float64), np.eye(4)
    )
    im2.to_filename(tmp_path / "IXI001-Guys-1234-DTI-01.nii.gz")

    _write_subject_dti_if_exists(
        bids_path=tmp_path, subject="IXI001", data_directory=tmp_path
    )
    dti_image = list(tmp_path.rglob(pattern="*dwi.nii.gz"))
    dti_json = list(tmp_path.rglob(pattern="*dwi.json"))

    assert (
        len(dti_image) == 1
        and dti_image[0]
        == tmp_path
        / "sub-IXI001"
        / "ses-M000"
        / "dwi"
        / "sub-IXI001_ses-M000_dwi.nii.gz"
    )
    assert (
        len(dti_json) == 1
        and dti_json[0]
        == tmp_path / "sub-IXI001" / "ses-M000" / "dwi" / "sub-IXI001_ses-M000_dwi.json"
    )


def test_write_dti_empty(tmp_path):
    _write_subject_dti_if_exists(
        bids_path=tmp_path, subject="IXI001", data_directory=tmp_path
    )
    dti_files = list(tmp_path.rglob(pattern="*dwi.nii.gz"))
    assert not dti_files


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
    assert assert_frame_equal(_get_img_data(tmp_path), input) is None


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
    assert (
        assert_series_equal(
            _get_mapping(tmp_path, ClinicalDataMapping.MARITAL),
            marital.set_index("ID")["MARITAL"],
        )
        is None
    )


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
    assert (
        assert_series_equal(
            _get_mapping(tmp_path, ClinicalDataMapping.ETHNIC),
            ethnic.set_index("ID")["ETHNIC"],
        )
        is None
    )


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
    assert (
        assert_series_equal(
            _get_mapping(tmp_path, ClinicalDataMapping.OCCUPATION),
            occup.set_index("ID")["OCCUPATION"],
        )
        is None
    )


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
    assert (
        assert_series_equal(
            _get_mapping(tmp_path, ClinicalDataMapping.QUALIFICATION),
            qualif.set_index("ID")["QUALIFICATION"],
        )
        is None
    )


def test_get_mapping_fileerror(tmp_path):
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
    qualif.to_excel(excel_writer=tmp_path / "IXI2.xls", sheet_name="Qualification")

    with pytest.raises(
        FileNotFoundError,
        match=f"Clinical data stored in the folder {tmp_path} is expected to be an excel file named 'IXI.xls'. "
        f"In case the file downloaded from the IXI website changed format, please do not hesitate to report to us !",
    ):
        _get_mapping(tmp_path, ClinicalDataMapping.QUALIFICATION)


def test_get_mapping_keyerror(tmp_path):
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
    qualif.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Qualification_error")

    with pytest.raises(
        ValueError,
        match=f"Qualification mapping is expected to be contained in a sheet called Qualification coming from the clinical data excel. "
        f"Possibilities are supposed to be described in a QUALIFICATION column associated to keys from the 'ID' column. "
        f"In case the file downloaded from the IXI website changed format, please do not hesitate to report to us !",
    ):
        _get_mapping(tmp_path, ClinicalDataMapping.QUALIFICATION)


def test_get_mapping_valueerror(tmp_path):
    qualif = pd.DataFrame(
        {
            "ID": [1, 2, 3, 4, 5],
            "QUALIFICATION_error": [
                "No qualifications",
                "O - levels, GCSEs, or CSEs",
                "A - levels",
                "Further education e.g.City & Guilds / NVQs",
                "University or Polytechnic degree",
            ],
        }
    )
    qualif.to_excel(excel_writer=tmp_path / "IXI.xls", sheet_name="Qualification")

    with pytest.raises(
        ValueError,
        match=f"Qualification mapping is expected to be contained in a sheet called Qualification coming from the clinical data excel. "
        f"Possibilities are supposed to be described in a QUALIFICATION column associated to keys from the 'ID' column. "
        f"In case the file downloaded from the IXI website changed format, please do not hesitate to report to us !",
    ):
        _get_mapping(tmp_path, ClinicalDataMapping.QUALIFICATION)


def test_write_json_image(tmp_path):
    _write_json_image(tmp_path / "test.json", hospital="Guys", field="1.5")
    with open(tmp_path / "test.json", "r") as f:
        data = json.load(f)
    assert data["InstitutionName"] == "Guys"
    assert data["MagneticFieldStrength (T)"] == "1.5"


def test_write_subject_no_dti(tmp_path):
    df = image_dataframe_builder(tmp_path)
    bids_dir = tmp_path / "BIDS"
    file_path = (
        bids_dir
        / f"sub-{df['subject'][0]}"
        / "ses-M000"
        / "anat"
        / f"sub-{df['subject'][0]}_ses-M000_T1w"
    )
    _write_subject_no_dti(df, bids_dir)
    json_files = list(bids_dir.rglob(f"sub-{df['subject'][0]}*.json"))
    nii_files = list(bids_dir.rglob(f"sub-{df['subject'][0]}*.nii.gz"))
    assert len(json_files) == 1 and json_files[0] == Path(f"{file_path}.json")
    assert len(nii_files) == 1 and nii_files[0] == Path(f"{file_path}.nii.gz")


def test_write_subject_no_dti_empty(tmp_path):
    bids_dir = tmp_path / "BIDS"
    bids_dir.mkdir()
    _write_subject_no_dti(pd.DataFrame(), bids_dir)
    assert not list(bids_dir.iterdir())


def test_find_subject_dti_data(tmp_path):
    (tmp_path / "IXI001-Guys-1234-T1.nii.gz").touch()
    (tmp_path / "IXI001-Guys-1234-DTI").touch()
    (tmp_path / "IXI001-Guys-DTI.nii.gz").touch()
    tmp_image = tmp_path / "IXI001-Guys-1234-DTI-01.nii.gz"
    tmp_image.touch()
    list_dti = _find_subject_dti_data(data_directory=tmp_path, subject="IXI001")
    assert len(list_dti) == 1 and list_dti[0] == tmp_image


def test_identify_expected_modalities(tmp_path):
    (tmp_path / "IXI-DTI").mkdir()
    (tmp_path / "IXIdti").mkdir()
    (tmp_path / "foo-bar").mkdir()
    assert _identify_expected_modalities(tmp_path) == ["DTI"]


def test_write_scans_not_empty(tmp_path):
    (tmp_path / "sub-IXI001" / "ses-M000" / "anat").mkdir(parents=True)
    (tmp_path / "sub-IXI001" / "ses-M000" / "anat" / "sub-IXI001_T1w.nii.gz").touch()
    write_scans(tmp_path, participant="IXI001")
    tsv_files = list(tmp_path.rglob("*.tsv"))
    file_path = tmp_path / "sub-IXI001" / "ses-M000" / "sub-IXI001_ses-M000_scans.tsv"
    assert len(tsv_files) == 1 and tsv_files[0] == file_path
    assert (
        assert_frame_equal(
            pd.read_csv(file_path, sep="\t"),
            pd.DataFrame({"filename": ["anat/sub-IXI001_T1w.nii.gz"]}),
        )
        is None
    )


def formatted_clinical_data_builder() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source_id": ["IXI001"],
            "participant_id": ["sub-IXI001"],
            "session_id": ["ses-M000"],
            "acq_time": ["2024-08-23"],
            "sex": ["female"],
            "ethnicity": ["Other"],
            "marital status": ["Single"],
            "occupation": ["Other"],
            "qualification": ["A - levels"],
            "date of birth": ["2000-01-01"],
            "weight": [80],
        }
    )


def test_write_sessions(tmp_path):
    clinical = formatted_clinical_data_builder()
    (tmp_path / "sub-IXI001").mkdir()
    write_sessions(tmp_path, clinical, "IXI001")
    tsv_files = list(tmp_path.rglob("*.tsv"))
    file_path = tmp_path / "sub-IXI001" / "sub-IXI001_sessions.tsv"
    assert len(tsv_files) == 1 and tsv_files[0] == file_path
    assert (
        assert_frame_equal(
            pd.read_csv(file_path, sep="\t"),
            clinical[["source_id", "session_id", "acq_time"]],
        )
        is None
    )


def test_write_participants(tmp_path):
    clinical = formatted_clinical_data_builder()
    expected = clinical.copy()
    write_participants(tmp_path, clinical, ["IXI001", "IXI002"])
    expected.drop(["acq_time", "session_id"], axis=1, inplace=True)
    expected = pd.concat(
        [expected, pd.DataFrame({col: ["n/a"] for col in expected.columns})]
    ).reset_index(drop=True)
    expected.loc[1, "source_id"] = "IXI002"
    expected.loc[0, "weight"] = str(expected.loc[0, "weight"])
    tsv_files = list(tmp_path.rglob("*.tsv"))
    assert len(tsv_files) == 1 and tsv_files[0] == tmp_path / "participants.tsv"
    assert (
        assert_frame_equal(
            pd.read_csv(tmp_path / "participants.tsv", sep="\t", na_filter=False),
            expected,
        )
        is None
    )


@patch("clinica.iotools.converters.ixi_to_bids.ixi_to_bids_utils.cprint")
def test_check_modalities(mock_cprint, tmp_path):
    (tmp_path / "IXI-DTI").mkdir()
    (tmp_path / "IXI-DTI" / "IXI001-DTI-00.nii.gz").touch()
    (tmp_path / "IXI-T1").mkdir()
    (tmp_path / "IXI-T1" / "IXI001-T1.nii.gz").touch()
    (tmp_path / "IXI-T1" / "IXI002-T1.nii.gz").touch()
    (tmp_path / "IXIT1").mkdir()

    message = (
        f"Modalities : dti , T1w were identified inside {tmp_path} for conversion.\n"
        f"Some subjects do not have data for the following modalities :\n"
        f"IXI002 : dti\n"
    )

    check_modalities(tmp_path, ["IXI001", "IXI002"])
    mock_cprint.assert_called_once_with(message)
