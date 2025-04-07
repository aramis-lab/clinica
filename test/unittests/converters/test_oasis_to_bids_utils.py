from pathlib import Path

import nibabel as nb
import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal


@pytest.fixture
def clinical_data_path(tmp_path: Path) -> Path:
    clinical_data_path = tmp_path / "clinical"
    _build_clinical_data(clinical_data_path)
    return clinical_data_path


def _build_clinical_data(clinical_data_path: Path) -> None:
    clinical_data_path.mkdir()

    df = pd.DataFrame(
        {
            "ID": ["OAS1_0001_MR1", "OAS1_0002_MR1"],
            "M/F": ["F", "M"],
            "Hand": ["R", "L"],
            "Age": [74, 67],
            "Educ": [2, 2],
            "SES": [3, 3],
            "MMSE": [29, 29],
            "CDR": [0, 0],
            "eTIV": [1344, 1344],
            "nWBV": [0.704, 0.645],
            "ASF": [1.306, 1.100],
            "Delay": [float("nan"), float("nan")],
        }
    )
    df.to_excel(
        clinical_data_path / "oasis_cross-sectional-5708aa0a98d82080.xlsx", index=False
    )


@pytest.fixture
def sessions_path_success(tmp_path: Path) -> Path:
    sessions_path_success = tmp_path / "spec"
    _build_spec_sessions_success(sessions_path_success)
    return sessions_path_success


def _build_spec_sessions_success(sessions_path_success: Path) -> None:
    sessions_path_success.mkdir()
    spec = pd.DataFrame(
        {
            "BIDS CLINICA": ["cdr_global", "MMS", "diagnosis", "foo"],
            "ADNI": [np.nan, np.nan, np.nan, "foo"],
            "OASIS": ["CDR", "MMSE", "CDR", np.nan],
            "OASIS location": [
                "oasis_cross-sectional-5708aa0a98d82080.xlsx",
                "oasis_cross-sectional-5708aa0a98d82080.xlsx",
                "oasis_cross-sectional-5708aa0a98d82080.xlsx",
                np.nan,
            ],
        }
    )
    spec.to_csv(sessions_path_success / "sessions.tsv", index=False, sep="\t")


@pytest.fixture
def sessions_path_error(tmp_path: Path) -> Path:
    sessions_path_error = tmp_path / "spec"
    _build_spec_sessions_error(sessions_path_error)
    return sessions_path_error


def _build_spec_sessions_error(sessions_path_error: Path) -> None:
    sessions_path_error.mkdir()
    spec = pd.DataFrame(
        {
            "BIDS CLINICA": ["foo"],
            "OASIS": ["foo"],
            "OASIS location": [
                "foo.csv",
            ],
        }
    )
    spec.to_csv(sessions_path_error / "sessions.tsv", index=False, sep="\t")


@pytest.fixture
def bids_dir(tmp_path: Path) -> Path:
    bids_dir = tmp_path / "BIDS"
    _build_bids_dir(bids_dir)
    return bids_dir


def _build_bids_dir(bids_dir: Path) -> None:
    (bids_dir / "sub-OASIS10001" / "ses-M000").mkdir(parents=True)
    (bids_dir / "sub-OASIS10001" / "ses-M006").mkdir(parents=True)
    (bids_dir / "sub-OASIS10002" / "ses-M000").mkdir(parents=True)


@pytest.fixture
def expected() -> pd.DataFrame:
    expected = {
        "sub-OASIS10001": {
            "session_id": "ses-M000",
            "cdr_global": 0,
            "MMS": 29,
            "diagnosis": "CN",
        },
        "sub-OASIS10002": {
            "session_id": "ses-M000",
            "cdr_global": 0,
            "MMS": 29,
            "diagnosis": "CN",
        },
    }

    expected = pd.DataFrame.from_dict(expected).T
    expected.index.names = ["BIDS ID"]

    return expected


def test_create_sessions_df_success(
    tmp_path,
    clinical_data_path: Path,
    sessions_path_success: Path,
    expected: pd.DataFrame,
):
    from clinica.converters.oasis_to_bids._utils import create_sessions_df

    result = create_sessions_df(
        clinical_data_path,
        sessions_path_success,
        ["sub-OASIS10001", "sub-OASIS10002"],
    )

    assert_frame_equal(expected, result, check_like=True, check_dtype=False)


def test_create_sessions_df_missing_clinical_data(
    tmp_path,
    clinical_data_path: Path,
    sessions_path_success: Path,
    expected: pd.DataFrame,
):
    from clinica.converters.oasis_to_bids._utils import create_sessions_df

    result = create_sessions_df(
        clinical_data_path,
        sessions_path_success,
        ["sub-OASIS10001", "sub-OASIS10002", "sub-OASIS10004"],
    )
    missing_line = pd.DataFrame.from_dict(
        {
            "sub-OASIS10004": {
                "session_id": "ses-M000",
                "diagnosis": "n/a",
                "cdr_global": "n/a",
                "MMS": "n/a",
            }
        }
    ).T
    missing_line.index.names = ["BIDS ID"]
    expected = pd.concat([expected, missing_line])

    assert_frame_equal(expected, result, check_like=True, check_dtype=False)


def test_create_sessions_df_file_not_found(
    tmp_path,
    clinical_data_path: Path,
    sessions_path_error: Path,
):
    from clinica.converters.oasis_to_bids._utils import create_sessions_df

    with pytest.raises(FileNotFoundError):
        create_sessions_df(
            clinical_data_path,
            sessions_path_error,
            ["sub-OASIS10001", "sub-OASIS10002"],
        )


def test_write_sessions_tsv(
    tmp_path,
    clinical_data_path: Path,
    bids_dir: Path,
    sessions_path_success: Path,
    expected: pd.DataFrame,
):
    from clinica.converters.oasis_to_bids._utils import (
        create_sessions_df,
        write_sessions_tsv,
    )

    sessions = create_sessions_df(
        clinical_data_path,
        sessions_path_success,
        ["sub-OASIS10001", "sub-OASIS10002"],
    )
    write_sessions_tsv(bids_dir, sessions)
    sessions_files = list(bids_dir.rglob("*.tsv"))

    assert len(sessions_files) == 2
    for file in sessions_files:
        assert_frame_equal(
            pd.read_csv(file, sep="\t").reset_index(drop=True),
            expected.loc[[file.parent.name]].reset_index(drop=True),
            check_like=True,
            check_dtype=False,
        )
        assert file.name == f"{file.parent.name}_sessions.tsv"


def test_write_scans_tsv(tmp_path, bids_dir: Path) -> None:
    from clinica.converters.oasis_to_bids._utils import write_scans_tsv

    image_path = (
        bids_dir
        / "sub-OASIS10001"
        / "ses-M000"
        / "anat"
        / "sub-OASIS10001_ses-M000_T1.nii.gz"
    )
    image_path.parent.mkdir(parents=True)
    image_path.touch()

    write_scans_tsv(bids_dir)

    for session_path in bids_dir.rglob("ses-M*"):
        tsv_path = list(session_path.rglob("*scans.tsv"))
        if session_path.name != "ses-M000":
            assert not tsv_path
        else:
            assert len(tsv_path) == 1
            sub = session_path.parent.name
            assert (
                tsv_path[0] == bids_dir / sub / "ses-M000" / f"{sub}_ses-M000_scans.tsv"
            )
            file = pd.read_csv(tsv_path[0], sep="\t")
            if sub == "sub-OASIS10001":
                assert file["filename"].loc[0] == f"anat/{image_path.name}"
            elif sub == "sub-OASIS10002":
                assert file.empty


def test_get_first_image(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import get_first_image

    folder_path = tmp_path / "folder"

    (folder_path).mkdir()
    (folder_path / "file_1.txt").touch()
    (folder_path / "file_2.img").touch()

    assert get_first_image(folder_path) == (folder_path / "file_2.img")


def test_get_first_image_not_found_error(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import get_first_image

    folder_path = tmp_path / "folder"

    (folder_path).mkdir()
    (folder_path / "file.txt").touch()

    with pytest.raises(
        FileNotFoundError, match=f"No file ending in .img found in {folder_path}."
    ):
        get_first_image(folder_path)


def test_get_image_with_good_orientation(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import get_image_with_good_orientation

    folder_path = tmp_path / "folder"
    folder_path.mkdir()

    # Define an empty array saved as a Nifti image
    image_path = folder_path / "image.nii.gz"
    zero_array = np.zeros((256, 256, 160, 1), dtype=np.float32)
    nifti_img = nb.Nifti1Image(zero_array, np.eye(4))
    nb.save(nifti_img, str(image_path))

    # Load the corrected Nifti image
    image = get_image_with_good_orientation(image_path)
    hdr = image.header

    affine = np.array([0, 0, -1, 80, 1, 0, 0, -128, 0, 1, 0, -128, 0, 0, 0, 1]).reshape(
        4, 4
    )
    hdr_sform = affine.astype(np.int16)

    # Check the values from header data and the image dimension
    assert hdr.get_data_shape() == (256, 256, 160)
    assert hdr.get_data_dtype() == np.int16
    assert hdr["bitpix"] == 16

    (sform_matrix, sform_code) = hdr.get_sform(coded=True)
    assert np.array_equal(hdr_sform, sform_matrix)
    assert sform_code == 1  # code 1 means scanner-based anatomical coordinates

    (qform_matrix, qform_code) = hdr.get_qform(coded=True)
    assert np.array_equal(hdr_sform, qform_matrix)
    assert qform_code == 1  # code 1 means scanner-based anatomical coordinates

    assert hdr["extents"] == 16384
    assert hdr["xyzt_units"] == 10

    assert len(image.shape) == 3  # the image dimension should be 3D


def test_get_subjects_list_from_data(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import get_subjects_list

    source_dir = tmp_path / "dataset"
    source_dir.mkdir()

    for filename in (
        "OAS1_0001_MR1",
        "OAS2_0002_MR1",
        "OAS1_0003_MR2",
        "OAS1_004_MR1",
        "foo",
    ):
        (source_dir / filename).mkdir()

    (source_dir / "OAS1_0005_MR1").touch()

    assert get_subjects_list(source_dir) == [source_dir / "OAS1_0001_MR1"]


def test_get_subjects_list_from_file(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import get_subjects_list

    source_dir = tmp_path / "dataset"
    source_dir.mkdir()

    lines = [
        "OAS1_0001_MR1\n",
        "OAS2_0002_MR1\n",
        "OAS1_0003_MR2\n",
        "OAS1_004_MR1\n",
        "foo\n",
    ]

    for filename in lines:
        (source_dir / filename[:-1]).mkdir()

    (source_dir / "OAS1_0005_MR1").touch()

    lines.append("OAS1_0005_MR1\n")

    subjects_list_dir = tmp_path / "subjects_list_dir"
    subjects_list_dir.mkdir()

    subjects_list = subjects_list_dir / "subjects_list.txt"
    subjects_list.touch()

    with open(str(subjects_list), "a") as file:
        file.writelines(lines)

    assert get_subjects_list(source_dir, subjects_list) == [
        source_dir / "OAS1_0001_MR1"
    ]


def test_mapping_diagnosis(tmp_path) -> None:
    from clinica.converters.oasis_to_bids._utils import mapping_diagnosis

    values = [0.0, 0.5, 1.0, np.nan, 2.0, 1.5, np.nan]

    mapped_values = [mapping_diagnosis(val) for val in values]

    expected_mapped_values = ["CN", "AD", "AD", "n/a", "AD", "AD", "n/a"]

    assert mapped_values == expected_mapped_values


@pytest.mark.parametrize(
    "cdr,diagnosis",
    [
        (0, "CN"),
        (12, "AD"),
        (-2, "n/a"),
        ("n/a", "n/a"),
        ("foo", "n/a"),
    ],
)
def test_convert_cdr_to_diagnosis(cdr, diagnosis: str):
    from clinica.converters.oasis_to_bids._utils import _convert_cdr_to_diagnosis

    assert diagnosis == _convert_cdr_to_diagnosis(cdr)
