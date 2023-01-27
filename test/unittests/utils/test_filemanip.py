import json
import os

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_equal

from clinica.utils.exceptions import ClinicaException


@pytest.mark.parametrize(
    "dictionary, expected",
    [
        (
            {
                "TotalReadoutTime": 1,
                "EstimatedTotalReadoutTime": "",
                "PhaseEncodingSteps": "",
                "PixelBandwidth": "",
                "PhaseEncodingDirection": "j+",
                "PhaseEncodingAxis": "",
            },
            [1, "j+"],
        ),
        (
            {
                "EstimatedTotalReadoutTime": 1,
                "PhaseEncodingAxis": "j",
            },
            [1, "j+"],
        ),
        (
            {
                "PhaseEncodingSteps": 1,
                "PixelBandwidth": 0.5,
                "PhaseEncodingDirection": "j-",
            },
            [2, "j-"],
        ),
    ],
)
def test_extract_metadata_from_json_dwi(tmp_path, dictionary, expected):
    """This function tests that the outputs of `extract_metadata_from_json` are what you'd expect in the case of DWI."""
    import json

    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        handle_missing_keys_dwi,
    )

    with open(tmp_path / "metadata.json", "w") as outfile:
        json.dump(dictionary, outfile)
    assert (
        extract_metadata_from_json(
            tmp_path / "metadata.json",
            [
                "TotalReadoutTime",
                "PhaseEncodingDirection",
            ],
            handle_missing_keys_dwi,
        )
        == expected
    )


@pytest.mark.parametrize(
    "input_list, dictionary, error_type, error_log",
    [
        (
            [
                "TotalReadoutTime",
                "PhaseEncodingDirection",
            ],
            {
                "PhaseEncodingSteps": 1,
                "PixelBandwidth": 0.5,
            },
            ClinicaException,
            "Could not recover the PhaseEncodingDirection from JSON file.",
        ),
        (
            [
                "TotalReadoutTime",
                "PhaseEncodingDirection",
            ],
            {
                "PhaseEncodingDirection": "j+",
            },
            ClinicaException,
            "Could not recover the TotalReadoutTime from JSON file.",
        ),
        (
            [
                "TotalReadoutTime",
                "PhaseEncodingDirection",
            ],
            {
                "PhaseEncodingSteps": 1,
                "PixelBandwidth": 0,
                "PhaseEncodingDirection": "j+",
            },
            ValueError,
            "Pixel Bandwidth value is not valid.",
        ),
        (
            [
                "blabla",
                "PhaseEncodingDirection",
            ],
            {
                "PhaseEncodingDirection": "j+",
            },
            ValueError,
            "Could not recover the missing keys {'blabla'} from JSON file.",
        ),
    ],
)
def test_extract_metadata_from_json_dwi_errors(
    tmp_path, input_list, dictionary, error_type, error_log
):
    """This function tests that `extract_metadata_from_json` errors as expected in the case of DWI."""
    import json

    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        handle_missing_keys_dwi,
    )

    with open(tmp_path / "metadata.json", "w") as outfile:
        json.dump(dictionary, outfile)
    with pytest.raises(
        error_type,
        match=error_log,
    ):

        extract_metadata_from_json(
            tmp_path / "metadata.json",
            input_list,
            handle_missing_keys_dwi,
        )


def test_zip_nii(tmp_path):
    from clinica.utils.filemanip import zip_nii

    assert zip_nii(None) is None
    assert zip_nii(tmp_path / "foo.gz", same_dir=True) == str(tmp_path / "foo.gz")
    assert not (tmp_path / "foo.gz").exists()
    with pytest.raises(FileNotFoundError):
        zip_nii(tmp_path / "foo.nii")
    (tmp_path / "foo.nii").touch()
    assert zip_nii(tmp_path / "foo.nii", same_dir=True) == str(tmp_path / "foo.nii.gz")
    for ext in ("", ".gz"):
        assert (tmp_path / f"foo.nii{ext}").exists()


def test_unzip_nii(tmp_path):
    from clinica.utils.filemanip import unzip_nii

    assert unzip_nii(None) is None
    assert unzip_nii(tmp_path / "foo.nii") == str(tmp_path / "foo.nii")
    assert not (tmp_path / "foo.nii").exists()


def test_zip_unzip_nii(tmp_path):
    """Test that unzip(zip(x)) == x."""
    from clinica.utils.filemanip import unzip_nii, zip_nii

    p = tmp_path / "foo.nii"
    with p.open(mode="w") as f:
        f.write("Test")
    zip_nii(p, same_dir=True)
    os.remove(p)
    assert (tmp_path / "foo.nii.gz").exists()
    unzip_nii(tmp_path / "foo.nii.gz", same_dir=True)
    assert (tmp_path / "foo.nii").exists()
    with (tmp_path / "foo.nii").open() as f:
        line = f.readline()
    assert line == "Test"


@pytest.fixture
def test_image(case):
    shapes = {
        "3d": (5, 6, 7),
        "4d_dummy": (5, 6, 7, 1),
        "4d": (5, 6, 7, 8),
        "5d": (5, 6, 7, 8, 9),
    }
    data = np.zeros(shapes[case])
    data[2:4, 1:5, 3:6] = 1
    affine = np.diag((4, 3, 2, 1))
    return nib.Nifti1Image(data, affine=affine)


@pytest.mark.parametrize("case", ["3d", "4d_dummy", "4d", "5d"])
def test_load_img_3d(tmp_path, case, test_image):
    from clinica.utils.filemanip import load_volume

    with pytest.raises(
        FileNotFoundError,
        match="No such file or no access: 'foo'",
    ):
        load_volume("foo")
    filepath = tmp_path / "foo.nii.gz"
    nib.save(test_image, filepath)
    if case in ("4d", "5d"):
        with pytest.raises(
            ValueError,
            match=f"The image is not 3D but {case.upper()}.",
        ):
            load_volume(filepath)
    else:
        img2 = load_volume(filepath)
        if case == "3d":
            assert_array_equal(test_image.get_fdata(), img2.get_fdata())
        elif case == "4d_dummy":
            assert_array_equal(test_image.get_fdata().squeeze(), img2.get_fdata())


def test_save_participants_sessions(tmp_path):
    from clinica.utils.filemanip import save_participants_sessions

    with pytest.raises(
        ValueError,
        match="The number of participant IDs is not equal to the number of session IDs.",
    ):
        save_participants_sessions(["sub-01"], ["ses-M000", "ses-M006"], tmp_path)
    save_participants_sessions(["sub-01"] * 2, ["ses-M000", "ses-M006"], tmp_path)
    assert (tmp_path / "participants.tsv").exists()
    df = pd.read_csv(tmp_path / "participants.tsv", sep="\t")
    for col in ("participant_id", "session_id"):
        assert col in df.columns
    assert_array_equal(df.participant_id.values, ["sub-01"] * 2)
    assert_array_equal(df.session_id.values, ["ses-M000", "ses-M006"])


@pytest.mark.parametrize(
    "filename",
    [
        "foo",
        "sub-01.json",
        "ses-M000.nii.gz",
        "sub-01_ses-M000.json",
        "sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
    ],
)
def test_get_subject_id_error(filename):
    from clinica.utils.filemanip import get_subject_id

    with pytest.raises(
        ValueError,
        match=(
            f"Input filename {filename} is not in a BIDS or CAPS compliant format. "
            "It does not contain the subject and session information."
        ),
    ):
        get_subject_id(filename)


@pytest.mark.parametrize(
    "filename,expected",
    [
        (
            "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
            "sub-01_ses-M000",
        ),
        (
            "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
            "sub-foo_ses-bar",
        ),
    ],
)
def test_get_subject_id(filename, expected):
    from clinica.utils.filemanip import get_subject_id

    assert get_subject_id(filename) == expected


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("foo.nii.gz", "foo"),
        ("sub-01/ses-M000/sub-01_ses-M000.tar.gz", "sub-01_ses-M000"),
        ("foo/bar/baz/foo-bar_baz.niml.dset", "foo-bar_baz"),
    ],
)
def test_get_filename_no_ext(filename, expected):
    from clinica.utils.filemanip import get_filename_no_ext

    assert get_filename_no_ext(filename) == expected


def test_extract_image_ids_error():
    from clinica.utils.filemanip import extract_image_ids

    with pytest.raises(
        ValueError,
        match=(
            "Input filename foo.bar is not in a BIDS or CAPS compliant format. "
            "It does not contain the subject and session information."
        ),
    ):
        extract_image_ids(["foo.bar"])


def test_extract_image_ids():
    from clinica.utils.filemanip import extract_image_ids

    assert (
        extract_image_ids(
            [
                "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
                "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
                "sub-01_ses-M000.tar.gz",
            ]
        )
        == ["sub-01_ses-M000"] * 3
    )


def test_extract_subjects_sessions_from_filename():
    from clinica.utils.filemanip import extract_subjects_sessions_from_filename

    assert (
        extract_subjects_sessions_from_filename(
            [
                "sub-01/ses-M000/pet/sub-01_ses-M000_trc-18FAV45_pet.nii.gz",
                "foo/bar/baz/sub-foo/ses-bar/foooo/sub-01_ses-M000_foo.json",
                "sub-01_ses-M000.tar.gz",
            ]
        )
    ) == (["sub-01", "sub-01", "sub-01"], ["ses-M000", "ses-M000", "ses-M000"])


def test_extract_crash_files_from_log_file_error():
    from clinica.utils.filemanip import extract_crash_files_from_log_file

    with pytest.raises(
        ValueError,
        match="extract_crash_files_from_log_file",
    ):
        extract_crash_files_from_log_file("foo.log")


def test_read_participant_tsv(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.filemanip import read_participant_tsv

    with pytest.raises(
        ClinicaException,
        match="The TSV file you gave is not a file.",
    ):
        read_participant_tsv(tmp_path / "foo.tsv")

    df = pd.DataFrame(
        {
            "participant_id": ["sub-01", "sub-01", "sub-02"],
            "session_id": ["ses-M000", "ses-M006", "ses-M000"],
        }
    )
    df.to_csv(tmp_path / "foo.tsv", sep="\t")
    assert read_participant_tsv(tmp_path / "foo.tsv") == (
        ["sub-01", "sub-01", "sub-02"],
        ["ses-M000", "ses-M006", "ses-M000"],
    )

    for column in ("participant_id", "session_id"):
        df.drop(column, axis=1).to_csv(tmp_path / "foo.tsv", sep="\t")
        with pytest.raises(
            ClinicaException,
            match=f"The TSV file does not contain {column} column",
        ):
            read_participant_tsv(tmp_path / "foo.tsv")


def test_extract_metadata_from_json(tmp_path):
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.filemanip import extract_metadata_from_json

    data = {"foo": "foo_val", "bar": "bar_val"}

    with pytest.raises(
        FileNotFoundError,
        match="Clinica could not open the following JSON file",
    ):
        extract_metadata_from_json(tmp_path / "foo.json", ["foo", "bar"])

    with open(tmp_path / "foo.json", "w") as fp:
        json.dump(data, fp)

    with pytest.raises(
        ClinicaException,
        match="Clinica could not find the following keys in the following JSON file",
    ):
        extract_metadata_from_json(tmp_path / "foo.json", ["foo", "baz"])

    assert extract_metadata_from_json(tmp_path / "foo.json", ["foo", "bar"]) == [
        "foo_val",
        "bar_val",
    ]
