import json
import os
import re
from pathlib import Path

import pytest

from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError
from clinica.utils.inputs import InvalidSubjectSession
from clinica.utils.testing_utils import (
    build_bids_directory,
    build_caps_directory,
)


@pytest.mark.parametrize(
    "input_subjects, input_sessions, to_remove, expected_subjects, expected_sessions",
    [
        (
            ["sub1", "sub1", "sub2"],
            ["ses1", "ses2", "ses1"],
            [
                InvalidSubjectSession("sub1", "ses1"),
                InvalidSubjectSession("sub3", "ses1"),
            ],
            ["sub1", "sub2"],
            ["ses2", "ses1"],
        ),
        (["sub1"], ["ses1"], [], ["sub1"], ["ses1"]),
        (
            ["sub1", "sub2", "sub2"],
            ["ses1", "ses1", "ses1"],
            [InvalidSubjectSession("sub2", "ses1")],
            ["sub1"],
            ["ses1"],
        ),
    ],
)
def test_remove_sub_ses_from_list(
    input_subjects, input_sessions, to_remove, expected_subjects, expected_sessions
):
    from clinica.utils.inputs import _remove_sub_ses_from_list

    result_subjects, result_sessions = _remove_sub_ses_from_list(
        input_subjects, input_sessions, to_remove
    )
    assert result_subjects == expected_subjects and result_sessions == expected_sessions
    assert input_subjects == input_subjects and input_sessions == input_sessions


def test_get_parent_path(tmp_path):
    from clinica.utils.inputs import _get_parent_path

    assert _get_parent_path(tmp_path / "bids" / "foo.txt") == str(tmp_path / "bids")


@pytest.mark.parametrize("extension", [".txt", ".tar.gz", ".nii.gz", ".foo.bar.baz"])
def test_get_extension(tmp_path, extension):
    from clinica.utils.inputs import _get_extension

    assert _get_extension(tmp_path / "bids" / f"foo{extension}") == extension


@pytest.mark.parametrize(
    "filename,expected_suffix",
    [
        ("foo.nii.gz", "_foo"),
        ("sub-01_bar.txt", "_bar"),
        ("sub-01_ses-M000_T1w.nii", "_T1w"),
        ("sub-01_ses-M000_run-123_pet.nii.gz", "_pet"),
    ],
)
def test_get_suffix(tmp_path, filename, expected_suffix):
    from clinica.utils.inputs import _get_suffix

    assert _get_suffix(tmp_path / "bids" / filename) == expected_suffix


@pytest.mark.parametrize(
    "filename,expected_run_number",
    [
        ("foo_run-01.txt", "01"),
        ("foo_run-00006_bar.txt.gz", "00006"),
        ("sub-01_ses-M000_run-03_pet.nii.gz", "03"),
    ],
)
def test_get_run_number(tmp_path, filename, expected_run_number):
    from clinica.utils.inputs import _get_run_number

    assert _get_run_number(str(tmp_path / "bids" / filename)) == expected_run_number


def test_select_run(tmp_path):
    from clinica.utils.inputs import _select_run

    files = [
        str(tmp_path / "bids" / "foo_run-01.txt"),
        str(tmp_path / "_run-00"),
        str(tmp_path / "bids" / "sub-01" / "sub-01_ses-M00_run-003_dwi.nii.gz"),
    ]

    assert _select_run(files) == files[-1]


def test_check_common_properties_of_files(tmp_path):
    from clinica.utils.inputs import _check_common_properties_of_files

    files = [
        tmp_path / "bids" / "foo_bar_baz.foo.bar",
        tmp_path / "bids" / "foo_pet.nii.gz",
        tmp_path / "caps" / "foo_sub-01_ses-M000_T1w.json",
        tmp_path / "foo" / "bar" / "foo_123.tar.gz",
    ]

    def first_entity_dummy_property_extractor(filename: Path) -> str:
        """Dummy extractor for testing purposes."""
        return filename.name.split("_")[0]

    assert (
        _check_common_properties_of_files(
            files,
            "first entity",
            first_entity_dummy_property_extractor,
        )
        == "foo"
    )


def test_check_common_properties_of_files_error(tmp_path):
    from clinica.utils.inputs import _check_common_properties_of_files

    files = [
        tmp_path / "bids" / "sub-01_ses-M000_pet.nii.gz",
        tmp_path / "caps" / "sub-02_ses-M123_pet.nii",
        tmp_path / "caps" / "foo.json",
    ]

    def first_letter_dummy_property_extractor(filename: Path) -> str:
        """Dummy extractor for testing purposes."""
        return filename.name[0]

    with pytest.raises(
        ValueError, match="The provided files do not share the same first letter."
    ):
        _check_common_properties_of_files(
            files,
            "first letter",
            first_letter_dummy_property_extractor,
        )


def test_get_entities(tmp_path):
    from clinica.utils.inputs import _get_entities

    files = [
        tmp_path / "bids" / "sub-01_ses-M000_run-01_pet.nii.gz",
        tmp_path / "caps" / "sub-02_ses-M000_pet.nii.gz",
        tmp_path / "foo" / "sub-01_ses-M000_run-02_pet.nii.gz",
        tmp_path / "bids" / "sub-03_ses-M000_pet.json",
    ]

    assert _get_entities(files, "_pet") == {
        "run": {"01", "02"},
        "ses": {"M000"},
        "sub": {"01", "02", "03"},
    }


@pytest.mark.parametrize(
    "files",
    [
        [],
        ["foo.txt"],
        ["foo.txt", "bar.json", "baz.png"],
        ["sub-01_ses-M000_run-01_pet.nii.gz"],
        [
            "sub-01_ses-M000_run-01_pet.nii.gz",
            "sub-01_ses-M000_run-02_pet.nii.gz",
            "sub-01_ses-M003_run-01_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_pet.nii.gz",
            "BIDS/sub-01/ses-M000/pet/sub-01_ses-M000_run-03_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_pet.nii",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-03_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_T1w.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-03_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M006_run-02_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-02_ses-M000_run-03_pet.nii.gz",
        ],
    ],
    ids=(
        "empty list of files",
        "single file without run entity",
        "no run entity anywhere",
        "single file",
        "one file has a different entity value",
        "parent paths are different",
        "extensions are different",
        "suffixes are different",
        "run numbers are not different",
        "entities are different",
    ),
)
def test_are_not_multiple_runs(files):
    from clinica.utils.inputs import _are_multiple_runs

    assert not _are_multiple_runs(files)


@pytest.mark.parametrize(
    "files",
    [
        [
            "sub-01_run-01_T1w.txt",
            "sub-01_run-02_T1w.txt",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-03_pet.nii.gz",
        ],
        [
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-01_pet.nii.gz",
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-02_desc-crop_pet.nii.gz",  # debatable that this works...
            "bids/sub-01/ses-M000/pet/sub-01_ses-M000_run-03_pet.nii.gz",
        ],
    ],
)
def test_are_multiple_runs(files):
    from clinica.utils.inputs import _are_multiple_runs

    assert _are_multiple_runs(files)


def test_insensitive_glob(tmp_path):
    from clinica.utils.inputs import insensitive_glob

    files = [
        "foo.py",
        "Bar.txt",
        "bAZ.py",
        "Fooo.PY",
        "folder_1",
        "folder_2",
        "folder_1/foo1.py",
        "folder_2/BaR2.PY",
        "folder_1/baz.txt",
    ]
    for file in files:
        d = tmp_path / file
        d.mkdir()
    python_files = insensitive_glob(str(tmp_path / "*.py"))
    assert set([Path(f).name for f in python_files]) == {"foo.py", "bAZ.py", "Fooo.PY"}
    text_files = insensitive_glob(str(tmp_path / "*.txt"))
    assert set([Path(f).name for f in text_files]) == {"Bar.txt"}
    assert len(insensitive_glob(str(tmp_path / "*.json"))) == 0
    all_python_files = insensitive_glob(str(tmp_path / "**/*.py"), recursive=True)
    assert set([Path(f).name for f in all_python_files]) == {
        "foo.py",
        "bAZ.py",
        "Fooo.PY",
        "foo1.py",
        "BaR2.PY",
    }


def test_check_bids_folder_no_subject_folder_error(tmp_path):
    from clinica.utils.inputs import check_bids_folder

    bids = tmp_path / "bids"
    bids.mkdir()
    (bids / "dataset_description.json").touch()
    (bids / "foo").mkdir()

    with pytest.raises(
        ClinicaBIDSError,
        match="Your BIDS directory does not contains a single folder whose name",
    ):
        check_bids_folder(bids)


def test_check_bids_folder(tmp_path):
    from clinica.utils.inputs import check_bids_folder

    bids = tmp_path / "bids"
    bids.mkdir()
    (bids / "dataset_description.json").touch()
    (bids / "sub-01").mkdir()

    assert check_bids_folder(bids) is None


def test_check_caps_folder(tmp_path):
    """Test function `check_caps_folder`."""
    from clinica.utils.inputs import check_caps_folder

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with pytest.raises(
        ClinicaCAPSError,
        match=re.escape(
            f"The derivative directory ({tmp_path}) you provided is missing a dataset_description.json file."
        ),
    ):
        check_caps_folder(tmp_path)
    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump(
            {"Name": "Example dataset", "BIDSVersion": "1.0.2", "CAPSVersion": "1.0.0"},
            fp,
        )
    assert check_caps_folder(tmp_path) is None
    (tmp_path / "sub-01").mkdir()
    with pytest.raises(
        ClinicaCAPSError,
        match="Your CAPS directory contains at least one folder whose name starts with 'sub-'.",
    ):
        check_caps_folder(tmp_path)


def test_find_images_path_error_no_file(tmp_path):
    """Test function `find_images_path`."""
    from clinica.utils.inputs import find_images_path

    (tmp_path / "sub-01" / "ses-M00" / "anat").mkdir(parents=True)
    errors, results = [], []

    find_images_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )

    assert len(results) == 0
    assert len(errors) == 1
    assert errors[0] == ("sub-01", "ses-M00")


def test_find_images_path_error_more_than_one_file(tmp_path):
    """Test function `find_images_path`."""
    from clinica.utils.inputs import find_images_path

    errors, results = [], []
    (tmp_path / "sub-01" / "ses-M00" / "anat" / "sub-01_ses-M00_T1w.nii.gz").mkdir(
        parents=True
    )
    (
        tmp_path / "sub-01" / "ses-M00" / "anat" / "sub-01_ses-M00_foo-bar_T1w.nii.gz"
    ).mkdir(parents=True)

    find_images_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )

    assert len(results) == 0
    assert len(errors) == 1
    assert errors[0] == ("sub-01", "ses-M00")


def test_find_images_path(tmp_path):
    """Test function `find_images_path`."""
    from clinica.utils.inputs import find_images_path

    (tmp_path / "sub-01" / "ses-M00" / "anat" / "sub-01_ses-M00_T1w.nii.gz").mkdir(
        parents=True
    )
    errors, results = [], []

    find_images_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )

    assert len(results) == 1
    assert len(errors) == 0
    assert Path(results[0]).relative_to(tmp_path) == Path(
        "sub-01/ses-M00/anat/sub-01_ses-M00_T1w.nii.gz"
    )


def test_find_images_path_multiple_runs(tmp_path):
    from clinica.utils.inputs import find_images_path

    errors, results = [], []
    (
        tmp_path
        / "sub-01"
        / "ses-M06"
        / "anat"
        / "sub-01_ses-M06_run-01_foo-bar_T1w.nii.gz"
    ).mkdir(parents=True)
    (
        tmp_path
        / "sub-01"
        / "ses-M06"
        / "anat"
        / "sub-01_ses-M06_run-02_foo-bar_T1w.nii.gz"
    ).mkdir(parents=True)

    find_images_path(
        tmp_path, "sub-01", "ses-M06", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )

    assert len(results) == 1
    assert len(errors) == 0
    assert Path(results[0]).relative_to(tmp_path) == Path(
        "sub-01/ses-M06/anat/sub-01_ses-M06_run-02_foo-bar_T1w.nii.gz"
    )


def test_check_information():
    """Test utility function `_check_information`."""
    from clinica.utils.inputs import _check_information

    with pytest.raises(
        TypeError,
        match="A dict or list of dicts must be provided for the argument 'information'",
    ):
        _check_information(42)  # noqa

    with pytest.raises(
        ValueError,
        match="'information' must contain the keys 'pattern' and 'description'",
    ):
        _check_information({})

    with pytest.raises(
        ValueError,
        match="'information' can only contain the keys 'pattern', 'description' and 'needed_pipeline'",
    ):
        _check_information({"pattern": "foo", "description": "bar", "foo": "bar"})

    with pytest.raises(
        ValueError,
        match="pattern argument cannot start with",
    ):
        _check_information({"pattern": "/foo", "description": "bar"})


def test_format_errors():
    """Test utility function `_format_errors`."""
    from clinica.utils.inputs import format_clinica_file_reader_errors

    information = {"description": "foo bar baz"}
    assert (
        format_clinica_file_reader_errors([], information)
        == "Clinica encountered 0 problem(s) while getting foo bar baz:\n"
    )
    information["needed_pipeline"] = ["pipeline_1", "pipeline_3"]
    assert format_clinica_file_reader_errors([], information) == (
        "Clinica encountered 0 problem(s) while getting foo bar baz:\n"
        "Please note that the following clinica pipeline(s) must have "
        "run to obtain these files: ['pipeline_1', 'pipeline_3']\n"
    )
    errors = [
        InvalidSubjectSession("sub1", "ses1"),
        InvalidSubjectSession("sub2", "ses1"),
        InvalidSubjectSession("sub3", "ses1"),
    ]
    assert format_clinica_file_reader_errors(errors, information) == (
        "Clinica encountered 3 problem(s) while getting foo bar baz:\n"
        "Please note that the following clinica pipeline(s) must have "
        "run to obtain these files: ['pipeline_1', 'pipeline_3']\n"
        "\t* (sub1 | ses1)\n\t* (sub2 | ses1)\n\t* (sub3 | ses1)\n"
        "Clinica could not identify which file to use (missing or too many) for these sessions. They will not be processed."
    )
    information.pop("needed_pipeline")
    assert format_clinica_file_reader_errors(errors, information) == (
        "Clinica encountered 3 problem(s) while getting foo bar baz:\n"
        "\t* (sub1 | ses1)\n\t* (sub2 | ses1)\n\t* (sub3 | ses1)\n"
        "Clinica could not identify which file to use (missing or too many) for these sessions. They will not be processed."
    )


@pytest.mark.parametrize("data_type", ["T1w", "flair"])
def test_clinica_file_reader_bids_directory(tmp_path, data_type):
    """Test reading from a BIDS directory with function `clinica_file_reader`."""
    from clinica.utils.inputs import clinica_file_reader

    config = {
        "sub-01": ["ses-M00"],
        "sub-02": ["ses-M00", "ses-M06"],
        "sub-06": ["ses-M00"],
    }

    build_bids_directory(tmp_path, config)

    desc = "T1w MRI" if data_type == "T1w" else "FLAIR T2w MRI"
    information = {
        "pattern": f"sub-*_ses-*_{data_type}.nii*",
        "description": desc,
    }

    with pytest.raises(
        ValueError,
        match="Subjects and sessions must have the same length.",
    ):
        clinica_file_reader(
            ["sub-02"],
            ["ses-M00", "ses-M06"],
            tmp_path,
            information,
            n_procs=1,
        )
    assert clinica_file_reader([], [], tmp_path, information, n_procs=1) == ([], [])
    results, errors = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, n_procs=1
    )
    assert len(results) == 1
    assert Path(results[0]).relative_to(tmp_path) == Path(
        f"sub-01/ses-M00/anat/sub-01_ses-M00_{data_type}.nii.gz"
    )
    assert not errors

    results, errors = clinica_file_reader(
        ["sub-01", "sub-02", "sub-02", "sub-06"],
        ["ses-M00", "ses-M00", "ses-M06", "ses-M00"],
        tmp_path,
        information,
        n_procs=4,
    )
    assert len(results) == 4
    assert not errors

    (
        tmp_path
        / "sub-01"
        / "ses-M00"
        / "anat"
        / f"sub-01_ses-M00_foo-bar_{data_type}.nii.gz"
    ).mkdir()
    results, errors = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, n_procs=1
    )
    assert len(results) == 0
    assert errors == [InvalidSubjectSession("sub-01", "ses-M00")]


def test_clinica_file_reader_caps_directory(tmp_path):
    """Test reading from a CAPS directory with function `clinica_file_reader`."""
    from clinica.utils.inputs import clinica_file_reader

    config = {
        "pipelines": {"t1_linear": {"uncropped_image": True}},
        "subjects": {
            "sub-01": ["ses-M00"],
            "sub-02": ["ses-M00", "ses-M06"],
            "sub-06": ["ses-M00"],
        },
    }

    build_caps_directory(tmp_path, config)

    information = {
        "pattern": "*space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz",
        "description": "T1w image registered in MNI152NLin2009cSym space using t1-linear pipeline",
        "needed_pipeline": "t1-linear",
    }

    with pytest.raises(
        ValueError,
        match="Subjects and sessions must have the same length.",
    ):
        clinica_file_reader(
            ["sub-01"],
            ["ses-M00", "ses-M06"],
            tmp_path,
            information,
            n_procs=1,
        )

    assert clinica_file_reader([], [], tmp_path, information, n_procs=1) == ([], [])

    results, errors = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, n_procs=1
    )
    assert len(results) == 1
    assert not errors

    results, errors = clinica_file_reader(
        ["sub-01", "sub-02", "sub-02", "sub-06"],
        ["ses-M00", "ses-M00", "ses-M06", "ses-M00"],
        tmp_path,
        information,
        n_procs=4,
    )
    assert len(results) == 4
    assert not errors

    (
        tmp_path
        / "subjects"
        / "sub-01"
        / "ses-M00"
        / "t1_linear"
        / "sub-01_ses-M00_foo-bar_T1w_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz"
    ).mkdir()
    results, errors = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, n_procs=1
    )
    assert len(results) == 0
    assert errors == [InvalidSubjectSession("sub-01", "ses-M00")]


def test_clinica_file_reader_dwi_dti_error(tmp_path):
    from clinica.utils.input_files import dwi_dti
    from clinica.utils.inputs import clinica_file_reader
    # todo : should be tested by check_caps_folder instead ?

    query = dwi_dti("FA", space="T1w")
    with pytest.raises(ClinicaCAPSError):
        clinica_file_reader(["sub-01"], ["ses-M000"], tmp_path, query)


def test_clinica_file_reader_dwi_dti(tmp_path):
    from clinica.utils.dwi import DTIBasedMeasure
    from clinica.utils.input_files import dwi_dti
    from clinica.utils.inputs import clinica_file_reader, clinica_list_of_files_reader

    dti_folder = (
        tmp_path
        / "subjects"
        / "sub-01"
        / "ses-M000"
        / "dwi"
        / "dti_based_processing"
        / "native_space"
    )
    dti_folder.mkdir(parents=True)
    with open(tmp_path / "dataset_description.json", "w") as fp:
        json.dump(
            {"Name": "Example dataset", "BIDSVersion": "1.0.2", "CAPSVersion": "1.0.0"},
            fp,
        )
    for measure in DTIBasedMeasure:
        (dti_folder / f"sub-01_ses-M000_space-T1w_{measure.value}.nii.gz").touch()
    query = dwi_dti("FA", space="T1w")
    found_files, _ = clinica_file_reader(["sub-01"], ["ses-M000"], tmp_path, query)
    assert found_files == [str(dti_folder / "sub-01_ses-M000_space-T1w_FA.nii.gz")]

    queries = [dwi_dti(measure) for measure in DTIBasedMeasure]
    found_files = clinica_list_of_files_reader(
        ["sub-01"], ["ses-M000"], tmp_path, queries, raise_exception=True
    )

    assert found_files == [
        [str(x)]
        for x in (
            dti_folder / f"sub-01_ses-M000_space-T1w_{m.value}.nii.gz"
            for m in DTIBasedMeasure
        )
    ]


def test_clinica_list_of_files_reader(tmp_path):
    from clinica.utils.inputs import clinica_list_of_files_reader

    config = {
        "sub-01": ["ses-M00"],
        "sub-02": ["ses-M00", "ses-M06"],
        "sub-06": ["ses-M00"],
    }

    build_bids_directory(tmp_path, config)

    information = [
        {
            "pattern": "sub-*_ses-*_t1w.nii*",
            "description": "T1w MRI",
        },
        {
            "pattern": "sub-*_ses-*_flair.nii*",
            "description": "FLAIR T2w MRI",
        },
    ]

    results = clinica_list_of_files_reader(
        ["sub-02", "sub-06", "sub-02"],
        ["ses-M00", "ses-M00", "ses-M06"],
        tmp_path,
        information,
        raise_exception=True,
    )
    assert len(results) == 2
    assert [len(r) for r in results] == [3, 3]

    (
        tmp_path
        / "sub-02"
        / "ses-M00"
        / "anat"
        / f"sub-02_ses-M00_foo-bar_flair.nii.gz"
    ).mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="Clinica faced",
    ):
        clinica_list_of_files_reader(
            ["sub-02", "sub-06", "sub-02"],
            ["ses-M00", "ses-M00", "ses-M06"],
            tmp_path,
            information,
            raise_exception=True,
        )
    results = clinica_list_of_files_reader(
        ["sub-02", "sub-06", "sub-02"],
        ["ses-M00", "ses-M00", "ses-M06"],
        tmp_path,
        information,
        raise_exception=False,
    )

    assert len(results) == 2
    assert len(results[0]) == 3
    assert len(results[1]) == 0


def test_clinica_group_reader(tmp_path):
    from clinica.utils.inputs import clinica_group_reader

    config = {
        "sub-01": ["ses-M00"],
        "sub-02": ["ses-M00", "ses-M06"],
        "sub-06": ["ses-M00"],
    }
    build_caps_directory(tmp_path, config)
    group_label = "UnitTest"
    information = {
        "pattern": os.path.join(
            f"group-{group_label}", "t1", f"group-{group_label}_template.nii*"
        ),
        "description": f"T1w template file of group {group_label}",
        "needed_pipeline": "t1-volume or t1-volume-create-dartel",
    }
    with pytest.raises(
        ClinicaCAPSError,
        match="Clinica encountered a problem while getting T1w template file of group UnitTest. No file was found",
    ):
        for raise_exception in [True, False]:
            clinica_group_reader(tmp_path, information, raise_exception=raise_exception)
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / f"group-{group_label}").mkdir()
    (tmp_path / "groups" / f"group-{group_label}" / "t1").mkdir()
    (
        tmp_path
        / "groups"
        / f"group-{group_label}"
        / "t1"
        / f"group-{group_label}_template.nii.gz"
    ).mkdir()
    result = clinica_group_reader(tmp_path, information, raise_exception=True)
    assert Path(result).relative_to(tmp_path) == Path(
        "groups/group-UnitTest/t1/group-UnitTest_template.nii.gz"
    )
    (
        tmp_path
        / "groups"
        / f"group-{group_label}"
        / "t1"
        / f"group-{group_label}_template.nii.gz2"
    ).mkdir()
    with pytest.raises(
        ClinicaCAPSError,
        match="Clinica encountered a problem while getting T1w template file of group UnitTest. 2 files were found",
    ):
        clinica_group_reader(tmp_path, information, raise_exception=True)
    result = clinica_group_reader(tmp_path, information, raise_exception=False)
    assert Path(result).stem == "group-UnitTest_template.nii"
