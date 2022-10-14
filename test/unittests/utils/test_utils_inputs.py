import os
from pathlib import Path
from typing import Dict

import pytest


def build_bids_directory(directory: os.PathLike, subjects_sessions: Dict) -> None:
    """Build a fake BIDS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : PathLike
        Path to folder where the fake BIDS dataset should be created.

    subjects_sessions : Dict
        Dictionary containing the subjects and their associated sessions.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    It only adds (for now...) T1W nifti images for all subjects and sessions.
    """
    directory = Path(directory)
    data_types = {"anat"}
    suffixes = {"T1w", "flair"}
    extensions = {"nii.gz"}
    for sub, sessions in subjects_sessions.items():
        (directory / sub).mkdir()
        for ses in sessions:
            (directory / sub / ses).mkdir()
            for data_type in data_types:
                (directory / sub / ses / data_type).mkdir()
                for suffix in suffixes:
                    for extension in extensions:
                        (
                            directory
                            / sub
                            / ses
                            / data_type
                            / f"{sub}_{ses}_{suffix}.{extension}"
                        ).mkdir()


def build_caps_directory(directory: os.PathLike, subjects_sessions: Dict) -> None:
    """Build a fake CAPS dataset at the specified location following the
    specified structure.

    Parameters
    ----------
    directory : PathLike
        Path to folder where the fake CAPS dataset should be created.

    subjects_sessions : Dict
        Dictionary containing the subjects and their associated sessions.

    Notes
    -----
    This function is a simple prototype for creating fake datasets for testing.
    It only adds (for now...) partial results of the t1-linear pipeline
    for all subjects and sessions.
    """
    directory = Path(directory)
    (directory / "subjects").mkdir()
    for sub, sessions in subjects_sessions.items():
        (directory / "subjects" / sub).mkdir()
        for ses in sessions:
            (directory / "subjects" / sub / ses).mkdir()
            (directory / "subjects" / sub / ses / "t1_linear").mkdir()
            (
                directory
                / "subjects"
                / sub
                / ses
                / "t1_linear"
                / f"{sub}_{ses}_T1w_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz"
            ).mkdir()


def rmtree(f: Path):
    if f.is_file():
        f.unlink()
    else:
        for child in f.iterdir():
            rmtree(child)
        f.rmdir()


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


def test_determine_caps_or_bids(tmp_path):
    from clinica.utils.inputs import determine_caps_or_bids

    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "subjects" / "sub-01").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / "foo.txt").mkdir()
    assert not determine_caps_or_bids(tmp_path)
    rmtree(tmp_path / "subjects")
    (tmp_path / "sub-01").mkdir()
    (tmp_path / "sub-01" / "foo.txt").mkdir()
    assert determine_caps_or_bids(tmp_path)


@pytest.mark.parametrize("folder_type", ["BIDS", "CAPS"])
def test_common_checks(folder_type):
    from clinica.utils.exceptions import ClinicaBIDSError, ClinicaCAPSError
    from clinica.utils.inputs import _common_checks

    with pytest.raises(
        ValueError,
        match="Argument you provided to ",
    ):
        _common_checks(1, folder_type)  # noqa

    error = ClinicaBIDSError if folder_type == "BIDS" else ClinicaCAPSError

    with pytest.raises(
        error,
        match=f"The {folder_type} directory you gave is not a folder.",
    ):
        _common_checks(Path("fooooo"), folder_type)


def test_check_bids_folder(tmp_path):
    from clinica.utils.exceptions import ClinicaBIDSError
    from clinica.utils.inputs import check_bids_folder

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="The BIDS directory",
    ):
        check_bids_folder(tmp_path)
    rmtree(tmp_path / "subjects")
    (tmp_path / "data").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="The BIDS directory you provided is empty.",
    ):
        check_bids_folder(tmp_path / "data")
    (tmp_path / "data" / "foo").mkdir()
    with pytest.raises(
        ClinicaBIDSError,
        match="Your BIDS directory does not contains a single folder whose name",
    ):
        check_bids_folder(tmp_path / "data")
    (tmp_path / "data" / "sub-01").mkdir()
    assert check_bids_folder(tmp_path / "data") is None


def test_check_caps_folder(tmp_path):
    """Test function `check_caps_folder`."""
    from clinica.utils.exceptions import ClinicaCAPSError
    from clinica.utils.inputs import check_caps_folder

    (tmp_path / "subjects").mkdir()
    (tmp_path / "subjects" / "foo.txt").mkdir()
    assert check_caps_folder(tmp_path) is None
    (tmp_path / "sub-01").mkdir()
    with pytest.raises(
        ClinicaCAPSError,
        match="Your CAPS directory contains at least one folder whose name starts with 'sub-'.",
    ):
        check_caps_folder(tmp_path)


def test_find_sub_ses_pattern_path(tmp_path):
    """Test function `find_sub_ses_pattern_path`."""
    from clinica.utils.inputs import find_sub_ses_pattern_path

    (tmp_path / "sub-01").mkdir()
    (tmp_path / "sub-01" / "ses-M00").mkdir()
    (tmp_path / "sub-01" / "ses-M00" / "anat").mkdir()

    errors, results = [], []
    find_sub_ses_pattern_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )
    assert len(results) == 0
    assert len(errors) == 1
    assert errors[0] == "\t* (sub-01 | ses-M00): No file found\n"

    (tmp_path / "sub-01" / "ses-M00" / "anat" / "sub-01_ses-M00_T1w.nii.gz").mkdir()

    errors, results = [], []
    find_sub_ses_pattern_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )
    assert len(results) == 1
    assert len(errors) == 0
    assert Path(results[0]).relative_to(tmp_path) == Path(
        "sub-01/ses-M00/anat/sub-01_ses-M00_T1w.nii.gz"
    )

    results = []
    (
        tmp_path / "sub-01" / "ses-M00" / "anat" / "sub-01_ses-M00_foo-bar_T1w.nii.gz"
    ).mkdir()
    find_sub_ses_pattern_path(
        tmp_path, "sub-01", "ses-M00", errors, results, True, "sub-*_ses-*_t1w.nii*"
    )
    assert len(results) == 0
    assert len(errors) == 1
    assert "\t*  (sub-01 | ses-M00): More than 1 file found:" in errors[0]


def test_check_information():
    """Test utility function `_check_information`."""
    from clinica.utils.inputs import _check_information

    with pytest.raises(
        TypeError,
        match="A dict must be provided for the argument 'information'",
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
    from clinica.utils.inputs import _format_errors

    information = {"description": "foo bar baz"}
    assert (
        _format_errors([], information)
        == "Clinica encountered 0 problem(s) while getting foo bar baz:\n"
    )
    information["needed_pipeline"] = ["pipeline_1", "pipeline_3"]
    assert _format_errors([], information) == (
        "Clinica encountered 0 problem(s) while getting foo bar baz:\n"
        "Please note that the following clinica pipeline(s) must have "
        "run to obtain these files: ['pipeline_1', 'pipeline_3']\n"
    )
    errors = ["error 1: foo", "error 2: bar", "error 3: baz"]
    assert _format_errors(errors, information) == (
        "Clinica encountered 3 problem(s) while getting foo bar baz:\n"
        "Please note that the following clinica pipeline(s) must have "
        "run to obtain these files: ['pipeline_1', 'pipeline_3']\n"
        "error 1: foo\nerror 2: bar\nerror 3: baz"
    )
    information.pop("needed_pipeline")
    assert _format_errors(errors, information) == (
        "Clinica encountered 3 problem(s) while getting foo bar baz:\n"
        "error 1: foo\nerror 2: bar\nerror 3: baz"
    )


@pytest.mark.parametrize("data_type", ["T1w", "flair"])
def test_clinica_file_reader_bids_directory(tmp_path, data_type):
    """Test reading from a BIDS directory with function `clinica_file_reader`."""
    from clinica.utils.exceptions import ClinicaBIDSError
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
            raise_exception=True,
            n_procs=1,
        )
    assert clinica_file_reader(
        [], [], tmp_path, information, raise_exception=True, n_procs=1
    ) == ([], "")
    results, error_msg = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, raise_exception=True, n_procs=1
    )
    assert len(results) == 1
    assert Path(results[0]).relative_to(tmp_path) == Path(
        f"sub-01/ses-M00/anat/sub-01_ses-M00_{data_type}.nii.gz"
    )
    assert error_msg == f"Clinica encountered 0 problem(s) while getting {desc}:\n"

    results, error_msg = clinica_file_reader(
        ["sub-01", "sub-02", "sub-02", "sub-06"],
        ["ses-M00", "ses-M00", "ses-M06", "ses-M00"],
        tmp_path,
        information,
        raise_exception=True,
        n_procs=4,
    )
    assert len(results) == 4
    assert error_msg == f"Clinica encountered 0 problem(s) while getting {desc}:\n"

    (
        tmp_path
        / "sub-01"
        / "ses-M00"
        / "anat"
        / f"sub-01_ses-M00_foo-bar_{data_type}.nii.gz"
    ).mkdir()
    results, error_msg = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, raise_exception=False, n_procs=1
    )
    assert len(results) == 0
    expected_msg = (
        f"Clinica encountered 1 problem(s) while getting {desc}:\n"
        "\t*  (sub-01 | ses-M00): More than 1 file found:\n\t\t"
    )
    assert expected_msg in error_msg
    with pytest.raises(
        ClinicaBIDSError,
    ):
        clinica_file_reader(
            ["sub-01"],
            ["ses-M00"],
            tmp_path,
            information,
            raise_exception=True,
            n_procs=1,
        )


def test_clinica_file_reader_caps_directory(tmp_path):
    """Test reading from a CAPS directory with function `clinica_file_reader`."""
    from clinica.utils.exceptions import ClinicaCAPSError
    from clinica.utils.inputs import clinica_file_reader

    config = {
        "sub-01": ["ses-M00"],
        "sub-02": ["ses-M00", "ses-M06"],
        "sub-06": ["ses-M00"],
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
            raise_exception=True,
            n_procs=1,
        )

    assert clinica_file_reader(
        [], [], tmp_path, information, raise_exception=True, n_procs=1
    ) == ([], "")

    results, error_msg = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, raise_exception=True, n_procs=1
    )
    assert len(results) == 1
    expected_error_msg = (
        "Clinica encountered 0 problem(s) while getting T1w image registered "
        "in MNI152NLin2009cSym space using t1-linear pipeline:\n"
        "Please note that the following clinica pipeline(s) must have run to "
        "obtain these files: t1-linear\n"
    )
    assert error_msg == expected_error_msg

    results, error_msg = clinica_file_reader(
        ["sub-01", "sub-02", "sub-02", "sub-06"],
        ["ses-M00", "ses-M00", "ses-M06", "ses-M00"],
        tmp_path,
        information,
        raise_exception=True,
        n_procs=4,
    )
    assert len(results) == 4
    assert error_msg == expected_error_msg

    (
        tmp_path
        / "subjects"
        / "sub-01"
        / "ses-M00"
        / "t1_linear"
        / "sub-01_ses-M00_foo-bar_T1w_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz"
    ).mkdir()
    results, error_msg = clinica_file_reader(
        ["sub-01"], ["ses-M00"], tmp_path, information, raise_exception=False, n_procs=1
    )
    assert len(results) == 0
    expected_msg = (
        "Clinica encountered 1 problem(s) while getting T1w image registered "
        "in MNI152NLin2009cSym space using t1-linear pipeline:\n"
        "Please note that the following clinica pipeline(s) must have run to "
        "obtain these files: t1-linear\n"
        "\t*  (sub-01 | ses-M00): More than 1 file found:\n"
    )
    assert expected_msg in error_msg
    with pytest.raises(
        ClinicaCAPSError,
    ):
        clinica_file_reader(
            ["sub-01"],
            ["ses-M00"],
            tmp_path,
            information,
            raise_exception=True,
            n_procs=1,
        )


def test_clinica_list_of_files_reader(tmp_path):
    from clinica.utils.exceptions import ClinicaBIDSError
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
    from clinica.utils.exceptions import ClinicaCAPSError
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
