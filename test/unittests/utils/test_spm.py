import os
import re
from pathlib import Path
from unittest import mock

import pytest

from clinica.utils.exceptions import ClinicaMissingDependencyError
from clinica.utils.spm import SPMTissue


@pytest.mark.parametrize("index", [-1, 0, 7, 10.2, 3.3, "foo", "3", "", None])
def test_get_spm_tissue_from_index_error(index):
    from clinica.utils.spm import get_spm_tissue_from_index

    with pytest.raises(
        ValueError,
        match=f"No SPM tissue matching index {index}.",
    ):
        get_spm_tissue_from_index(index)


@pytest.mark.parametrize(
    "index,expected,expected_value",
    [
        (1, SPMTissue.GRAY_MATTER, "graymatter"),
        (2, SPMTissue.WHITE_MATTER, "whitematter"),
        (3, SPMTissue.CSF, "csf"),
        (4, SPMTissue.BONE, "bone"),
        (5, SPMTissue.SOFT_TISSUE, "softtissue"),
        (6, SPMTissue.BACKGROUND, "background"),
    ],
)
def test_get_spm_tissue_from_index_error(index, expected, expected_value):
    from clinica.utils.spm import get_spm_tissue_from_index

    assert get_spm_tissue_from_index(index) == expected
    assert get_spm_tissue_from_index(index).value == expected_value


def test_spm_standalone_is_available_no_env_variable_error():
    from clinica.utils.spm import use_spm_standalone_if_available

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="Clinica could not find spm software: the SPM_HOME variable is not set.",
    ):
        use_spm_standalone_if_available()


def test_spm_standalone_is_available_warning(tmp_path):
    from clinica.utils.spm import use_spm_standalone_if_available

    (tmp_path / "spm_home_folder").mkdir()
    with pytest.warns(
        UserWarning,
        match=re.escape(
            "SPM standalone is not available on this system. The pipeline will try to use SPM and Matlab instead. "
            "If you want to rely on spm standalone, please make sure to set the following environment variables: "
            "$SPMSTANDALONE_HOME, and $MCR_HOME"
        ),
    ):
        with mock.patch.dict(
            os.environ, {"SPM_HOME": str(tmp_path / "spm_home_folder")}
        ):
            assert not use_spm_standalone_if_available()


def test_spm_standalone_is_available(tmp_path, mocker):
    from clinica.utils.spm import use_spm_standalone_if_available

    mocker.patch(
        "clinica.utils.spm.configure_nipype_interface_to_work_with_spm_standalone",
        return_value=None,
    )
    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(tmp_path),
            "MCR_HOME": str(tmp_path),
        },
    ):
        assert use_spm_standalone_if_available()


def test_use_spm_standalone_if_available_error(tmp_path):
    from clinica.utils.exceptions import ClinicaEnvironmentVariableError
    from clinica.utils.spm import use_spm_standalone_if_available

    non_existent_folder = tmp_path / "foo"
    non_existent_folder.touch()

    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(non_existent_folder),
            "MCR_HOME": str(non_existent_folder),
        },
    ):
        with pytest.raises(
            ClinicaEnvironmentVariableError,
            match=re.escape(
                "The SPMSTANDALONE_HOME environment variable you gave is not a folder"
            ),
        ):
            use_spm_standalone_if_available()


@pytest.mark.parametrize(
    "platform,expected_command",
    [
        ("Darwin", "cd /foo/bar && ./run_spm12.sh /foo/bar/baz script"),
        ("Linux", "/foo/bar/run_spm12.sh /foo/bar/baz script"),
    ],
)
def test_get_platform_dependant_matlab_command_for_spm_standalone(
    mocker, platform, expected_command
):
    from clinica.utils.spm import (
        _get_platform_dependant_matlab_command_for_spm_standalone,
    )

    mocker.patch("platform.system", return_value=platform)
    mocker.patch(
        "clinica.utils.spm._get_real_spm_standalone_file", return_value="run_spm12.sh"
    )

    assert (
        _get_platform_dependant_matlab_command_for_spm_standalone(
            Path("/foo/bar"), Path("/foo/bar/baz")
        )
        == expected_command
    )


def test_get_platform_dependant_matlab_command_for_spm_standalone_error(mocker):
    from clinica.utils.spm import (
        _get_platform_dependant_matlab_command_for_spm_standalone,
    )

    mocker.patch("platform.system", return_value="foo")
    mocker.patch(
        "clinica.utils.spm._get_real_spm_standalone_file", return_value="run_spm12.sh"
    )

    with pytest.raises(
        SystemError,
        match="Clinica only support macOS and Linux. Your system is foo.",
    ):
        _get_platform_dependant_matlab_command_for_spm_standalone(
            Path("/foo/bar"), Path("/foo/bar/baz")
        )


def test_get_real_spm_standalone_file_no_file_error(tmp_path):
    from clinica.utils.spm import _get_real_spm_standalone_file

    with pytest.raises(
        FileNotFoundError,
        match=f"There is no or several 'run_spmXX.sh' in your SPMSTANDALONE_HOME {tmp_path} : ",
    ):
        _get_real_spm_standalone_file(tmp_path)


def test_get_real_spm_standalone_file_several_files_error(tmp_path):
    from clinica.utils.spm import _get_real_spm_standalone_file

    (tmp_path / "run_spm1.sh").touch()
    (tmp_path / "run_spm2.sh").touch()

    with pytest.raises(
        FileNotFoundError,
        match=f"There is no or several 'run_spmXX.sh' in your SPMSTANDALONE_HOME {tmp_path} : {tmp_path/"run_spm1.sh"} ; {tmp_path/"run_spm2.sh"}",
    ):
        _get_real_spm_standalone_file(tmp_path)


@pytest.mark.parametrize(
    "spm_file",
    [
        "run_spm.sh",
        "run_spm12.sh",
        "run_spm24.sh",
        "run_spm25.01.sh",
        "run_spmrc24.sh",
        "run_spm_25.sh",
    ],
)
def test_get_real_spm_standalone_file_error(tmp_path, spm_file):
    from clinica.utils.spm import _get_real_spm_standalone_file

    (tmp_path / spm_file).touch()
    assert spm_file == _get_real_spm_standalone_file(tmp_path)
