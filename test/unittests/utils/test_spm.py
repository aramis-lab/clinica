import os
import re
from pathlib import Path
from unittest import mock

import pytest


def test_spm_standalone_is_available_no_env_variable():
    from clinica.utils.spm import use_spm_standalone_if_available

    with pytest.warns(
        UserWarning,
        match=re.escape(
            "SPM standalone is not available on this system. The pipeline will try to use SPM and Matlab instead. "
            "If you want to rely on spm standalone, please make sure to set the following environment variables: "
            "$SPMSTANDALONE_HOME, and $MCR_HOME"
        ),
    ):
        assert not use_spm_standalone_if_available()


def test_spm_standalone_is_available(tmp_path, mocker):
    from clinica.utils.spm import use_spm_standalone_if_available

    mocker.patch("clinica.utils.spm._configure_spm_nipype_interface", return_value=None)
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
def test_get_platform_dependant_matlab_command(mocker, platform, expected_command):
    from clinica.utils.spm import _get_platform_dependant_matlab_command  # noqa

    mocker.patch("platform.system", return_value=platform)

    assert (
        _get_platform_dependant_matlab_command(Path("/foo/bar"), Path("/foo/bar/baz"))
        == expected_command
    )


def test_get_platform_dependant_matlab_command_error(mocker):
    from clinica.utils.spm import _get_platform_dependant_matlab_command  # noqa

    mocker.patch("platform.system", return_value="foo")

    with pytest.raises(
        SystemError,
        match="Clinica only support macOS and Linux. Your system is foo.",
    ):
        _get_platform_dependant_matlab_command(Path("/foo/bar"), Path("/foo/bar/baz"))
