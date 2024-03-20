import os
import re
from unittest import mock

import pytest


def test_spm_standalone_is_available_no_env_variable():
    from clinica.utils.spm import spm_standalone_is_available

    assert not spm_standalone_is_available()


def test_spm_standalone_is_available_error(tmp_path):
    from clinica.utils.spm import spm_standalone_is_available

    not_a_folder = tmp_path / "file.txt"
    not_a_folder.touch()

    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(not_a_folder),
            "MCR_HOME": str(not_a_folder),
        },
    ):
        with pytest.raises(
            FileNotFoundError,
            match=re.escape(
                "[Error] $SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder"
            ),
        ):
            spm_standalone_is_available()


def test_spm_standalone_is_available(tmp_path):
    from clinica.utils.spm import spm_standalone_is_available

    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(tmp_path),
            "MCR_HOME": str(tmp_path),
        },
    ):
        assert spm_standalone_is_available()


def test_use_spm_standalone(tmp_path):
    from clinica.utils.spm import use_spm_standalone

    non_existent_folder = tmp_path / "foo"

    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(non_existent_folder),
            "MCR_HOME": str(non_existent_folder),
        },
    ):
        with pytest.raises(
            FileNotFoundError,
            match=re.escape(
                "$SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder"
            ),
        ):
            use_spm_standalone()


@pytest.mark.parametrize(
    "platform,expected_command",
    [
        ("Darwin", "cd /foo/bar && ./run_spm12.sh /foo/bar/baz script"),
        ("Linux", "/foo/bar/run_spm12.sh /foo/bar/baz script"),
    ],
)
def test_get_platform_dependant_matlab_command(mocker, platform, expected_command):
    from clinica.utils.spm import _get_platform_dependant_matlab_command

    mocker.patch("platform.system", return_value=platform)

    assert (
        _get_platform_dependant_matlab_command("/foo/bar", "/foo/bar/baz")
        == expected_command
    )


def test_get_platform_dependant_matlab_command_error(mocker):
    from clinica.utils.spm import _get_platform_dependant_matlab_command

    mocker.patch("platform.system", return_value="foo")

    with pytest.raises(
        SystemError,
        match="Clinica only support macOS and Linux. Your system is foo.",
    ):
        _get_platform_dependant_matlab_command("/foo/bar", "/foo/bar/baz")
