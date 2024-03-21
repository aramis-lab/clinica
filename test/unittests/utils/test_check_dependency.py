import os
import re
from unittest import mock

import pytest

from clinica.utils.check_dependency import (
    SoftwareEnvironmentVariable,
    ThirdPartySoftware,
)
from clinica.utils.exceptions import (
    ClinicaEnvironmentVariableError,
    ClinicaMissingDependencyError,
)


def test_is_binary_present():
    from clinica.utils.check_dependency import is_binary_present

    assert not is_binary_present("foo")
    assert is_binary_present("ls")


def test_check_binary():
    from clinica.utils.check_dependency import check_binary

    check_binary("ls")
    with pytest.raises(
        ClinicaMissingDependencyError,
        match=(
            "Command foo is unknown on your system. "
            "Please verify that you have correctly installed the corresponding software."
        ),
    ):
        check_binary("foo")


def test_check_environment_variable_not_set_error():
    from clinica.utils.check_dependency import check_environment_variable

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="Clinica could not find fsl software: the FOO variable is not set.",
    ):
        check_environment_variable(
            SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL)
        )


def test_check_environment_variable_not_a_folder_error(tmp_path):
    from clinica.utils.check_dependency import check_environment_variable

    not_a_folder = tmp_path / "foo.txt"
    not_a_folder.touch()

    with mock.patch.dict(os.environ, {"FOO": str(not_a_folder)}):
        with pytest.raises(
            ClinicaEnvironmentVariableError,
            match="The FOO environment variable you gave is not a folder",
        ):
            check_environment_variable(
                SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL)
            )


def test_check_environment_variable(tmp_path):
    from clinica.utils.check_dependency import check_environment_variable

    with mock.patch.dict(os.environ, {"FOO": str(tmp_path)}):
        assert (
            check_environment_variable(
                SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL)
            )
            == tmp_path
        )


@pytest.mark.parametrize(
    "binaries,env,expected_msg",
    [
        (
            ("foo",),
            None,
            "Clinica could not find fsl software: the foo command is not present in your PATH environment.",
        ),
        (
            ("ls",),
            (SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL),),
            "Clinica could not find fsl software: the FOO variable is not set.",
        ),
    ],
)
def test_check_software_errors(binaries, env, expected_msg):
    from clinica.utils.check_dependency import _check_software  # noqa

    with pytest.raises(
        ClinicaMissingDependencyError,
        match=expected_msg,
    ):
        _check_software(
            ThirdPartySoftware.FSL, binaries=binaries, environment_variables=env
        )


@pytest.mark.parametrize(
    "binaries,env,complementary,error_class,expected_msg",
    [
        (
            ("ls",),
            ("FOO", "foo"),
            None,
            ClinicaEnvironmentVariableError,
            "The FOO environment variable you gave is not a folder",
        ),
        (
            ("ls", "bar"),
            ("FOO", ""),
            "Foo bar baz",
            ClinicaMissingDependencyError,
            (
                "Clinica could not find fsl software: the bar command is "
                "not present in your PATH environment. Foo bar baz"
            ),
        ),
    ],
)
def test_check_software_env_variable_errors(
    tmp_path, binaries, env, complementary, error_class, expected_msg
):
    from clinica.utils.check_dependency import _check_software  # noqa

    with mock.patch.dict(os.environ, {env[0]: str(tmp_path / env[1])}):
        with pytest.raises(error_class, match=expected_msg):
            _check_software(
                ThirdPartySoftware.FSL,
                binaries=binaries,
                environment_variables=(
                    SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL),
                ),
                complementary_info=complementary,
            )


def test_check_software(tmp_path):
    from clinica.utils.check_dependency import _check_software  # noqa

    assert _check_software("foo") is None
    assert _check_software("foo", binaries=("ls",)) is None

    (tmp_path / "foo").mkdir()
    with mock.patch.dict(os.environ, {"FOO": str(tmp_path / "foo")}):
        assert (
            _check_software(
                ThirdPartySoftware.FSL,
                binaries=("ls",),
                environment_variables=(
                    SoftwareEnvironmentVariable("FOO", ThirdPartySoftware.FSL),
                ),
            )
            is None
        )


def test_check_spm_error(tmp_path):
    from clinica.utils.check_dependency import _check_spm  # noqa

    with pytest.raises(
        ClinicaMissingDependencyError,
        match=re.escape(
            "Clinica could not find the SPM software (regular or standalone)."
        ),
    ):
        _check_spm()


@pytest.mark.parametrize("env_name", ["SPMSTANDALONE_HOME", "MCR_HOME"])
def test_check_spm_standalone_error_incomplete_config(tmp_path, env_name):
    from clinica.utils.check_dependency import _check_spm  # noqa

    with mock.patch.dict(os.environ, {env_name: str(tmp_path)}):
        with pytest.raises(
            ClinicaMissingDependencyError,
            match=re.escape(
                "Clinica could not find the SPM software (regular or standalone)."
            ),
        ):
            _check_spm()


def test_check_spm_standalone(tmp_path):
    from clinica.utils.check_dependency import _check_spm  # noqa

    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(tmp_path),
            "MCR_HOME": str(tmp_path),
        },
    ):
        _check_spm()


def test_check_spm_alone_error_matlab_not_installed(tmp_path):
    from clinica.utils.check_dependency import _check_spm  # noqa

    with mock.patch.dict(os.environ, {"SPM_HOME": str(tmp_path)}):
        with pytest.raises(
            ClinicaMissingDependencyError,
            match=re.escape(
                "[Error] Clinica could not find spm software: the matlab "
                "command is not present in your PATH environment."
            ),
        ):
            _check_spm()


def test_check_spm_alone(tmp_path, mocker):
    from clinica.utils.check_dependency import _check_spm  # noqa

    mocker.patch("clinica.utils.check_dependency.is_binary_present", return_value=True)
    with mock.patch.dict(os.environ, {"SPM_HOME": str(tmp_path)}):
        _check_spm()
