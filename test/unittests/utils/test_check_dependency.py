import os

import pytest


def test_check_dependency():
    from clinica.utils.check_dependency import is_binary_present

    assert not is_binary_present("foo")
    assert is_binary_present("ls")


def test_check_environment_variable():
    from clinica.utils.check_dependency import check_environment_variable
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="Clinica could not find foo software: the FOO variable is not set.",
    ):
        check_environment_variable("FOO", "foo")
    os.environ["FOO"] = "bar"
    with pytest.raises(
        ClinicaMissingDependencyError,
        match="The FOO environment variable you gave is not a folder",
    ):
        check_environment_variable("FOO", "foo")
    os.environ["FOO"] = "."
    assert check_environment_variable("FOO", "foo") == "."
    os.environ.pop("FOO")


@pytest.mark.parametrize(
    "binaries,env,expected_msg",
    [
        (
            ["foo"],
            None,
            "Clinica could not find foo software: the foo command is not present in your PATH environment.",
        ),
        (
            ["ls"],
            ("FOO", "foo"),
            "Clinica could not find foo software: the FOO variable is not set.",
        ),
    ],
)
def test_check_software_errors(binaries, env, expected_msg):
    from clinica.utils.check_dependency import _check_software
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    with pytest.raises(
        ClinicaMissingDependencyError,
        match=expected_msg,
    ):
        _check_software("foo", binaries=binaries, env=env)


@pytest.mark.parametrize(
    "binaries,env,complementary,expected_msg",
    [
        (
            ["ls"],
            ("FOO", "foo"),
            None,
            "The FOO environment variable you gave is not a folder",
        ),
        (
            ["ls", "bar"],
            ("FOO", "."),
            "Foo bar baz",
            (
                "Clinica could not find foo software: the bar command is "
                "not present in your PATH environment. Foo bar baz"
            ),
        ),
    ],
)
def test_check_software_env_variable_errors(binaries, env, complementary, expected_msg):
    from clinica.utils.check_dependency import _check_software
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    name, value = env
    os.environ[name] = value
    with pytest.raises(
        ClinicaMissingDependencyError,
        match=expected_msg,
    ):
        _check_software(
            "foo",
            binaries=binaries,
            env=("FOO", "foo"),
            complementary_info=complementary,
        )
    os.environ.pop(name)


def test_check_software():
    from clinica.utils.check_dependency import _check_software

    assert _check_software("foo") is None
    assert _check_software("foo", binaries=["ls"]) is None
    os.environ["FOO"] = "."
    assert _check_software("foo", binaries=["ls"], env=("FOO", "foo")) is None
    os.environ.pop("FOO")
