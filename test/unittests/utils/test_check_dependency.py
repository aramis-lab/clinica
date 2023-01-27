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


def test_check_software():
    from clinica.utils.check_dependency import _check_software
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    assert _check_software("foo") is None
    with pytest.raises(
        ClinicaMissingDependencyError,
        match=(
            "Clinica could not find foo software: the foo "
            "command is not present in your PATH environment."
        ),
    ):
        _check_software("foo", binaries=["foo"])

    assert _check_software("foo", binaries=["ls"]) is None
    with pytest.raises(
        ClinicaMissingDependencyError,
        match="Clinica could not find foo software: the FOO variable is not set.",
    ):
        _check_software("foo", binaries=["ls"], env=("FOO", "foo"))

    os.environ["FOO"] = "foo"

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="The FOO environment variable you gave is not a folder",
    ):
        _check_software("foo", binaries=["ls"], env=("FOO", "foo"))

    os.environ["FOO"] = "."
    assert _check_software("foo", binaries=["ls"], env=("FOO", "foo")) is None
    with pytest.raises(
        ClinicaMissingDependencyError,
        match=(
            "Clinica could not find foo software: the bar command is not "
            "present in your PATH environment. Foo bar baz"
        ),
    ):
        _check_software(
            "foo",
            binaries=["ls", "bar"],
            env=("FOO", "foo"),
            complementary_info="Foo bar baz",
        )
    os.environ.pop("FOO")
