import os
import re
from test.unittests.iotools.converters.adni_to_bids.modality_converters.test_adni_fmap import (
    expected,
)
from unittest import mock

import pytest
from packaging.version import Version

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


@pytest.mark.parametrize(
    "software,expected",
    [
        (ThirdPartySoftware.FREESURFER, Version("6.0.0")),
        (ThirdPartySoftware.FSL, Version("5.0.5")),
        (ThirdPartySoftware.ANTS, Version("2.5.0")),
        (ThirdPartySoftware.DCM2NIIX, Version("1.0.20240202")),
        (ThirdPartySoftware.MRTRIX, Version("3.0.3")),
        (ThirdPartySoftware.CONVERT3D, Version("1.0.0")),
        (ThirdPartySoftware.MATLAB, Version("9.2.0.556344")),
        (ThirdPartySoftware.SPM, Version("12.7219")),
        (ThirdPartySoftware.MCR, Version("9.0.1")),
        (ThirdPartySoftware.SPMSTANDALONE, Version("12.7219")),
        (ThirdPartySoftware.PETPVC, Version("0.0.0")),
    ],
)
def test_get_software_min_version_supported(software: str, expected: Version):
    from clinica.utils.check_dependency import get_software_min_version_supported

    assert get_software_min_version_supported(software) == expected


def test_get_freesurfer_version(mocker):
    from clinica.utils.check_dependency import get_software_version

    mocker.patch("nipype.interfaces.freesurfer.Info.looseversion", return_value="1.2.3")

    assert get_software_version("freesurfer") == Version("1.2.3")


def test_get_spm_version():
    from clinica.utils.check_dependency import get_software_version

    class SPMVersionMock:
        version: str = "12.6789"

    with mock.patch(
        "nipype.interfaces.spm.SPMCommand", wraps=SPMVersionMock
    ) as spm_mock:
        assert get_software_version("spm") == Version("12.6789")
        spm_mock.assert_called_once()


def test_get_spm_standalone_version(tmp_path):
    from clinica.utils.check_dependency import get_software_version

    class SPMStandaloneVersionMock:
        version: str = "13.234"

        def set_mlab_paths(self, matlab_cmd: str, use_mcr: bool):
            pass

    with mock.patch.dict(
        os.environ,
        {
            "SPM_HOME": str(tmp_path / "spm_home_mock"),
            "MCR_HOME": str(tmp_path / "mcr_home_mock"),
        },
    ):
        with mock.patch(
            "nipype.interfaces.spm.SPMCommand", wraps=SPMStandaloneVersionMock
        ) as spm_mock:
            with mock.patch.object(
                SPMStandaloneVersionMock, "set_mlab_paths", return_value=None
            ) as mock_method:
                assert get_software_version("spm standalone") == Version("13.234")
                spm_mock.set_mlab_paths.assert_called()
                mock_method.assert_called_once_with(
                    matlab_cmd=f"{tmp_path / 'spm_home_mock' / 'run_spm12.sh'} {tmp_path / 'mcr_home_mock'} script",
                    use_mcr=True,
                )


def test_get_fsl_version(mocker):
    from clinica.utils.check_dependency import get_software_version

    mocker.patch("nipype.interfaces.fsl.Info.version", return_value="3.2.1:9e026117")

    assert get_software_version("fsl") == Version("3.2.1")


mcr_version_test_suite = [
    ("2024a", Version("24.1")),
    ("2023b", Version("23.2")),
    ("2023a", Version("9.14")),
    ("2022b", Version("9.13")),
    ("2022a", Version("9.12")),
    ("2021b", Version("9.11")),
    ("2021a", Version("9.10")),
    ("2020b", Version("9.9")),
    ("2020a", Version("9.8")),
    ("2019b", Version("9.7")),
    ("2019a", Version("9.6")),
    ("2018b", Version("9.5")),
    ("2018a", Version("9.4")),
    ("2017b", Version("9.3")),
    ("2017a", Version("9.2")),
    ("2016b", Version("9.1")),
    ("2016a", Version("9.0.1")),
    ("2015b", Version("0.0.0")),
    ("1789a", Version("0.0.0")),
]


@pytest.mark.parametrize("version,expected", mcr_version_test_suite)
def test_map_mcr_release_to_version_number(version, expected):
    from clinica.utils.check_dependency import _map_mcr_release_to_version_number

    assert _map_mcr_release_to_version_number(version) == expected


@pytest.mark.parametrize("version,expected", mcr_version_test_suite)
def test_get_mcr_version(tmp_path, version, expected):
    from clinica.utils.check_dependency import get_software_version

    with mock.patch.dict(
        os.environ,
        {
            "MCR_HOME": str(tmp_path / "mcr_home_mock" / version / "v95"),
        },
    ):
        assert get_software_version("MCR") == expected
