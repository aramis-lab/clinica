import os
import re
from functools import partial
from test.unittests.iotools.converters.adni_to_bids.modality_converters.test_adni_fmap import (
    expected,
)
from typing import Optional, Union
from unittest import mock

import pytest
from packaging.specifiers import SpecifierSet
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


def test_get_spm_version(tmp_path):
    from clinica.utils.check_dependency import get_software_version

    class SPMVersionMock:
        version: str = "12.6789"

    spm_home_mock = tmp_path / "spm_home"
    spm_home_mock.mkdir()

    with mock.patch.dict(
        os.environ,
        {
            "SPM_HOME": str(spm_home_mock),
        },
    ):
        with mock.patch(
            "nipype.interfaces.spm.SPMCommand", wraps=SPMVersionMock
        ) as spm_mock:
            assert get_software_version("spm") == Version("12.6789")
            spm_mock.assert_called_once()


@pytest.mark.parametrize("platform", ["linux", "darwin"])
def test_get_spm_standalone_version(tmp_path, mocker, platform):
    from clinica.utils.check_dependency import get_software_version

    class SPMStandaloneVersionMock:
        version: str = "13.234"

        def set_mlab_paths(self, matlab_cmd: str, use_mcr: bool):
            pass

    spm_standalone_home_mock = tmp_path / "spm_standalone_home"
    spm_standalone_home_mock.mkdir()
    mcr_home_mock = tmp_path / "mcr_home"
    mcr_home_mock.mkdir()

    mocker.patch("platform.system", return_value=platform)
    with mock.patch.dict(
        os.environ,
        {
            "SPMSTANDALONE_HOME": str(spm_standalone_home_mock),
            "MCR_HOME": str(mcr_home_mock),
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
                    matlab_cmd=(
                        f"cd {tmp_path / 'spm_standalone_home'} && ./run_spm12.sh {tmp_path / 'mcr_home'} script"
                        if platform == "darwin"
                        else f"{tmp_path / 'spm_standalone_home' / 'run_spm12.sh'} {tmp_path / 'mcr_home'} script"
                    ),
                    use_mcr=True,
                )


def test_get_fsl_version(mocker):
    from clinica.utils.check_dependency import get_software_version

    mocker.patch("nipype.interfaces.fsl.Info.version", return_value="3.2.1:9e026117")

    assert get_software_version("fsl") == Version("3.2.1")


def ants_version_mock(executable: str, two_dashes: bool = True) -> str:
    return "ANTs Version: 2.5.0.post7-g46ab4e7\nCompiled: Sep  8 2023 14:35:39"


def test_get_ants_version():
    from clinica.utils.check_dependency import get_software_version

    with mock.patch(
        "clinica.utils.check_dependency._run_command", wraps=ants_version_mock
    ) as ants_mock:
        assert get_software_version("ants") == Version("2.5.0")
        ants_mock.assert_called_once_with("antsRegistration", two_dashes=True)


def dcm2niix_version_mock(executable: str, two_dashes: bool = True) -> str:
    return (
        "Compression will be faster with 'pigz' installed http://macappstore.org/pigz/\n"
        "Chris Rorden's dcm2niiX version v1.0.20240202  Clang15.0.0 ARM (64-bit MacOS)\n"
        "v1.0.20240202"
    )


def test_get_dcm2niix_version():
    from clinica.utils.check_dependency import get_software_version

    with mock.patch(
        "clinica.utils.check_dependency._run_command", wraps=dcm2niix_version_mock
    ) as dcm2niix_mock:
        assert get_software_version("dcm2niix") == Version("1.0.20240202")
        dcm2niix_mock.assert_called_once_with("dcm2niix", two_dashes=True)


def mrtrix_version_mock(executable: str, two_dashes: bool = True) -> str:
    return (
        "== mrtransform 3.0.3 ==\n"
        "64 bit release version with nogui, built Oct 22 2021, using Eigen 3.3.7\n"
        "Author(s): J-Donald Tournier (jdtournier@gmail.com) and David Raffelt "
        "(david.raffelt@florey.edu.au) and Max Pietsch (maximilian.pietsch@kcl.ac.uk)\n"
        "Copyright (c) 2008-2021 the MRtrix3 contributors.\n\n"
        "This Source Code Form is subject to the terms of the Mozilla Public\n"
        "License, v. 2.0. If a copy of the MPL was not distributed with this\n"
        "file, You can obtain one at http://mozilla.org/MPL/2.0/.\n\n"
        "Covered Software is provided under this License on an 'as is'\n"
        "basis, without warranty of any kind, either expressed, implied, or\n"
        "statutory, including, without limitation, warranties that the\n"
        "Covered Software is free of defects, merchantable, fit for a\n"
        "particular purpose or non-infringing.\n"
        "See the Mozilla Public License v. 2.0 for more details.\n\n"
        "For more details, see http://www.mrtrix.org/."
    )


def test_get_mrtrix_version():
    from clinica.utils.check_dependency import get_software_version

    with mock.patch(
        "clinica.utils.check_dependency._run_command", wraps=mrtrix_version_mock
    ) as mrtrix_mock:
        assert get_software_version("mrtrix") == Version("3.0.3")
        mrtrix_mock.assert_called_once_with("mrtransform", two_dashes=True)


def convert3d_version_mock(executable: str, two_dashes: bool = True) -> str:
    return "Version 1.0.0"


def test_get_convert3d_version():
    from clinica.utils.check_dependency import get_software_version

    with mock.patch(
        "clinica.utils.check_dependency._run_command", wraps=convert3d_version_mock
    ) as c3d_mock:
        assert get_software_version("convert3d") == Version("1.0.0")
        c3d_mock.assert_called_once_with("c3d", two_dashes=False)


def matlab_version_mock() -> str:
    return (
        "\x1b[?1h\x1b=\n                                                                                          "
        "< M A T L A B (R) >\n                                                                                "
        "Copyright 1984-2017 The MathWorks, Inc.\n"
        "                                                                                 R2017b (9.3.0.713579) 64-bit (glnxa64)\n"
        "September 14, 2017\n\n \nFor online documentation, see http://www.mathworks.com/support\n"
        "For product information, visit www.mathworks.com.\n \n\x1b[?1l\x1b>"
    )


def test_get_matlab_version():
    from clinica.utils.check_dependency import get_software_version

    with mock.patch(
        "clinica.utils.check_dependency._get_matlab_start_session_message",
        wraps=matlab_version_mock,
    ) as matlab_mock:
        assert get_software_version("matlab") == Version("9.3.0.713579")
        matlab_mock.assert_called_once()


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
    from clinica.utils.check_dependency import _map_mcr_release_to_version_number  # noqa

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


def test_get_petpvc_version():
    from clinica.utils.check_dependency import get_software_version

    assert get_software_version("petpvc") == Version("0.0.0")


def test_check_software_version(mocker):
    from clinica.utils.check_dependency import _check_software_version  # noqa
    from clinica.utils.stream import LoggingLevel

    mocker.patch(
        "clinica.utils.check_dependency.get_software_version",
        return_value=Version("1.0.1"),
    )
    mocker.patch(
        "clinica.utils.check_dependency.get_software_min_version_supported",
        return_value=Version("1.1.0"),
    )

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="ants version is 1.0.1. We strongly recommend to have ants >=1.1.0.",
    ):
        _check_software_version(ThirdPartySoftware.ANTS, log_level=LoggingLevel.ERROR)
    _check_software_version(
        ThirdPartySoftware.ANTS,
        log_level=LoggingLevel.ERROR,
        specifier=SpecifierSet("==1.0.1"),
    )
    _check_software_version(
        ThirdPartySoftware.ANTS,
        log_level=LoggingLevel.ERROR,
        specifier=SpecifierSet("<2"),
    )
    _check_software_version(
        ThirdPartySoftware.ANTS,
        log_level=LoggingLevel.ERROR,
        specifier=SpecifierSet(">1.0"),
    )

    with pytest.raises(
        ClinicaMissingDependencyError,
        match="ants version is 1.0.1. We strongly recommend to have ants ==1.2.3.",
    ):
        _check_software_version(
            ThirdPartySoftware.ANTS,
            log_level=LoggingLevel.ERROR,
            specifier=SpecifierSet("==1.2.3"),
        )

    with pytest.warns(
        UserWarning,
        match="ants version is 1.0.1. We strongly recommend to have ants >=1.1.0.",
    ):
        _check_software_version(ThirdPartySoftware.ANTS, log_level=LoggingLevel.WARNING)

    with pytest.warns(
        UserWarning,
        match="ants version is 1.0.1. We strongly recommend to have ants <=0.23.",
    ):
        _check_software_version(
            ThirdPartySoftware.ANTS,
            log_level=LoggingLevel.WARNING,
            specifier=SpecifierSet("<=0.23"),
        )


def version_mock(
    software: Union[str, ThirdPartySoftware], value_for_mock: Optional[str] = None
) -> Version:
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if value_for_mock:
        return Version(value_for_mock)
    raise ClinicaMissingDependencyError("Not installed !")


software_dependency_test_cases = [
    ("2.3.4", "<1.2.3", False),
    ("2.3.4", ">1.2.3", True),
    ("0.1.0", "==0.1.0", True),
    ("0.1.0", "==0.0.1", False),
    ("0.0.0", "<0.0.0", False),
    ("0.0.0", "<=0.0.0", True),
    ("1.0.1", ">=0.129.1234", True),
    ("0.0.0", "", True),
    ("1.0.0", "", True),
    (None, ">=0.0.1", False),
]


@pytest.mark.parametrize(
    "installed_version,constraint,satisfied", software_dependency_test_cases
)
def test_software_dependency_string_constructor(
    installed_version: str, constraint: str, satisfied: bool
):
    from clinica.utils.check_dependency import SoftwareDependency

    with mock.patch(
        "clinica.utils.check_dependency.get_software_version",
        wraps=partial(version_mock, value_for_mock=installed_version),
    ) as version_mock_:
        dependency = SoftwareDependency.from_strings("ants", constraint)
        version_mock_.assert_called_once_with(ThirdPartySoftware.ANTS)
        assert dependency.name == ThirdPartySoftware.ANTS
        assert dependency.version_constraint == SpecifierSet(
            ">=0.0.0" if constraint == "" else constraint
        )
        if installed_version:
            assert dependency.installed_version == Version(installed_version)
        else:
            assert dependency.installed_version is None
        assert dependency.is_satisfied() is satisfied
        assert dependency.to_dict() == {
            "name": "ants",
            "version_constraint": ">=0.0.0" if constraint == "" else constraint,
            "installed_version": installed_version or "",
        }


@pytest.mark.parametrize(
    "installed_version,constraint,satisfied", software_dependency_test_cases
)
def test_software_dependency_dict_constructor(
    installed_version: str, constraint: str, satisfied: bool
):
    from clinica.utils.check_dependency import SoftwareDependency

    with mock.patch(
        "clinica.utils.check_dependency.get_software_version",
        wraps=partial(version_mock, value_for_mock=installed_version),
    ) as version_mock_:
        dependency = SoftwareDependency.from_dict(
            {"name": "ants", "version": constraint, "foo": "bar"}
        )
        version_mock_.assert_called_once_with(ThirdPartySoftware.ANTS)
        assert dependency.name == ThirdPartySoftware.ANTS
        assert dependency.version_constraint == SpecifierSet(
            ">=0.0.0" if constraint == "" else constraint
        )
        if installed_version:
            assert dependency.installed_version == Version(installed_version)
        else:
            assert dependency.installed_version is None
        assert dependency.is_satisfied() is satisfied
        assert dependency.to_dict() == {
            "name": "ants",
            "version_constraint": ">=0.0.0" if constraint == "" else constraint,
            "installed_version": installed_version or "",
        }
