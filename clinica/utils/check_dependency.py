"""This module contains utilities to check dependencies before running Clinica.

These functions can check binaries, software (e.g. FreeSurfer) or toolboxes (e.g. SPM).
"""
import functools
import os
from enum import Enum
from functools import partial
from pathlib import Path
from typing import Optional, Tuple, Union

from packaging.specifiers import SpecifierSet
from packaging.version import Version

from clinica.utils.exceptions import ClinicaMissingDependencyError
from clinica.utils.stream import log_and_raise, log_and_warn

__all__ = [
    "ThirdPartySoftware",
    "SoftwareEnvironmentVariable",
    "get_fsl_home",
    "get_freesurfer_home",
    "get_mcr_home",
    "get_spm_standalone_home",
    "get_spm_home",
    "is_binary_present",
    "check_binary",
    "check_environment_variable",
    "check_software",
    "get_software_min_version_supported",
    "get_software_version",
]


class ThirdPartySoftware(str, Enum):
    """Possible third party software that clinica depends on."""

    ANTS = "ants"
    CONVERT3D = "convert3d"
    DCM2NIIX = "dcm2niix"
    FREESURFER = "freesurfer"
    FSL = "fsl"
    MATLAB = "matlab"
    MCR = "MCR"
    MRTRIX = "mrtrix"
    PETPVC = "petpvc"
    SPM = "spm"
    SPMSTANDALONE = "spm standalone"


class SoftwareEnvironmentVariable:
    """Represents an environment variable name linked to a specific software.

    Attributes
    ----------
    name : str
        The name of the environment variable.

    software : ThirdPartySoftware
        The name of the software related to the environment variable.
    """

    def __init__(self, name: str, software: Union[str, ThirdPartySoftware]):
        self.name = name
        self.software = ThirdPartySoftware(software)


def get_fsl_home() -> Path:
    """Return the path to the home directory of FSL."""
    return check_environment_variable(
        SoftwareEnvironmentVariable("FSLDIR", ThirdPartySoftware.FSL)
    )


def get_freesurfer_home() -> Path:
    """Return the path to the home directory of Freesurfer."""
    return check_environment_variable(
        SoftwareEnvironmentVariable("FREESURFER_HOME", ThirdPartySoftware.FREESURFER)
    )


def get_mcr_home() -> Path:
    """Return the path to the home directory of Matlab MCR."""
    return check_environment_variable(
        SoftwareEnvironmentVariable("MCR_HOME", ThirdPartySoftware.MCR)
    )


def get_spm_standalone_home() -> Path:
    """Return the path to the home directory of SPM standalone"""
    return check_environment_variable(
        SoftwareEnvironmentVariable(
            "SPMSTANDALONE_HOME", ThirdPartySoftware.SPMSTANDALONE
        )
    )


def get_spm_home() -> Path:
    """Return the path to the home directory of SPM."""
    return check_environment_variable(
        SoftwareEnvironmentVariable("SPM_HOME", ThirdPartySoftware.SPM)
    )


def is_binary_present(name: str) -> bool:
    """Check if a binary is present.

    Parameters
    ----------
    name : str
        The name of the program.

    Returns
    -------
    True if the binary is present, False otherwise.

    Warnings
    --------
    Do not use this function with a binary GUI, it will open the GUI.

    References
    ----------
    Taken from:
    https://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script

    Examples
    --------
    >>> is_binary_present("ls")
    True
    >>> is_binary_present("foo")
    False
    """
    import os
    import subprocess

    from clinica.compat import errno

    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def check_binary(name: str):
    if not is_binary_present(name):
        raise ClinicaMissingDependencyError(
            f"Command {name} is unknown on your system. "
            "Please verify that you have correctly installed the corresponding software."
        )


def check_environment_variable(variable: SoftwareEnvironmentVariable) -> Path:
    """Check if the provided environment variable is set and returns its value.

    Parameters
    ----------
    variable : SoftwareEnvironmentVariable
        The software environment variable to check.

    Returns
    -------
    Path :
        The path value associated to the environment variable.

    Raises
    ------
    ClinicaMissingDependencyError :
        If the variable is not set.

    ClinicaEnvironmentVariableError :
        If the variable associated value is not a directory.

    Examples
    --------
    >>> check_environment_variable(SoftwareEnvironmentVariable("ANTSPATH", ThirdPartySoftware.ANTS))
    '/opt/ANTs/bin/'
    """
    from .exceptions import ClinicaEnvironmentVariableError

    if not (content_var := os.environ.get(variable.name, "")):
        raise ClinicaMissingDependencyError(
            f"Clinica could not find {variable.software.value} software: "
            f"the {variable.name} variable is not set."
        )
    if not os.path.isdir(content_var):
        raise ClinicaEnvironmentVariableError(
            f"The {variable.name} environment variable "
            f"you gave is not a folder (content: {content_var})."
        )
    return Path(content_var)


def _check_software(
    name: ThirdPartySoftware,
    binaries: Optional[Tuple[str]] = None,
    environment_variables: Optional[Tuple[SoftwareEnvironmentVariable]] = None,
    complementary_info: Optional[str] = None,
) -> None:
    """Check if the software is available.

    Parameters
    ----------
    name : ThirdPartySoftware
        Name of the software.

    binaries : tuple of str, optional
        Tuple of associated binaries to check.
        If None, nothing is checked.

    environment_variables : tuple of SoftwareEnvironmentVariable, optional
        These environment variables will be checked, meaning
        that they should be set and point to existing folders.
        If None, nothing is checked.

    complementary_info : str, optional
        Information specific to the software that should be
        added at the end of the error message when checking
        for binaries.

    Raises
    ------
    ClinicaMissingDependencyError
        If the software checks fail.
    """
    if environment_variables:
        for variable in environment_variables:
            check_environment_variable(variable)
    binaries = binaries or ()
    complementary_info = complementary_info or ""
    for binary in binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"[Error] Clinica could not find {name.value} software: "
                f"the {binary} command is not present in your PATH "
                f"environment. {complementary_info}"
            )


_check_dcm2niix = functools.partial(
    _check_software,
    name=ThirdPartySoftware.DCM2NIIX,
    binaries=("dcm2niix",),
    complementary_info=(
        "This software can be downloaded and installed "
        "from https://github.com/rordenlab/dcm2niix."
    ),
)

_check_ants = functools.partial(
    _check_software,
    name=ThirdPartySoftware.ANTS,
    binaries=("N4BiasFieldCorrection", "antsRegistrationSyNQuick.sh"),
)

_check_convert3d = functools.partial(
    _check_software,
    name=ThirdPartySoftware.CONVERT3D,
    binaries=("c3d_affine_tool", "c3d"),
)

_check_freesurfer = functools.partial(
    _check_software,
    name=ThirdPartySoftware.FREESURFER,
    binaries=("mri_convert", "recon-all"),
    environment_variables=(
        SoftwareEnvironmentVariable("FREESURFER_HOME", ThirdPartySoftware.FREESURFER),
    ),
    complementary_info=(
        "Do you have the line `source $FREESURFER_HOME/SetUpFreeSurfer.sh` "
        "in your configuration file?"
    ),
)

_check_mrtrix = functools.partial(
    _check_software,
    name=ThirdPartySoftware.MRTRIX,
    binaries=("transformconvert", "mrtransform", "dwi2response", "tckgen"),
)

_check_petpvc = functools.partial(
    _check_software,
    name=ThirdPartySoftware.PETPVC,
    binaries=(
        "petpvc",
        "pvc_diy",
        "pvc_gtm",
        "pvc_iy",
        "pvc_labbe",
        "pvc_make4d",
        "pvc_mg",
        "pvc_mtc",
        "pvc_rbv",
        "pvc_relabel",
        "pvc_rl",
        "pvc_simulate",
        "pvc_stc",
        "pvc_vc",
    ),
)


def _check_spm():
    """Check that SPM is installed, either regular with Matlab or as a standalone."""
    try:
        _check_spm_standalone()
    except ClinicaMissingDependencyError as e1:
        try:
            _check_spm_alone()
        except ClinicaMissingDependencyError as e2:
            raise ClinicaMissingDependencyError(
                "Clinica could not find the SPM software (regular or standalone).\n"
                "Please make sure you have installed SPM in one of these two ways "
                "and have set the required environment variables.\n"
                f"Full list of errors: \n- {e1}\n- {e2}"
            )


_check_spm_standalone = functools.partial(
    _check_software,
    name=ThirdPartySoftware.SPMSTANDALONE,
    environment_variables=(
        SoftwareEnvironmentVariable(
            "SPMSTANDALONE_HOME", ThirdPartySoftware.SPMSTANDALONE
        ),
        SoftwareEnvironmentVariable("MCR_HOME", ThirdPartySoftware.MCR),
    ),
)

_check_spm_alone = functools.partial(
    _check_software,
    name=ThirdPartySoftware.SPM,
    environment_variables=(
        SoftwareEnvironmentVariable("SPM_HOME", ThirdPartySoftware.SPM),
    ),
    binaries=("matlab",),
)

_check_matlab = functools.partial(
    _check_software,
    name=ThirdPartySoftware.MATLAB,
    binaries=("matlab",),
)

_check_fsl = functools.partial(
    _check_software,
    name=ThirdPartySoftware.FSL,
    binaries=("bet", "flirt", "fast", "first"),
    environment_variables=(
        SoftwareEnvironmentVariable("FSLDIR", ThirdPartySoftware.FSL),
    ),
)


def get_software_min_version_supported(
    software: Union[str, ThirdPartySoftware],
) -> Version:
    """Return the minimum version of the provided third-party software required by Clinica.

    Parameters
    ----------
    software : str or ThirdPartySoftware
        One of the third-party software of Clinica.

    Returns
    -------
    Version :
        The minimum version number of the software required by Clinica.

    Examples
    --------
    >>> from clinica.utils.check_dependency import get_software_min_version_supported
    >>> get_software_min_version_supported("ants")
    <Version('2.5.0')>
    """
    software = ThirdPartySoftware(software)
    if software == ThirdPartySoftware.FREESURFER:
        return Version("6.0.0")
    if software == ThirdPartySoftware.FSL:
        return Version("5.0.5")
    if software == ThirdPartySoftware.ANTS:
        return Version("2.5.0")
    if software == ThirdPartySoftware.DCM2NIIX:
        return Version("1.0.20240202")
    if software == ThirdPartySoftware.MRTRIX:
        return Version("3.0.3")
    if software == ThirdPartySoftware.CONVERT3D:
        return Version("1.0.0")
    if software == ThirdPartySoftware.MATLAB:
        return Version("9.2.0.556344")
    if software == ThirdPartySoftware.SPM:
        return Version("12.7219")
    if software == ThirdPartySoftware.MCR:
        return Version("9.0.1")
    if software == ThirdPartySoftware.SPMSTANDALONE:
        return Version("12.7219")
    if software == ThirdPartySoftware.PETPVC:
        return Version("0.0.0")


def get_software_version(software: Union[str, ThirdPartySoftware]) -> Version:
    """Return the version of the provided third-party software.

    Parameters
    ----------
    software : str or ThirdPartySoftware
        One of the third-party software of Clinica.

    Returns
    -------
    Version :
        The version number of the installed software.

    Notes
    -----
    This function assumes the software are correctly installed.
    It doesn't run any check and directly try to infer the version number by calling an
    underlying executable.

    Examples
    --------
    >>> from clinica.utils.check_dependency import get_software_version
    >>> get_software_version("freesurfer")
    <Version('7.2.0')>
    """
    software = ThirdPartySoftware(software)
    if software == ThirdPartySoftware.FREESURFER:
        return _get_freesurfer_version()
    if software == ThirdPartySoftware.FSL:
        return _get_fsl_version()
    if software == ThirdPartySoftware.ANTS:
        return _get_ants_version()
    if software == ThirdPartySoftware.DCM2NIIX:
        return _get_dcm2niix_version()
    if software == ThirdPartySoftware.MRTRIX:
        return _get_mrtrix_version()
    if software == ThirdPartySoftware.CONVERT3D:
        return _get_convert3d_version()
    if software == ThirdPartySoftware.MATLAB:
        return _get_matlab_version()
    if software == ThirdPartySoftware.SPM:
        return _get_spm_version()
    if software == ThirdPartySoftware.MCR:
        return _get_mcr_version()
    if software == ThirdPartySoftware.SPMSTANDALONE:
        return _get_spm_standalone_version()
    if software == ThirdPartySoftware.PETPVC:
        return Version("0.0.0")


def _get_freesurfer_version() -> Version:
    from nipype.interfaces import freesurfer

    return Version(str(freesurfer.Info.looseversion()))


def _get_spm_version() -> Version:
    from nipype.interfaces import spm

    return Version(spm.SPMCommand().version)


def _get_spm_standalone_version() -> Version:
    import os
    from pathlib import Path

    from nipype.interfaces import spm

    spm_path = Path(os.environ["SPM_HOME"])
    matlab_cmd = f"{spm_path / 'run_spm12.sh'} {os.environ['MCR_HOME']} script"
    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
    return Version(spm.SPMCommand().version)


def _get_fsl_version() -> Version:
    import re

    from nipype.interfaces import fsl

    from clinica.utils.stream import cprint

    raw_output = str(fsl.Info.version())
    try:
        return Version(re.search(r"\s*([\d.]+)", raw_output).group(1))
    except Exception as e:
        cprint(msg=str(e), lvl="error")


def _get_software_version_from_command_line(
    executable: str, prepend_with_v: bool = True, two_dashes: bool = True
) -> Version:
    import re

    from clinica.utils.stream import log_and_raise

    try:
        return Version(
            re.search(
                rf"{'v' if prepend_with_v else ''}\s*([\d.]+)",
                _run_command(executable, two_dashes=two_dashes),
            )
            .group(1)
            .strip(".")
        )
    except Exception as e:
        log_and_raise(str(e), RuntimeError)


def _run_command(executable: str, two_dashes: bool = True) -> str:
    import subprocess

    return subprocess.run(
        [executable, f"-{'-' if two_dashes else ''}version"], stdout=subprocess.PIPE
    ).stdout.decode("utf-8")


_get_ants_version = partial(
    _get_software_version_from_command_line,
    executable="antsRegistration",
    prepend_with_v=False,
)
_get_dcm2niix_version = partial(
    _get_software_version_from_command_line, executable="dcm2niix"
)
_get_mrtrix_version = partial(
    _get_software_version_from_command_line,
    executable="mrtransform",
    prepend_with_v=False,
)
_get_convert3d_version = partial(
    _get_software_version_from_command_line,
    executable="c3d",
    prepend_with_v=False,
    two_dashes=False,
)


def _get_matlab_version() -> Version:
    import re

    from clinica.utils.stream import log_and_raise

    try:
        return Version(
            re.search(r"\(\s*([\d.]+)\)", _get_matlab_start_session_message())
            .group(1)
            .strip(".")
        )
    except Exception as e:
        log_and_raise(str(e), RuntimeError)


def _get_matlab_start_session_message() -> str:
    """Start Matlab, get the message displayed at the beginning of the session, and quit."""
    import subprocess

    return subprocess.run(
        ["matlab", "-r", "quit", "-nojvm"], stdout=subprocess.PIPE
    ).stdout.decode("utf-8")


def _get_mcr_version() -> Version:
    import os

    mcr_home_path = Path(os.environ.get("MCR_HOME"))
    raw_version = mcr_home_path.parent.name
    return _map_mcr_release_to_version_number(raw_version)


def _map_mcr_release_to_version_number(mcr_release: str) -> Version:
    """Map the MCR release to a proper version number.

    If the found version is older than the minimum supported version, we don't bother mapping it
    to its version number, 0.0.0 is returned.
    If the release is >= 2024a, we use the fact that Matlab is now using proper calendar versioning.
    If the release is in-between, we use a hardcoded mapping.
    """
    mcr_versions_mapping = {
        "2023b": Version("23.2"),
        "2023a": Version("9.14"),
        "2022b": Version("9.13"),
        "2022a": Version("9.12"),
        "2021b": Version("9.11"),
        "2021a": Version("9.10"),
        "2020b": Version("9.9"),
        "2020a": Version("9.8"),
        "2019b": Version("9.7"),
        "2019a": Version("9.6"),
        "2018b": Version("9.5"),
        "2018a": Version("9.4"),
        "2017b": Version("9.3"),
        "2017a": Version("9.2"),
        "2016b": Version("9.1"),
        "2016a": Version("9.0.1"),
    }
    if int(mcr_release[:4]) >= 2024:
        return Version(f"{mcr_release[2:4]}.{'1' if mcr_release[-1] == 'a' else '2'}")
    return mcr_versions_mapping.get(mcr_release, Version("0.0.0"))


class SeverityLevel(str, Enum):
    ERROR = "error"
    WARNING = "warning"


def _check_software_version(
    software: ThirdPartySoftware,
    severity: Optional[SeverityLevel] = None,
    specifier: Optional[SpecifierSet] = None,
):
    """Check that the installed version of the software is satisfying the constraints imposed by Clinica."""
    from clinica.utils.stream import cprint

    severity = severity or SeverityLevel.WARNING
    if specifier is None:
        specifier = SpecifierSet(f">={get_software_min_version_supported(software)}")
    if (installed_version := get_software_version(software)) not in specifier:
        (log_and_raise if severity == SeverityLevel.ERROR else log_and_warn)(
            f"{software.value} version is {installed_version}. We strongly recommend to have {software.value} {specifier}.",
            (
                ClinicaMissingDependencyError
                if severity == SeverityLevel.ERROR
                else UserWarning
            ),
        )
    cprint(
        f"Found installation of {software.value} with version {installed_version}, satisfying {specifier}.",
        lvl="info",
    )


def check_software(
    software: Union[str, ThirdPartySoftware],
    specifier: Optional[Union[str, SpecifierSet]] = None,
):
    """Run some checks on the given software.

    These checks are of two types:
        - checks to verify that the executable is present in the PATH.
          Also check configurations made through environment variables.
        - checks on the installed version. The installed version has to
          satisfy the constraints imposed by Clinica.

    Parameters
    ----------
    software : str or ThirdPartySoftware
        One of the third-party software of Clinica.

    specifier : str or SpecifierSet, optional
        A version constraint for the software (i.e. '>=2.0.1', '<1.0.0', or '==0.10.9').
        If not provided, the specifier will be '>= minimum_version' where 'minimum_version'
        is the minimum version number of the software required by Clinica.

    Raises
    ------
    ClinicaMissingDependencyError :
        If an issue is found with the installation of the provided software.
        In some cases where the software is correctly installed, but the version is considered
        a bit old compared to the one recommended by Clinica, a warning could be given instead of
        an error. In this situation, it is strongly recommended to upgrade the dependency even
        though Clinica does not raise.

    Examples
    --------
    >>> from clinica.utils.check_dependency import check_software
    >>> check_software("dcm2niix")
    >>> check_software("ants", ">=2.5")
    """
    software = ThirdPartySoftware(software)
    if specifier:
        specifier = SpecifierSet(specifier)
    severity = SeverityLevel.WARNING
    if software == ThirdPartySoftware.ANTS:
        _check_ants()
    if software == ThirdPartySoftware.FSL:
        _check_fsl()
        severity = SeverityLevel.ERROR
    if software == ThirdPartySoftware.FREESURFER:
        _check_freesurfer()
        severity = SeverityLevel.ERROR
    if (
        software == ThirdPartySoftware.SPM
        or software == ThirdPartySoftware.SPMSTANDALONE
        or software == ThirdPartySoftware.MCR
    ):
        _check_spm()
        severity = SeverityLevel.ERROR
    if software == ThirdPartySoftware.MATLAB:
        _check_matlab()
    if software == ThirdPartySoftware.DCM2NIIX:
        _check_dcm2niix()
    if software == ThirdPartySoftware.PETPVC:
        _check_petpvc()
    if software == ThirdPartySoftware.MRTRIX:
        _check_mrtrix()
    if software == ThirdPartySoftware.CONVERT3D:
        _check_convert3d()
    _check_software_version(software, severity, specifier)
