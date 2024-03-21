"""This module contains utilities to check dependencies before running Clinica.

These functions can check binaries, software (e.g. FreeSurfer) or toolboxes (e.g. SPM).
"""
import functools
import os
from enum import Enum
from pathlib import Path
from typing import Optional, Tuple, Union

from clinica.utils.exceptions import ClinicaMissingDependencyError

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


def _check_fsl_above_version_five() -> None:
    """Check FSL software."""
    import nipype.interfaces.fsl as fsl

    from clinica.utils.stream import cprint

    _check_fsl()
    try:
        if fsl.Info.version().split(".") < ["5", "0", "5"]:
            raise ClinicaMissingDependencyError(
                "FSL version must be greater than 5.0.5"
            )
    except Exception as e:
        cprint(msg=str(e), lvl="error")


def _check_freesurfer_above_version_six() -> None:
    """Check FreeSurfer software >= 6.0.0."""
    import nipype.interfaces.freesurfer as freesurfer

    from clinica.utils.stream import cprint

    _check_freesurfer()
    try:
        if freesurfer.Info().version().split("-")[3].split(".") < ["6", "0", "0"]:
            raise ClinicaMissingDependencyError(
                "FreeSurfer version must be greater than 6.0.0"
            )
    except Exception as e:
        cprint(msg=str(e), lvl="error")


def check_software(software: Union[str, ThirdPartySoftware]):
    software = ThirdPartySoftware(software)
    if software == ThirdPartySoftware.ANTS:
        return _check_ants()
    if software == ThirdPartySoftware.FSL:
        return _check_fsl_above_version_five()
    if software == ThirdPartySoftware.FREESURFER:
        return _check_freesurfer_above_version_six()
    if (
        software == ThirdPartySoftware.SPM
        or software == ThirdPartySoftware.SPMSTANDALONE
        or software == ThirdPartySoftware.MCR
    ):
        return _check_spm()
    if software == ThirdPartySoftware.MATLAB:
        return _check_matlab()
    if software == ThirdPartySoftware.DCM2NIIX:
        return _check_dcm2niix()
    if software == ThirdPartySoftware.PETPVC:
        return _check_petpvc()
    if software == ThirdPartySoftware.MRTRIX:
        return _check_mrtrix()
    if software == ThirdPartySoftware.CONVERT3D:
        return _check_convert3d()
