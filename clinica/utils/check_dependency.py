"""This module contains utilities to check dependencies before running Clinica.

These functions can check binaries, software (e.g. FreeSurfer) or toolboxes (e.g. SPM).
"""
import functools
from typing import List, Optional, Tuple

from clinica.utils.exceptions import ClinicaMissingDependencyError


def is_binary_present(binary: str) -> bool:
    """Check if a binary is present.

    Parameters
    ----------
    binary : str
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
        subprocess.Popen([binary], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def check_environment_variable(environment_variable: str, software_name: str) -> str:
    """Check if the provided environment variable is set and returns its value.

    Parameters
    ----------
    environment_variable : str
        The name of the environment variable to check.

    software_name : str
        The name of the software related to the environment variable.

    Returns
    -------
    str :
        The value associated to the environment variable.

    Raises
    ------
    ClinicaMissingDependencyError
        If the variable is not set.
        If the variable associated value is not a directory.

    Examples
    --------
    >>> check_environment_variable("ANTSPATH", "ANTs")
    '/opt/ANTs/bin/'
    """
    import os

    content_var = os.environ.get(environment_variable, "")
    if not content_var:
        raise ClinicaMissingDependencyError(
            f"Clinica could not find {software_name} software: "
            f"the {environment_variable} variable is not set."
        )
    if not os.path.isdir(content_var):
        raise ClinicaMissingDependencyError(
            f"The {environment_variable} environment variable "
            f"you gave is not a folder (content: {content_var})."
        )
    return content_var


def _check_software(
    name: str,
    binaries: Optional[List[str]] = None,
    env: Optional[Tuple[str, str]] = None,
    complementary_info: Optional[str] = None,
) -> None:
    """Check if the software is available.

    Parameters
    ----------
    name : str
        Name of the software.

    binaries : list of str, optional
        List of associated binaries to check.
        If None, nothing is checked.

    env : (str, str), optional
        Tuple (environment variable, software name).
        This environment variable will be checked, meaning
        that it should be set and point to an existing folder.
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
    if env:
        check_environment_variable(*env)
    binaries = binaries or []
    complementary_info = complementary_info or ""
    for binary in binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"[Error] Clinica could not find {name} software: "
                f"the {binary} command is not present in your PATH "
                f"environment. {complementary_info}"
            )


check_dcm2niix = functools.partial(
    _check_software,
    name="dcm2niix",
    binaries=["dcm2niix"],
    complementary_info=(
        "This software can be downloaded and installed "
        "from https://github.com/rordenlab/dcm2niix."
    ),
)

check_ants = functools.partial(
    _check_software,
    name="ANTs",
    binaries=["N4BiasFieldCorrection", "antsRegistrationSyNQuick.sh"],
)

check_convert3d = functools.partial(
    _check_software,
    name="Convert3D",
    binaries=["c3d_affine_tool", "c3d"],
)

_check_freesurfer = functools.partial(
    _check_software,
    name="FreeSurfer",
    binaries=["mri_convert", "recon-all"],
    env=("FREESURFER_HOME", "FreeSurfer"),
    complementary_info=(
        "Do you have the line `source $FREESURFER_HOME/SetUpFreeSurfer.sh` "
        "in your configuration file?"
    ),
)

check_mrtrix = functools.partial(
    _check_software,
    name="MRtrix",
    binaries=["transformconvert", "mrtransform", "dwi2response", "tckgen"],
)

check_petpvc = functools.partial(
    _check_software,
    name="PETPVC",
    binaries=[
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
    ],
)

check_spm = functools.partial(
    _check_software,
    name="SPM",
    env=("SPM_HOME", "SPM"),
)

check_matlab = functools.partial(
    _check_software,
    name="Matlab",
    binaries=["matlab"],
)

_check_fsl = functools.partial(
    _check_software,
    name="FSL",
    binaries=["bet", "flirt", "fast", "first"],
    env=("FSLDIR", "FSL"),
)


def check_fsl() -> None:
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


def check_freesurfer() -> None:
    """Check FreeSurfer software."""
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
