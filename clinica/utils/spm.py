"""This module contains SPM utilities."""
import warnings
from enum import Enum
from os import PathLike
from pathlib import Path
from typing import Union

__all__ = [
    "SPMTissue",
    "get_spm_tissue_from_index",
    "get_tpm",
    "use_spm_standalone_if_available",
    "configure_nipype_interface_to_work_with_spm",
    "configure_nipype_interface_to_work_with_spm_standalone",
]


class SPMTissue(str, Enum):
    GRAY_MATTER = "graymatter"
    WHITE_MATTER = "whitematter"
    CSF = "csf"
    BONE = "bone"
    SOFT_TISSUE = "softtissue"
    BACKGROUND = "background"


def get_spm_tissue_from_index(index: Union[int, SPMTissue]) -> SPMTissue:
    if isinstance(index, SPMTissue):
        return index
    if index == 1:
        return SPMTissue.GRAY_MATTER
    if index == 2:
        return SPMTissue.WHITE_MATTER
    if index == 3:
        return SPMTissue.CSF
    if index == 4:
        return SPMTissue.BONE
    if index == 5:
        return SPMTissue.SOFT_TISSUE
    if index == 6:
        return SPMTissue.BACKGROUND
    raise ValueError(f"No SPM tissue matching index {index}.")


def get_tpm() -> PathLike:
    """Get Tissue Probability Map (TPM) from SPM.
    Returns
    -------
    PathLike :
        TPM.nii path from SPM
    """
    from glob import glob

    from .check_dependency import get_spm_home

    spm_home = get_spm_home()
    tpm_file_glob = glob(str(spm_home / "**/TPM.nii"), recursive=True)
    if len(tpm_file_glob) == 0:
        raise RuntimeError(f"No file found for TPM.nii in your $SPM_HOME in {spm_home}")
    if len(tpm_file_glob) > 1:
        error_str = f"Multiple files found for TPM.nii in your SPM_HOME {spm_home}:"
        for file in tpm_file_glob:
            error_str += "\n\t" + file
        raise RuntimeError(error_str)

    return tpm_file_glob[0]


def use_spm_standalone_if_available() -> bool:
    """Use SPM Standalone with MATLAB Common Runtime if it can be used on the user system.

    If there is something wrong with either SPM Standalone or the MCR environment variables
    configuration, a warning is given to the user, and False is returned.
    It is thus possible to use this function to try using the standalone version of SPM12
    while potentially defaulting to the standard SPM12 install.

    Returns
    -------
    bool :
        True if spm standalone was found and successfully configured.
        False otherwise.

    Raises
    ------
    ClinicaEnvironmentVariableError :
        If the environment variables are set to non-existent folders.
    """
    from .exceptions import ClinicaMissingDependencyError
    from .stream import log_and_warn

    try:
        configure_nipype_interface_to_work_with_spm_standalone()
        return True
    except ClinicaMissingDependencyError:
        log_and_warn(
            (
                "SPM standalone is not available on this system. "
                "The pipeline will try to use SPM and Matlab instead. "
                "If you want to rely on spm standalone, please make sure "
                "to set the following environment variables: "
                "$SPMSTANDALONE_HOME, and $MCR_HOME"
            ),
            UserWarning,
        )
        configure_nipype_interface_to_work_with_spm()
        return False


def configure_nipype_interface_to_work_with_spm() -> None:
    import nipype.interfaces.matlab as mlab

    from clinica.utils.stream import cprint

    from .check_dependency import get_spm_home

    cprint(f"Setting SPM path to {get_spm_home()}", lvl="info")
    mlab.MatlabCommand.set_default_paths(f"{get_spm_home()}")


def _get_platform_dependant_matlab_command_for_spm_standalone(
    spm_standalone_home: Path, mcr_home: Path
) -> str:
    import platform

    user_system = platform.system().lower()
    if user_system.startswith("darwin"):
        return f"cd {spm_standalone_home} && ./run_spm12.sh {mcr_home} script"
    if user_system.startswith("linux"):
        return f"{spm_standalone_home / 'run_spm12.sh'} {mcr_home} script"
    raise SystemError(
        f"Clinica only support macOS and Linux. Your system is {user_system}."
    )


def configure_nipype_interface_to_work_with_spm_standalone() -> None:
    from nipype.interfaces import spm

    from clinica.utils.stream import cprint

    from .check_dependency import get_mcr_home, get_spm_standalone_home

    spm_standalone_home = get_spm_standalone_home()
    mcr_home = get_mcr_home()
    cprint(
        f"SPM standalone has been found at {spm_standalone_home}, "
        f"with an MCR at {mcr_home} and will be used in this pipeline"
    )
    spm.SPMCommand.set_mlab_paths(
        matlab_cmd=_get_platform_dependant_matlab_command_for_spm_standalone(
            spm_standalone_home, mcr_home
        ),
        use_mcr=True,
    )
    cprint(f"Using SPM standalone version {spm.SPMCommand().version}", lvl="info")
