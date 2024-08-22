"""This module contains SPM utilities."""
import warnings
from os import PathLike
from pathlib import Path

__all__ = [
    "INDEX_TISSUE_MAP",
    "get_tpm",
    "use_spm_standalone_if_available",
]

INDEX_TISSUE_MAP = {
    1: "graymatter",
    2: "whitematter",
    3: "csf",
    4: "bone",
    5: "softtissue",
    6: "background",
}


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
    from clinica.utils.stream import cprint

    from .check_dependency import get_mcr_home, get_spm_home, get_spm_standalone_home
    from .exceptions import ClinicaMissingDependencyError

    try:
        spm_standalone_home = get_spm_standalone_home()
        mcr_home = get_mcr_home()
        cprint(
            f"SPM standalone has been found at {spm_standalone_home}, "
            f"with an MCR at {mcr_home} and will be used in this pipeline"
        )
        matlab_command = _get_platform_dependant_matlab_command(
            spm_standalone_home, mcr_home
        )
        _configure_spm_nipype_interface(matlab_command)
        return True
    except ClinicaMissingDependencyError:
        warnings.warn(
            "SPM standalone is not available on this system. "
            "The pipeline will try to use SPM and Matlab instead. "
            "If you want to rely on spm standalone, please make sure "
            "to set the following environment variables: "
            "$SPMSTANDALONE_HOME, and $MCR_HOME"
        )
        import nipype.interfaces.matlab as mlab

        cprint(f"Setting SPM path to {get_spm_home()}", lvl="info")
        mlab.MatlabCommand.set_default_paths(f"{get_spm_home()}")
        return False


def _get_platform_dependant_matlab_command(
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


def _configure_spm_nipype_interface(matlab_command: str):
    from nipype.interfaces import spm

    from clinica.utils.stream import cprint

    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_command, use_mcr=True)
    cprint(f"Using SPM standalone version {spm.SPMCommand().version}")
