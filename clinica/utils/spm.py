import warnings
from os import PathLike

"""This module contains SPM utilities."""

INDEX_TISSUE_MAP = {
    1: "graymatter",
    2: "whitematter",
    3: "csf",
    4: "bone",
    5: "softtissue",
    6: "background",
}


def check_spm_home():
    """Check and get SPM_HOME environment variable if present."""
    import os
    import platform

    from .check_dependency import check_environment_variable
    from .exceptions import ClinicaMissingDependencyError

    spm_home = check_environment_variable("SPM_HOME", "SPM")

    spm_standalone_home = os.environ.get("SPMSTANDALONE_HOME", "")
    if spm_standalone_home:
        if not os.path.isdir(spm_standalone_home):
            raise ClinicaMissingDependencyError(
                "The SPMSTANDALONE_HOME environment variable you "
                f"gave is not a folder (content: {spm_standalone_home})."
            )
        if platform.system() == "Darwin":
            spm_home = os.path.join(
                spm_standalone_home, "spm12.app", "Contents", "MacOS", "spm12_mcr"
            )
        else:
            spm_home = os.path.join(spm_standalone_home, "spm12_mcr")

    return spm_home


def get_tpm() -> PathLike:
    """Get Tissue Probability Map (TPM) from SPM.
    Returns
    -------
    PathLike :
        TPM.nii path from SPM
    """
    import os
    from glob import glob
    from os.path import join

    spm_home = os.getenv("SPM_HOME")

    if not spm_home:
        spm_home = os.getenv("SPMSTANDALONE_HOME")

    if not spm_home:
        raise RuntimeError(
            "Could not determine location of your SPM installation. Neither $SPM_HOME "
            "or $SPMSTANDALONE_HOME are present in your environment"
        )

    tpm_file_glob = glob(join(spm_home, "**/TPM.nii"), recursive=True)

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

    Returns
    -------
    bool :
        True if spm standalone was found and successfully configured.
        False otherwise.

    Raises
    ------
    FileNotFoundError :
        If the environment variables are set to non-existent folders.
    """
    import os
    import warnings

    from clinica.utils.stream import cprint

    if all(elem in os.environ.keys() for elem in ("SPMSTANDALONE_HOME", "MCR_HOME")):
        if os.path.isdir(os.path.expandvars("$SPMSTANDALONE_HOME")) and os.path.isdir(
            os.path.expandvars("$MCR_HOME")
        ):
            spm_standalone_home = os.getenv("SPMSTANDALONE_HOME")
            mcr_home = os.getenv("MCR_HOME")
            cprint(
                f"SPM standalone has been found at {spm_standalone_home}, "
                f"with an MCR at {mcr_home} and will be used in this pipeline"
            )
            matlab_command = _get_platform_dependant_matlab_command(
                spm_standalone_home, mcr_home
            )
            _configure_spm_nipype_interface(matlab_command)
            return True
        raise FileNotFoundError(
            "[Error] $SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder"
        )
    warnings.warn(
        "SPM standalone is not available on this system. "
        "The pipeline will try to use SPM and Matlab instead. "
        "If you want to rely on spm standalone, please make sure "
        "to set the following environment variables: "
        "$SPMSTANDALONE_HOME, $MCR_HOME, and $SPM_HOME."
    )
    return False


def _get_platform_dependant_matlab_command(
    spm_standalone_home: str, mcr_home: str
) -> str:
    import os
    import platform

    user_system = platform.system().lower()
    if user_system.startswith("darwin"):
        return f"cd {spm_standalone_home} && ./run_spm12.sh {mcr_home} script"
    if user_system.startswith("linux"):
        return f"{os.path.join(spm_standalone_home, 'run_spm12.sh')} {mcr_home} script"
    raise SystemError(
        f"Clinica only support macOS and Linux. Your system is {user_system}."
    )


def _configure_spm_nipype_interface(matlab_command: str) -> str:
    from nipype.interfaces import spm

    from clinica.utils.stream import cprint

    spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_command, use_mcr=True)
    cprint(f"Using SPM standalone version {spm.SPMCommand().version}")
