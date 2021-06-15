# coding: utf8


"""This module contains utilities to check dependencies before running Clinica.

These functions can check binaries, software (e.g. FreeSurfer) or toolboxes (e.g. SPM).
"""


def is_binary_present(binary):
    """Check if a binary is present.

    This function checks if the program is present. Do not use this function
    with a binary GUI, it will open the GUI.

    Taken from:
    https://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script

    Args:
        binary (str): Name of the program.

    Returns:
        True if the binary is present, False otherwise.
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


def check_environment_variable(environment_variable, software_name):
    import os

    from colorama import Fore

    from .exceptions import ClinicaMissingDependencyError

    content_var = os.environ.get(environment_variable, "")
    if not content_var:
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}\n[Error] Clinica could not find {software_name} software: "
            f"the {environment_variable} variable is not set.{Fore.RESET}"
        )
    if not os.path.isdir(content_var):
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}\n[Error] The {environment_variable} environment variable "
            f"you gave is not a folder (content: {content_var}).{Fore.RESET}"
        )
    return content_var


def check_software_requirements(current_version, version_requirements, software_name):
    from string import punctuation

    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    comparison_operator = "".join(
        [c for c in version_requirements if c in punctuation.replace(".", "")]
    )
    required_version = version_requirements.replace(comparison_operator, "")

    satisfy_version = eval(
        f"LooseVersion('{current_version}') {comparison_operator} LooseVersion('{required_version}')"
    )

    if not satisfy_version:
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}\n[Error] Your {software_name} version ({current_version}) "
            f"does not satisfy version requirements ({version_requirements}).{Fore.RESET}"
        )


def check_dcm2nii():
    """Check dcm2nii software."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present("dcm2nii"):
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}\n[Error] Clinica could not find dcm2nii tool from MRIcron in your PATH environment: "
            f"this can be downloaded from https://www.nitrc.org/frs/?group_id=152 (choose the 2016 version).{Fore.RESET}"
        )


def check_dcm2niix():
    """Check dcm2niix software."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present("dcm2niix"):
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}\n[Error] Clinica could not find dcm2niix software in your PATH environment: "
            f"this can be downloaded or installed from https://github.com/rordenlab/dcm2niix.{Fore.RESET}"
        )


def check_ants(version_requirements=None):
    """Check ANTs software."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable("ANTSPATH", "ANTs")

    list_binaries = ["N4BiasFieldCorrection", "antsRegistrationSyNQuick.sh"]
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}\n[Error] Clinica could not find ANTs software: "
                f"the {binary} command is not present in your PATH environment.{Fore.RESET}"
            )


def check_freesurfer(version_requirements=None):
    """Check FreeSurfer software."""
    import nipype.interfaces.freesurfer as freesurfer
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable("FREESURFER_HOME", "FreeSurfer")

    list_binaries = ["mri_convert", "recon-all"]
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}\n[Error] Clinica could not find FreeSurfer software: "
                f"the {binary} command is not present in your PATH environment: "
                f"did you have the line source ${{FREESURFER_HOME}}/SetUpFreeSurfer.sh "
                f"in your configuration file?{Fore.RESET}"
            )


def check_fsl(version_requirements=None):
    """Check FSL software."""
    import nipype.interfaces.fsl as fsl
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError
    from clinica.utils.stream import cprint

    check_environment_variable("FSLDIR", "FSL")

    try:
        if fsl.Info.version().split(".") < ["5", "0", "5"]:
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}FSL version must be greater than 5.0.5{Fore.RESET}"
            )
    except Exception as e:
        cprint(str(e))

    list_binaries = ["bet", "flirt", "fast", "first"]
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}\n[Error] Clinica could not find FSL software: "
                f"the {binary} command is not present in your PATH environment.{Fore.RESET}"
            )


def check_mrtrix(version_requirements=None):
    """Check MRtrix software."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable("MRTRIX_HOME", "MRtrix")

    list_binaries = ["transformconvert", "mrtransform", "dwi2response", "tckgen"]
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}\n[Error] Clinica could not find MRtrix software: "
                f"the {binary} command is not present in your PATH environment.{Fore.RESET}"
            )


def check_petpvc(version_requirements=None):
    """Check PETPVC software."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    list_binaries = [
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
    ]
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                f"{Fore.RED}\n[Error] Clinica could not find PETPVC software: the "
                f"{binary} command is not present in your PATH environment.{Fore.RESET}"
            )


def check_spm(version_requirements=None):
    """Check SPM software."""
    check_environment_variable("SPM_HOME", "SPM")

    # list_binaries = ['matlab']
    # for binary in list_binaries:
    #     if not is_binary_present(binary):
    #         raise RuntimeError(
    #             f"{binary} from SPM Software is not present in your PATH environment.
    #           )


def check_matlab():
    """Check Matlab toolbox."""
    from colorama import Fore

    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present("matlab"):
        raise ClinicaMissingDependencyError(
            f"{Fore.RED}Matlab was not found in PATH environment. Did you add it?{Fore.RESET}"
        )
