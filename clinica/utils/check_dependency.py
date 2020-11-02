# coding: utf8


"""
This module contains utilities to check dependencies before running Clinica.

These functions can check binaries, software (e.g. FreeSurfer) or toolboxes (e.g. SPM).
"""


def is_binary_present(binary):
    """
    Check if a binary is present.

    This function checks if the program is present. Do not use this function
    with a binary GUI, it will open the GUI.

    Taken from:
    https://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script

    Args:
        binary (str): Name of the program.

    Returns:
        True if the binary is present, False otherwise.
    """
    import subprocess
    import os
    import sys
    if sys.version_info[0] >= 3 and sys.version_info[1] >= 7:
        import errno
    else:
        errno = os.errno

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

    content_var = os.environ.get(environment_variable, '')
    if not content_var:
        raise ClinicaMissingDependencyError(
            '%s\n[Error] Clinica could not find %s software: the %s variable is not set.%s'
            % (Fore.RED, software_name, environment_variable, Fore.RESET))
    if not os.path.isdir(content_var):
        raise ClinicaMissingDependencyError(
            '%s\n[Error] The %s environment variable you gave is not a folder (content: %s).%s'
            % (Fore.RED, environment_variable, content_var, Fore.RESET))
    return content_var


def check_software_requirements(current_version, version_requirements, software_name):
    from distutils.version import LooseVersion
    from string import punctuation
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    comparison_operator = ''.join([c for c in version_requirements if c in punctuation.replace('.', '')])
    required_version = version_requirements.replace(comparison_operator, '')

    satisfy_version = eval('LooseVersion(\'%s\') %s LooseVersion(\'%s\')' %
                           (current_version, comparison_operator, required_version))

    if not satisfy_version:
        raise ClinicaMissingDependencyError(
            '%s\n[Error] Your %s version (%s) does not version requirements (%s).%s' %
            (Fore.RED, software_name, current_version, version_requirements, Fore.RESET))


def check_dcm2nii():
    """Check dcm2nii software."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present('dcm2nii'):
        raise ClinicaMissingDependencyError(
            '%s\n[Error] Clinica could not find dcm2nii tool from MRIcron in your PATH environment: '
            'this can be downloaded from https://www.nitrc.org/frs/?group_id=152 (choose the 2016 version).%s' %
            (Fore.RED, Fore.RESET))


def check_dcm2niix():
    """Check dcm2niix software."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present('dcm2niix'):
        raise ClinicaMissingDependencyError(
            '%s\n[Error] Clinica could not find dcm2niix software in your PATH environment: '
            'this can be downloaded or installed from https://github.com/rordenlab/dcm2niix.%s' %
            (Fore.RED, Fore.RESET))


def check_ants(version_requirements=None):
    """Check ANTs software."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable('ANTSPATH', 'ANTs')

    list_binaries = ['N4BiasFieldCorrection', 'antsRegistrationSyNQuick.sh']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] Clinica could not find ANTs software: the %s command is not present '
                'in your PATH environment.%s' % (Fore.RED, binary, Fore.RESET))


def check_freesurfer(version_requirements=None):
    """Check FreeSurfer software."""
    from colorama import Fore
    import nipype.interfaces.freesurfer as freesurfer
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable('FREESURFER_HOME', 'FreeSurfer')

    list_binaries = ['mri_convert', 'recon-all']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] Clinica could not find FreeSurfer software: the %s command is not present in your PATH '
                'environment: did you have the line source ${FREESURFER_HOME}/'
                'SetUpFreeSurfer.sh in your configuration file?%s' % (Fore.RED, binary, Fore.RESET))


def check_fsl(version_requirements=None):
    """Check FSL software."""
    import nipype.interfaces.fsl as fsl
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError
    from clinica.utils.stream import cprint

    check_environment_variable('FSLDIR', 'FSL')

    try:
        if fsl.Info.version().split(".") < ['5', '0', '5']:
            raise ClinicaMissingDependencyError(
                '%sFSL version must be greater than 5.0.5%s'
                % (Fore.RED, Fore.RESET))
    except Exception as e:
        cprint(str(e))

    list_binaries = ['bet', 'flirt', 'fast', 'first']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] Clinica could not find FSL software: the %s command is not present in your '
                'PATH environment.%s' % (Fore.RED, binary, Fore.RESET))


def check_mrtrix(version_requirements=None):
    """Check MRtrix software."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    check_environment_variable('MRTRIX_HOME', 'MRtrix')

    list_binaries = ['transformconvert', 'mrtransform', 'dwi2response', 'tckgen']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] Clinica could not find MRtrix software: the %s command is not present in your '
                'PATH environment.%s' % (Fore.RED, binary, Fore.RESET))


def check_petpvc(version_requirements=None):
    """Check PETPVC software."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    list_binaries = ['petpvc', 'pvc_diy', 'pvc_gtm', 'pvc_iy', 'pvc_labbe', 'pvc_make4d', 'pvc_mg',
                     'pvc_mtc', 'pvc_rbv', 'pvc_relabel', 'pvc_rl', 'pvc_simulate', 'pvc_stc', 'pvc_vc']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] Clinica could not find PETPVC software: the %s command is not present in your '
                'PATH environment.%s' % (Fore.RED, binary, Fore.RESET))


def check_spm(version_requirements=None):
    """Check SPM software."""
    check_environment_variable('SPM_HOME', 'SPM')

    # list_binaries = ['matlab']
    # for binary in list_binaries:
    #     if not is_binary_present(binary):
    #         raise RuntimeError(
    #             '%s from SPM Software is not present in your '
    #             'PATH environment.' % binary)


def check_matlab():
    """Check Matlab toolbox."""
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaMissingDependencyError

    if not is_binary_present("matlab"):
        raise ClinicaMissingDependencyError(
            '%sMatlab was not found in PATH environment. Did you add it?%s'
            % (Fore.RED, Fore.RESET))
