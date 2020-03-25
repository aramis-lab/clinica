# coding: utf8


"""This module contains utilities to check dependencies of the different
neuroimaging tools."""


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

    try:
        devnull = open(os.devnull)
        subprocess.Popen([binary], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
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

    check_environment_variable('MRTRIX_HOME', 'MRtrix')

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


def check_cat12():
    """
    Check installation of CAT12 (mostly used to provide atlases)
    """
    import os
    from colorama import Fore
    from .exceptions import ClinicaMissingDependencyError
    from .spm import check_spm_home

    spm_home = check_spm_home()

    cat12 = os.path.join(spm_home, 'toolbox', 'cat12')
    if not os.path.exists(cat12):
        raise ClinicaMissingDependencyError(
            '%s\n[Error] CAT12 is not installed in your SPM folder. It should be in %s.%s'
            % (Fore.RED, cat12, Fore.RESET))
    return cat12


def verify_cat12_atlases(atlas_list):
    """ If the user wants to use any of the atlases of CAT12 and has not installed it, we just remove it from the list
        of the computed atlases

    :param atlas_list: list of atlases given by the user
    :return: atlas_list updated with elements removed if CAT12 is not found, or return atlas_list if
             CAT12 has been found. If none of the CAT12 atlases are in atlas_list, function returns atlas_list
    """
    import sys
    from colorama import Fore
    from time import sleep
    from .exceptions import ClinicaMissingDependencyError

    cat12_atlases = ['Hammers', 'Neuromorphometrics', 'LPBA40']
    if any(atlas in cat12_atlases for atlas in atlas_list):
        try:
            check_cat12()
            atlas_list_updated = atlas_list
        except ClinicaMissingDependencyError:
            print(Fore.YELLOW + '[Warning] CAT12 is not installed in your system. Atlas statistics computed at the '
                  + 'end of this pipeline will not include any of Hammers, Neuromorphometrics and LPBA40'
                  + ' atlases.' + Fore.RESET)
            print('Pipeline will start in 15 sec. You can cancel it right now using Ctrl + C')
            print('You can download CAT12 at this address: ' + Fore.BLUE + 'http://www.neuro.uni-jena.de/cat/'
                  + Fore.RESET)
            atlas_list_updated = [atl for atl in atlas_list if atl not in cat12_atlases]
            try:
                sleep(15)
            except KeyboardInterrupt:
                print('\nClinica is now exiting.')
                sys.exit()
    else:
        atlas_list_updated = atlas_list
    return atlas_list_updated
