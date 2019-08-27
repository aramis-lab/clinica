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
#    cprint("Binary '%s' has been detected." % binary)
    return True


def check_ants():
    """
    Check ANTs software.

    This function checks if ANTs is present (ANTSPATH, binaries, & scripts).
    """
    import os
    from clinica.utils.stream import cprint

    try:
        antspath = os.environ.get('ANTSPATH', '')
        if not antspath:
            raise RuntimeError('ANTSPATH variable is not set')
    except Exception as e:
        cprint(str(e))

    list_binaries = ['N4BiasFieldCorrection', 'antsRegistrationSyNQuick.sh']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise RuntimeError(
                '%s from ANTs Software is not present '
                'in your PATH environment.' % binary)

    cprint('ANTs has been detected')


def check_freesurfer():
    """
    Check FreeSurfer software.

    This function checks if FreeSurfer is present (FREESURFER_HOME & binaries).
    """
    import os
    from clinica.utils.stream import cprint

    try:
        freesurfer_home = os.environ.get('FREESURFER_HOME', '')
        if not freesurfer_home:
            raise RuntimeError('FREESURFER_HOME variable is not set')
    except Exception as e:
        cprint(str(e))

    list_binaries = ['mri_convert', 'recon-all']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise RuntimeError(
                '%s from FreeSurfer Software is not present in your PATH '
                'environment: did you source ${FREESURFER_HOME}/'
                'SetUpFreeSurfer.sh ?' % binary)

    cprint('FreeSurfer has been detected')


def check_fsl():
    """
    Check FSL software.

    This function checks if FSL is present (FSLDIR & binaries) and if the
    version of FSL is recent.
    """
    import os
    import nipype.interfaces.fsl as fsl

    from clinica.utils.stream import cprint

    try:
        fsl_dir = os.environ.get('FSLDIR', '')
        if not fsl_dir:
            raise RuntimeError('FSLDIR variable is not set')
    except Exception as e:
        cprint(str(e))

    try:
        if fsl.Info.version().split(".") < ['5', '0', '5']:
            raise RuntimeError('FSL version must be greater than 5.0.5')
    except Exception as e:
        cprint(str(e))

    list_binaries = ['bet', 'flirt', 'fast', 'first']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise RuntimeError(
                '%s from FSL Software is not present in your '
                'PATH environment.' % binary)

    cprint('FSL has been detected')


def check_mrtrix():
    """
    Check MRtrix software.

    This function checks if MRtrix is present (MRTRIX_HOME & binaries).
    """
    import os
    from clinica.utils.stream import cprint

    try:
        mrtrix_home = os.environ.get('MRTRIX_HOME', '')
        if not mrtrix_home:
            raise RuntimeError('MRTRIX_HOME variable is not set')
    except Exception as e:
        print(str(e))

    list_binaries = ['transformconvert', 'mrtransform', 'dwi2response', 'tckgen']
    for binary in list_binaries:
        if not is_binary_present(binary):
            raise RuntimeError(
                '%s from MRtrix Software is not present in your '
                'PATH environment.' % binary)

    cprint('MRtrix has been detected')


def check_spm():
    """
    Check SPM software.

    This function checks if SPM is present (SPM_HOME).
    """
    import os
    from clinica.utils.stream import cprint

    try:
        spm_home = os.environ.get('SPM_HOME', '')
        if not spm_home:
            raise RuntimeError('SPM_HOME variable is not set')
    except Exception as e:
        print(str(e))

    # list_binaries = ['matlab']
    # for binary in list_binaries:
    #     if not is_binary_present(binary):
    #         raise RuntimeError(
    #             '%s from SPM Software is not present in your '
    #             'PATH environment.' % binary)

    cprint('SPM has been detected')


def check_matlab():
    """
    Check matlab toolbox.

    This function checks if matlab is present (matlab for linux and MATLABCMD for mac).
    """
    import os
    import sys
    from clinica.utils.stream import cprint

    if not is_binary_present("matlab"):
        raise RuntimeError('Matlab was not found in PATH environment. Did you add it?')
    cprint('Matlab has been detected')


def check_cat12():
    """
    Check installation of cat12 (mostly used to provide atlases)

    """
    import os
    import platform
    from clinica.utils.stream import cprint

    if "SPMSTANDALONE_HOME" in os.environ:
        if platform.system() == 'Darwin':
            SPM_HOME = os.environ['SPMSTANDALONE_HOME'] + "/spm12.app/Contents/MacOS/spm12_mcr"
        else:
            SPM_HOME = os.environ['SPMSTANDALONE_HOME'] + "/spm12_mcr"
    elif "SPM_HOME" in os.environ:
        SPM_HOME = os.environ['SPM_HOME']
    else:
        raise RuntimeError('SPM not installed. $SPM_HOME variable not found in environnement')

    if not os.path.isdir(SPM_HOME):
        raise RuntimeError('SPM and cat12 are not installed.' + SPM_HOME + ' does not exist')

    if not os.path.exists(os.path.join(SPM_HOME, 'toolbox', 'cat12')):
        raise RuntimeError('cat12 is not installed in your SPM folder :'
                           + str(os.path.join(SPM_HOME, 'toolbox', 'cat12')))
    cprint('cat12 has been detected')


def verify_cat12_atlases(atlas_list):
    """ If the user wants to use any of the atlases of cat12 and has not installed it, we just remove it from the list
        of the computed atlases

    :param atlas_list: list of atlases given by the user
    :return: atlas_list updated with elements removed if cat12 is not found, or return atlas_list if
             cat12 has been found. If none of the cat12 atlases are in atlas_list, function returns atlas_list
    """
    import sys
    from colorama import Fore
    from time import sleep

    cat12_atlases = ['Hammers', 'Neuromorphometrics', 'LPBA40']
    if any(atlas in cat12_atlases for atlas in atlas_list):
        try:
            check_cat12()
            atlas_list_updated = atlas_list
        except RuntimeError:
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
