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

    list_binaries = ['transformconvert', 'mrtransform']
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
    import os, sys
    from clinica.utils.stream import cprint
    # here, we check out the os, basically, clinica works for linux and MAC OS X.
    if sys.platform.startswith('linux'):
        list_binaries = 'matlab'
        if not is_binary_present(list_binaries):
            raise RuntimeError('Your platform is linux, the default command line for Matlab(matlab_cmd) is matlab. You can also export a variable MATLABCMD in your .zshrc or .bashrc file to use a specific version of Matlab.')
        else:
            cprint('matlab has been detected')
    elif sys.platform.startswith('darwin'):
        try:
            if not 'MATLABCMD' in os.environ:
                raise RuntimeError(
                    "Your platform is MAC OS X, the default command line for Matlab(matlab_cmd) is matlab, but it does not work on OS X, you mush export a variable MATLABCMD in your .zshrc or .bashrc, which points to your matlal. Besides, Mac os x will always choose to use OpengGl hardware mode.")
            else:
                cprint('matlab has been detected')
        except Exception as e:
            cprint(str(e))
            exit(1)
    else:
        cprint("Clinica will not work on windows platform ")

