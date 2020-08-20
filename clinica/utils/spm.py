# coding: utf8

"""
This module contains SPM utilities.
"""

INDEX_TISSUE_MAP = {
    1: 'graymatter',
    2: 'whitematter',
    3: 'csf',
    4: 'bone',
    5: 'softtissue',
    6: 'background'
}


def check_spm_home():
    """Check and get SPM_HOME environment variable if present."""
    import os
    import platform
    from colorama import Fore
    from .check_dependency import check_environment_variable
    from .exceptions import ClinicaMissingDependencyError

    spm_home = check_environment_variable('SPM_HOME', 'SPM')

    spm_standalone_home = os.environ.get('SPMSTANDALONE_HOME', '')
    if spm_standalone_home:
        if not os.path.isdir(spm_standalone_home):
            raise ClinicaMissingDependencyError(
                '%s\n[Error] The SPMSTANDALONE_HOME environment variable you gave is not a folder (content: %s).%s'
                % (Fore.RED, spm_standalone_home, Fore.RESET))
        if platform.system() == 'Darwin':
            spm_home = os.path.join(spm_standalone_home, 'spm12.app', 'Contents', 'MacOS', 'spm12_mcr')
        else:
            spm_home = os.path.join(spm_standalone_home, 'spm12_mcr')

    return spm_home


def get_tpm():
    """
    Get Tissue Probability Map (TPM) from SPM
    Args
        None
    Returns
        TPM.nii from SPM
    """
    import os
    from glob import glob
    from os.path import join

    spm_home = os.getenv('SPM_HOME')

    if not spm_home:
        # Try MCR to get a hint on SPM location
        spm_home = os.getenv('SPMSTANDALONE_HOME')

    if not spm_home:
        raise RuntimeError('Could not determine location of your SPM installation. '
                           + 'Neither $SPM_HOME or $SPMSTANDALONE_HOME are present '
                             'in your environment')

    tpm_file_glob = glob(join(spm_home, '**/TPM.nii'), recursive=True)
    if len(tpm_file_glob) == 0:
        raise RuntimeError('No file found for TPM.nii in your $SPM_HOME in ' + spm_home)
    if len(tpm_file_glob) > 1:
        error_str = 'Multiple files found for TPM.nii in your SPM_HOME ' + spm_home + ':'
        for file in tpm_file_glob:
            error_str += '\n\t' + file
        raise RuntimeError(error_str)
    return tpm_file_glob[0]


def spm_standalone_is_available():
    """
        Tells if SPM standalone can be used
    Returns:
        True if SPM standalone is detected, False otherwise. Note that it does not guarentee that SPM (classical) is
        up and running in the system.
    """
    import os
    from os.path import isdir, expandvars

    use_spm_stand = False
    if all(elem in os.environ.keys() for elem in ['SPMSTANDALONE_HOME', 'MCR_HOME']):
        if isdir(expandvars('$SPMSTANDALONE_HOME')) and isdir(expandvars('$MCR_HOME')):
            use_spm_stand = True
        else:
            raise FileNotFoundError('[Error] $SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder')
    return use_spm_stand


def use_spm_standalone():
    """
    Use SPM Standalone with MATLAB Common Runtime
    """
    import os
    from colorama import Fore
    import platform
    from clinica.utils.stream import cprint
    from nipype.interfaces import spm
    # This section of code determines whether to use SPM standalone or not
    if all(elem in os.environ.keys() for elem in ['SPMSTANDALONE_HOME', 'MCR_HOME']):
        spm_standalone_home = os.getenv('SPMSTANDALONE_HOME')
        mcr_home = os.getenv('MCR_HOME')
        if os.path.exists(spm_standalone_home) and os.path.exists(mcr_home):
            cprint(Fore.GREEN + 'SPM standalone has been found and will be used in this pipeline' + Fore.RESET)
            if platform.system().lower().startswith('darwin'):
                matlab_cmd = ('cd %s && ./run_spm12.sh %s script' %
                              (spm_standalone_home, mcr_home))
            elif platform.system().lower().startswith('linux'):
                matlab_cmd = ('%s %s script' %
                              (os.path.join(spm_standalone_home, 'run_spm12.sh'), mcr_home))
            else:
                raise SystemError('Clinica only support macOS and Linux')
            spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
            cprint("Using SPM standalone version %s" % spm.SPMCommand().version)
        else:
            raise FileNotFoundError('$SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder ')
