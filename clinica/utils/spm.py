# coding: utf8


def check_spm_home():
    """Check and get SPM_HOME environment variable if present."""
    import os
    import platform
    from .check_dependency import check_environment_variable

    spm_home = check_environment_variable('SPM_HOME', 'SPM')

    spm_standalone_home = os.environ.get('SPMSTANDALONE_HOME', '')
    if spm_standalone_home:
        assert(os.path.isdir(spm_standalone_home))
        if platform.system() == 'Darwin':
            spm_home = os.path.join(spm_standalone_home, 'spm12.app', 'Contents', 'MacOS', 'spm12_mcr')
        else:
            spm_home = os.path.join(spm_standalone_home, 'spm12_mcr')

    return spm_home


def get_tpm():
    """
    Get Tissue Probability Map (TPM) from SPM
    """
    import os
    import platform
    import nipype.interfaces.matlab as matlab
    import nipype.interfaces.spm as spm

    spm_home = os.getenv("SPM_HOME")
    matlab_home = os.getenv("MATLABCMD")
    matlab.MatlabCommand.set_default_matlab_cmd(matlab_home)
    matlab.MatlabCommand.set_default_paths(spm_home)

    if 'SPMSTANDALONE_HOME' in os.environ:
        if 'MCR_HOME' in os.environ:
            matlab_cmd = (os.path.join(os.environ['SPMSTANDALONE_HOME'], 'run_spm12.sh')
                          + ' '
                          + os.environ['MCR_HOME']
                          + ' script')
            spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
            version = spm.SPMCommand().version
        else:
            raise EnvironmentError(
                'MCR_HOME variable not in environment. Although, '
                + 'SPMSTANDALONE_HOME has been found')
    else:
        version = spm.Info.getinfo()

    if version:
        if isinstance(version, dict):
            spm_path = version['path']
            if version['name'] == 'SPM8':
                print(
                    'You are using SPM version 8. The recommended version to use with Clinica is SPM 12. '
                    + 'Please upgrade your SPM toolbox.')
                tissue_map = os.path.join(spm_path, 'toolbox', 'Seg', 'TPM.nii')
            elif version['name'] == 'SPM12':
                tissue_map = os.path.join(spm_path, 'tpm', 'TPM.nii')
            else:
                raise RuntimeError(
                    'SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
        if isinstance(version, str):
            if float(version) >= 12.7169:
                if platform.system() == 'Darwin':
                    tissue_map = os.path.join(str(spm_home), 'spm12.app', 'Contents', 'MacOS', 'spm12_mcr', 'spm12', 'spm12', 'tpm', 'TPM.nii')
                else:
                    # Path depends on version of SPM Standalone
                    if os.path.exists(os.path.join(str(spm_home), 'spm12_mcr', 'spm', 'spm12', 'tpm')):
                        tissue_map = os.path.join(str(spm_home), 'spm12_mcr', 'spm', 'spm12', 'tpm', 'TPM.nii')
                    else:
                        tissue_map = os.path.join(str(spm_home), 'spm12_mcr', 'spm12', 'spm12', 'tpm', 'TPM.nii')
            else:
                raise RuntimeError(
                    'SPM standalone version not supported. Please upgrade SPM standalone.')
    else:
        raise RuntimeError(
            'SPM could not be found. Please verify your SPM_HOME environment variable.')
    return tissue_map
