def get_tpm():
    """Get Tissue Probability Map (TPM) from SPM.

    Returns:
        str: TPM.nii from SPM
    """
    import os
    from glob import glob
    from os.path import join

    spm_home = os.getenv("SPM_HOME")

    if not spm_home:
        # Try MCR to get a hint on SPM location
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
