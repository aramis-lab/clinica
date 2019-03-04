# coding: utf8


def ants_registration_syn_quick(fix_image, moving_image):
    import os
    import os.path as op

    image_warped = op.abspath('SyN_QuickWarped.nii.gz')
    affine_matrix = op.abspath('SyN_Quick0GenericAffine.mat')
    warp = op.abspath('SyN_Quick1Warp.nii.gz')
    inverse_warped = op.abspath('SyN_QuickInverseWarped.nii.gz')
    inverse_warp = op.abspath('SyN_Quick1InverseWarp.nii.gz')

    cmd = 'antsRegistrationSyNQuick.sh -t br -d 3 -f %s  -m %s -o SyN_Quick' \
          % (fix_image, moving_image)
    os.system(cmd)

    return image_warped, affine_matrix, warp, inverse_warped, inverse_warp


def ants_combine_transform(in_file, transforms_list, reference):
    import os
    import os.path as op

    out_warp = op.abspath('out_warp.nii.gz')

    transforms = ""
    for trans in transforms_list:
        transforms += " " + trans
    cmd = 'antsApplyTransforms -o [out_warp.nii.gz,1] -i %s -r %s -t %s' \
          % (in_file, reference, transforms)
    os.system(cmd)

    return out_warp


def dwi_container_from_filename(bids_dwi_filename):
    """ Generate subjects/sub-<participant_id>/ses-<session_id> folder
    from BIDS filename.
    """
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', bids_dwi_filename)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format. '
            + 'It does not contain the subject and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session)


def rename_into_caps(in_bids_dwi,
                     fname_dwi, fname_bval, fname_bvec, fname_brainmask):
    """
    Rename the outputs of the pipelines into CAPS format namely:
    <source_file>_space-T1w_preproc[.nii.gz|bval|bvec]

    Args:
        in_bids_dwi (str): Input BIDS DWI to extract the <source_file>
        fname_dwi (str): Preprocessed DWI file.
        fname_bval (str): Preprocessed bval.
        fname_bvec (str): Preprocessed bvec.
        fname_brainmask (str): B0 mask.

    Returns:
        The different outputs in CAPS format
    """
    from nipype.utils.filemanip import split_filename
    from nipype.interfaces.utility import Rename
    import os

    # Extract <source_file> in format sub-CLNC01_ses-M00[_acq-label]_dwi
    _, source_file_dwi, _ = split_filename(in_bids_dwi)

    # Extract base path from fname:
    base_dir_dwi, _, _ = split_filename(fname_dwi)
    base_dir_bval, _, _ = split_filename(fname_bval)
    base_dir_bvec, _, _ = split_filename(fname_bvec)
    base_dir_brainmask, _, _ = split_filename(fname_brainmask)

    # Rename into CAPS DWI:
    rename_dwi = Rename()
    rename_dwi.inputs.in_file = fname_dwi
    rename_dwi.inputs.format_string = source_file_dwi + "_space-T1w_preproc.nii.gz"
    out_caps_dwi = rename_dwi.run()

    # Rename into CAPS bval:
    rename_bval = Rename()
    rename_bval.inputs.in_file = fname_bval
    rename_bval.inputs.format_string = source_file_dwi + "_space-T1w_preproc.bval"
    out_caps_bval = rename_bval.run()

    # Rename into CAPS bvec:
    rename_bvec = Rename()
    rename_bvec.inputs.in_file = fname_bvec
    rename_bvec.inputs.format_string = source_file_dwi + "_space-T1w_preproc.bvec"
    out_caps_bvec = rename_bvec.run()

    # Rename into CAPS DWI:
    rename_brainmask = Rename()
    rename_brainmask.inputs.in_file = fname_brainmask
    rename_brainmask.inputs.format_string = source_file_dwi + "_space-T1w_brainmask.nii.gz"
    out_caps_brainmask = rename_brainmask.run()

    return out_caps_dwi.outputs.out_file, out_caps_bval.outputs.out_file, \
        out_caps_bvec.outputs.out_file, out_caps_brainmask.outputs.out_file


def print_begin_pipeline(in_bids_or_caps_file):
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s Running pipeline for %s...' % (
        Fore.BLUE, now, Fore.RESET, m.group(0)))


def print_end_pipeline(in_bids_or_caps_file, final_file):
    """
    """
    from clinica.utils.stream import cprint
    import re
    import datetime
    from colorama import Fore

    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)',
                  in_bids_or_caps_file)
    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format.')
    now = datetime.datetime.now().strftime('%H:%M:%S')

    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, m.group(0)))
