# coding: utf8


def dwi_container_from_filename(dwi_filename):
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', dwi_filename)

    if m is None:
        raise ValueError(
            'Input filename is not in a BIDS or CAPS compliant format. '
            + 'It does not contain the subject and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session)


def parameters_from_dwi_metadata(dwi_json):
    """
    Extract BIDS metadata (EffectiveEchoSpacing, etc.) from a DWI .json file.

    Args:
        dwi_json: Path to a BIDS compliant .json file

    Returns:

    """
    import json
    import os

    from clinica.utils.stream import cprint

    if not os.path.exists(dwi_json):
        raise IOError('DWI .json file (%s) does not exist.' % dwi_json)

    with open(dwi_json) as data_frame:
        data = json.load(data_frame)

    effective_echo_spacing = data['EffectiveEchoSpacing']
    phase_encoding_direction = data['PhaseEncodingDirection']

    cprint(dwi_container_from_filename(dwi_json) + " (DWI JSON): EffectiveEchoSpacing=" + str(effective_echo_spacing) +
           ", PhaseEncodingDirection=" + str(phase_encoding_direction) + " (fname = " + dwi_json + ")")

    return [effective_echo_spacing, phase_encoding_direction]


def delta_echo_time_from_bids_fmap(fmap_json):
    """
    Extract BIDS metadata (EchoTime1, etc.) from a fmap .json file.

    Args:
        fmap_json: Path to a BIDS compliant .json file

    Returns:

    """
    import json
    import os

    from clinica.utils.stream import cprint

    if not os.path.exists(fmap_json):
        raise IOError('fmap .json file (%s) does not exist.' % fmap_json)

    with open(fmap_json) as data_frame:
        data = json.load(data_frame)

    echo_time_1 = data['EchoTime1']
    echo_time_2 = data['EchoTime2']
    delta_echo_time = echo_time_2 - echo_time_1

    cprint(dwi_container_from_filename(fmap_json) + " (FMap JSON): EchoTime1=" + str(echo_time_1) +
           ", EchoTime2=" + str(echo_time_2) + " (Delta = " + str(delta_echo_time) + ", fname=" + fmap_json + ")")

    return delta_echo_time


def rename_into_caps(in_bids_dwi,
                     fname_dwi, fname_bval, fname_bvec, fname_brainmask):
    """
    Rename the outputs of the pipelines into CAPS format namely:
    <source_file>_space-T1w_preproc[.nii.gz|bval|bvec]

    Args:
        in_bids_dwi (str): Input BIDS DWI to extract the <source_file>
        fname_dwi (str): Preprocessed DWI.
        fname_bval (str): Preprocessed DWI.
        fname_bvec (str): Preprocessed DWI.
        fname_brainmask (str): B0 mask.

    Returns:
        The different outputs in CAPS format
    """
    from nipype.utils.filemanip import split_filename
    from nipype.interfaces.utility import Rename
    import os

    # Extract <source_file> in format sub-CLNC01_ses-M00_[acq-label]_dwi
    _, source_file_dwi, _ = split_filename(in_bids_dwi)

    # Extract base path from fname:
    base_dir_dwi, _, _ = split_filename(fname_dwi)
    base_dir_bval, _, _ = split_filename(fname_bval)
    base_dir_bvec, _, _ = split_filename(fname_bvec)
    base_dir_brainmask, _, _ = split_filename(fname_brainmask)

    # Rename into CAPS DWI :
    rename_dwi = Rename()
    rename_dwi.inputs.in_file = fname_dwi
    rename_dwi.inputs.format_string = os.path.join(
        base_dir_dwi, source_file_dwi + "_space-b0_preproc.nii.gz")
    out_caps_dwi = rename_dwi.run()

    # Rename into CAPS bval :
    rename_bval = Rename()
    rename_bval.inputs.in_file = fname_bval
    rename_bval.inputs.format_string = os.path.join(
        base_dir_bval, source_file_dwi + "_space-b0_preproc.bval")
    out_caps_bval = rename_bval.run()

    # Rename into CAPS DWI :
    rename_bvec = Rename()
    rename_bvec.inputs.in_file = fname_bvec
    rename_bvec.inputs.format_string = os.path.join(
        base_dir_bvec, source_file_dwi + "_space-b0_preproc.bvec")
    out_caps_bvec = rename_bvec.run()

    # Rename into CAPS DWI :
    rename_brainmask = Rename()
    rename_brainmask.inputs.in_file = fname_brainmask
    rename_brainmask.inputs.format_string = os.path.join(
        base_dir_brainmask, source_file_dwi + "_space-b0_brainmask.nii.gz")
    out_caps_brainmask = rename_brainmask.run()

    return out_caps_dwi.outputs.out_file, out_caps_bval.outputs.out_file, \
        out_caps_bvec.outputs.out_file, out_caps_brainmask.outputs.out_file
