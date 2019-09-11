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


def rename_into_caps(in_bids_dwi,
                     fname_dwi, fname_bval, fname_bvec, fname_brainmask,
                     fname_magnitude, fname_fmap, fname_smoothed_fmap):
    """
    Rename the outputs of the pipelines into CAPS format namely:
    <source_file>_space-b0_preproc{.nii.gz|bval|bvec}
    <source_file>_space-b0_brainmask.nii.gz
    <source_file>_space-b0_magnitude1.nii.gz
    <source_file>_space-b0[_fwhm-<fwhmm>]_fmap.nii.gz

    Args:
        in_bids_dwi (str): Input BIDS DWI to extract the <source_file>
        fname_dwi (str): Preprocessed DWI.
        fname_bval (str): Preprocessed DWI.
        fname_bvec (str): Preprocessed DWI.
        fname_brainmask (str): B0 mask.
        fname_smoothed_fmap (str): Smoothed (calibrated) fmap on b0 space.
        fname_fmap (str): Calibrated fmap on b0 space.
        fname_magnitude (str): Magnitude image on b0 space.

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
    base_dir_smoothed_fmap, _, _ = split_filename(fname_smoothed_fmap)
    base_dir_calibrated_fmap, _, _ = split_filename(fname_fmap)
    base_dir_magnitude, _, _ = split_filename(fname_magnitude)

    # Rename into CAPS DWI:
    rename_dwi = Rename()
    rename_dwi.inputs.in_file = fname_dwi
    rename_dwi.inputs.format_string = os.path.join(
        base_dir_dwi, source_file_dwi + "_space-b0_preproc.nii.gz")
    out_caps_dwi = rename_dwi.run()

    # Rename into CAPS bval:
    rename_bval = Rename()
    rename_bval.inputs.in_file = fname_bval
    rename_bval.inputs.format_string = os.path.join(
        base_dir_bval, source_file_dwi + "_space-b0_preproc.bval")
    out_caps_bval = rename_bval.run()

    # Rename into CAPS bvec:
    rename_bvec = Rename()
    rename_bvec.inputs.in_file = fname_bvec
    rename_bvec.inputs.format_string = os.path.join(
        base_dir_bvec, source_file_dwi + "_space-b0_preproc.bvec")
    out_caps_bvec = rename_bvec.run()

    # Rename into CAPS brainmask:
    rename_brainmask = Rename()
    rename_brainmask.inputs.in_file = fname_brainmask
    rename_brainmask.inputs.format_string = os.path.join(
        base_dir_brainmask, source_file_dwi + "_space-b0_brainmask.nii.gz")
    out_caps_brainmask = rename_brainmask.run()

    # Rename into CAPS magnitude:
    rename_magnitude = Rename()
    rename_magnitude.inputs.in_file = fname_magnitude
    rename_magnitude.inputs.format_string = os.path.join(
        base_dir_magnitude, source_file_dwi + "_space-b0_magnitude1.nii.gz")
    out_caps_magnitude = rename_magnitude.run()

    # Rename into CAPS fmap:
    rename_calibrated_fmap = Rename()
    rename_calibrated_fmap.inputs.in_file = fname_fmap
    rename_calibrated_fmap.inputs.format_string = os.path.join(
        base_dir_calibrated_fmap, source_file_dwi + "_space-b0_fmap.nii.gz")
    out_caps_fmap = rename_calibrated_fmap.run()

    # Rename into CAPS smoothed fmap:
    rename_smoothed_fmap = Rename()
    rename_smoothed_fmap.inputs.in_file = fname_smoothed_fmap
    rename_smoothed_fmap.inputs.format_string = os.path.join(
        base_dir_smoothed_fmap, source_file_dwi + "_space-b0_fwhm-4_fmap.nii.gz")
    out_caps_smoothed_fmap = rename_smoothed_fmap.run()

    return (out_caps_dwi.outputs.out_file,
            out_caps_bval.outputs.out_file,
            out_caps_bvec.outputs.out_file,
            out_caps_brainmask.outputs.out_file,
            out_caps_magnitude.outputs.out_file,
            out_caps_fmap.outputs.out_file,
            out_caps_smoothed_fmap.outputs.out_file,)


def get_grad_fsl(bvec, bval):
    grad_fsl = (bvec, bval)
    return grad_fsl


def init_input_node(dwi, bvec, bval, dwi_json,
                    fmap_magnitude, fmap_phasediff, fmap_phasediff_json):
    """Extract "sub-<participant_id>_ses-<session_label>" from input node and print begin message."""
    import datetime
    from colorama import Fore
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.stream import cprint
    from clinica.utils.io import get_subject_id, extract_metadata_from_json
    from clinica.utils.dwi import check_dwi_volume
    from clinica.utils.epi import bids_dir_to_fsl_dir

    # Check the image IDs for each file are the same
    image_ids = [
        get_subject_id(dwi),
        get_subject_id(bvec),
        get_subject_id(bval),
        get_subject_id(dwi_json),
        get_subject_id(fmap_magnitude),
        get_subject_id(fmap_phasediff),
        get_subject_id(fmap_phasediff_json)
    ]
    if not len(set(image_ids)) == 1:
        raise ClinicaException('<image_id> from input files mismatch (found: %s)' % image_ids)

    # Check that the number of DWI, bvec & bval are the same
    check_dwi_volume(dwi, bvec, bval)

    # Read metadata from DWI JSON file:
    [total_readout_time, phase_encoding_direction] = \
        extract_metadata_from_json(dwi_json, ['TotalReadoutTime', 'PhaseEncodingDirection'])
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    # Read metadata from PhaseDiff JSON file:
    [echo_time_1, echo_time_2] = \
        extract_metadata_from_json(fmap_phasediff_json, ['EchoTime1', 'EchoTime2'])
    delta_echo_time = abs(echo_time_2 - echo_time_1)

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s Running pipeline for %s '
           '(TotalReadoutTime = %s, PhaseEncodingDirection = %s, DeltaEchoTime = %s)' %
           (Fore.BLUE, now, Fore.RESET, image_ids[0].replace('_', '|'),
            total_readout_time, phase_encoding_direction, delta_echo_time))
    return (image_ids[0], dwi, bvec, bval, total_readout_time, phase_encoding_direction,
            fmap_magnitude, fmap_phasediff, delta_echo_time)


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.stream import cprint
    import datetime
    from colorama import Fore

    now = datetime.datetime.now().strftime('%H:%M:%S')
    cprint('%s[%s]%s ...%s has completed.' % (
        Fore.GREEN, now, Fore.RESET, image_id.replace('_', '|')))
