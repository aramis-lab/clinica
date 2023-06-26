def rename_into_caps(
    in_bids_dwi: str,
    fname_dwi: str,
    fname_bval: str,
    fname_bvec: str,
    fname_brainmask: str,
    fname_magnitude: str,
    fname_fmap: str,
    fname_smoothed_fmap: str,
) -> tuple:
    """Rename the outputs of the pipelines into CAPS.

    Parameters
    ----------
    in_bids_dwi : str
        Path to input BIDS DWI to extract the <source_file>

    fname_dwi : str
        Name of preprocessed DWI file.

    fname_bval : str
        Name of preprocessed bval file.

    fname_bvec : str
        Name of preprocessed bvec file.

    fname_brainmask : str
        Name of B0 mask file.

    fname_smoothed_fmap : str
        Name of smoothed (calibrated) fmap file on b0 space.

    fname_fmap : str
        Name of calibrated fmap file on b0 space.

    fname_magnitude : str
        Name of magnitude image file on b0 space.

    Returns
    -------
    Tuple[str, str, str, str, str, str, str] :
        The different outputs in CAPS format.
    """
    from clinica.utils.dwi import rename_files

    return rename_files(
        in_bids_dwi,
        {
            fname_dwi: "_space-b0_preproc.nii.gz",
            fname_bval: "_space-b0_preproc.bval",
            fname_bvec: "_space-b0_preproc.bval",
            fname_brainmask: "_space-b0_brainmask.nii.gz",
            fname_magnitude: "_space-b0_magnitude1.nii.gz",
            fname_fmap: "_space-b0_fmap.nii.gz",
            fname_smoothed_fmap: "_space-b0_fwhm-4_fmap.nii.gz",
        },
    )


def get_grad_fsl(b_vectors_filename: str, b_values_filename) -> tuple:
    return b_vectors_filename, b_values_filename


def init_input_node(
    dwi, bvec, bval, dwi_json, fmap_magnitude, fmap_phasediff, fmap_phasediff_json
):
    """Initialize pipeline (read JSON, check files and print begin message)."""
    import datetime

    import nibabel as nib

    from clinica.utils.dwi import bids_dir_to_fsl_dir, check_dwi_volume
    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        get_subject_id,
        handle_missing_keys_dwi,
    )
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(dwi)

    # Check that the number of DWI, bvec & bval are the same
    try:
        check_dwi_volume(dwi, bvec, bval)
    except ValueError as e:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_msg = f"[{now}] Error: Number of DWIs, b-vals and b-vecs mismatch for {image_id.replace('_', ' | ')}"
        cprint(error_msg, lvl="error")
        raise ValueError(e)

    # Check that PhaseDiff and magnitude1 have the same header
    # Otherwise, FSL in FugueExtrapolationFromMask will crash
    img_phasediff = nib.load(fmap_phasediff)
    img_magnitude = nib.load(fmap_magnitude)
    if img_phasediff.shape != img_magnitude.shape:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_msg = (
            f"[{now}] Error: Headers of PhaseDiff and Magnitude1 are not the same "
            f"for {image_id.replace('_', ' | ')} ({img_phasediff.shape} vs {img_magnitude.shape})"
        )
        cprint(error_msg, lvl="error")
        raise NotImplementedError(error_msg)

    # Read metadata from DWI JSON file:
    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        dwi_json,
        [
            "TotalReadoutTime",
            "PhaseEncodingDirection",
        ],
        handle_missing_keys=handle_missing_keys_dwi,
    )
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    # Read metadata from PhaseDiff JSON file:
    [echo_time1, echo_time2] = extract_metadata_from_json(
        fmap_phasediff_json, ["EchoTime1", "EchoTime2"]
    )
    delta_echo_time = abs(echo_time2 - echo_time1)

    # Print begin message
    print_begin_image(
        image_id,
        ["TotalReadoutTime", "PhaseEncodingDirection", "DeltaEchoTime"],
        [str(total_readout_time), phase_encoding_direction, str(delta_echo_time)],
    )

    return (
        image_id,
        dwi,
        bvec,
        bval,
        total_readout_time,
        phase_encoding_direction,
        fmap_magnitude,
        fmap_phasediff,
        delta_echo_time,
    )


def print_end_pipeline(image_id, final_file):
    """Display end message for `image_id` when `final_file` is connected."""
    from clinica.utils.ux import print_end_image

    print_end_image(image_id)


def rads2hz(in_file: str, delta_te: float, out_file: str = None) -> str:
    """Convert input phase difference map to Hz.

    Parameters
    ----------
    in_file : str
        Path to the phase difference map image.

    delta_te : float
        The DeltaEchoTime from BIDS specifications.

    out_file : str, optional
        Path to output file. If None, the output image
        will be written in current folder. The file name
        will have the same base name as in_file, but with
        the "_radsec.nii.gz" suffix.

    Returns
    -------
    out_file : str
        Path to output file.
    """
    import math
    import os

    import nibabel as nb
    import numpy as np

    if out_file is None:
        fname, fext = os.path.splitext(os.path.basename(in_file))
        if fext == ".gz":
            fname, _ = os.path.splitext(fname)
        out_file = os.path.abspath(f"./{fname}_radsec.nii.gz")

    im = nb.load(in_file)
    data = im.get_data().astype(np.float32) * (1.0 / (float(delta_te) * 2 * math.pi))
    nb.Nifti1Image(data, im.affine, im.header).to_filename(out_file)

    return out_file
