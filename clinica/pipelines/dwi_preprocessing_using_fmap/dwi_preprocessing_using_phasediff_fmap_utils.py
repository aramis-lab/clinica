def rename_into_caps(
    dwi_filename: str,
    dwi_preproc_filename: str,
    b_values_preproc_filename: str,
    b_vectors_preproc_filename: str,
    b0_brain_mask_filename: str,
    calibrated_magnitude_image_filename: str,
    calibrated_field_map_image_filename: str,
    calibrated_smoothed_field_map_image_filename: str,
) -> tuple:
    """Rename the outputs of the pipelines into CAPS.

    Parameters
    ----------
    dwi_filename : str
        The path to input BIDS DWI to extract the <source_file>.

    dwi_preproc_filename : str
        The path to the preprocessed DWI file.

    b_values_preproc_filename : str
        The path to the preprocessed b-values file.

    b_vectors_preproc_filename : str
        The path to the preprocessed b-vectors file.

    b0_brain_mask_filename : str
        The path to the B0 mask file.

    calibrated_magnitude_image_filename : str
        The path to the magnitude image file on b0 space.

    calibrated_field_map_image_filename : str
        The path to the calibrated fmap file on b0 space.

    calibrated_smoothed_field_map_image_filename : str
        The path to the smoothed (calibrated) fmap file on b0 space.

    Returns
    -------
    Tuple[str, str, str, str, str, str, str] :
        The different outputs in CAPS format.
    """
    from clinica.utils.dwi import rename_files

    return rename_files(
        dwi_filename,
        {
            dwi_preproc_filename: "_space-b0_preproc.nii.gz",
            b_values_preproc_filename: "_space-b0_preproc.bval",
            b_vectors_preproc_filename: "_space-b0_preproc.bval",
            b0_brain_mask_filename: "_space-b0_brainmask.nii.gz",
            calibrated_magnitude_image_filename: "_space-b0_magnitude1.nii.gz",
            calibrated_field_map_image_filename: "_space-b0_fmap.nii.gz",
            calibrated_smoothed_field_map_image_filename: "_space-b0_fwhm-4_fmap.nii.gz",
        },
    )


def get_grad_fsl(b_vectors_filename: str, b_values_filename) -> tuple:
    return b_vectors_filename, b_values_filename


def init_input_node(
    dwi_filename: str,
    b_vectors_filename: str,
    b_values_filename: str,
    dwi_json_filename: str,
    fmap_magnitude_filename: str,
    fmap_phasediff_filename: str,
    fmap_phasediff_json_filename: str,
) -> tuple:
    """Initialize pipeline (read JSON, check files and print begin message).

    Parameters
    ----------
    dwi_filename : str
        The path to the DWI image.

    b_vectors_filename : str
        The path to the b-vectors file.

    b_values_filename : str
        The path to the b-values file.

    dwi_json_filename : str
        Path to the JSON metadata file for the DWI image.

    fmap_magnitude_filename : str
        The path to the (1st) magnitude image in BIDS format.

    fmap_phasediff_filename : str
        The path of the phase difference image in BIDS format.

    fmap_phasediff_json_filename : str
        The path of the phase difference JSON file in BIDS format
        and containing EchoTime1 & EchoTime2 metadata (see BIDS specifications).

    Returns
    -------
    image_id : str
        The subject ID extracted from the dwi image path.

    dwi_filename : str
        The path to the DWI image.

    b_vectors_filename : str
        The path to the b-vectors file.

    b_values_filename : str
        The path to the b-values file.

    total_readout_time : str
        The total readout time extracted from the dwi JSON file.

    phase_encoding_direction : str
        The phase encoding direction for the dwi image, extracted
        from the dwi JSON file.

    fmap_magnitude_filename : str
        The path to the (1st) magnitude image in BIDS format.

    fmap_phasediff_filename : str
        The path of the phase difference image in BIDS format.

    delta_echo_time : str
        The path to the JSON metadata file for the field map image.
    """
    import datetime

    import nibabel as nib

    from clinica.utils.dwi import DWIDataset, bids_dir_to_fsl_dir, check_dwi_volume
    from clinica.utils.filemanip import (
        extract_metadata_from_json,
        get_subject_id,
        handle_missing_keys_dwi,
    )
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(dwi_filename)

    try:
        check_dwi_volume(
            DWIDataset(dwi_filename, b_values_filename, b_vectors_filename)
        )
    except ValueError as e:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_msg = f"[{now}] Error: Number of DWIs, b-vals and b-vecs mismatch for {image_id.replace('_', ' | ')}"
        cprint(error_msg, lvl="error")
        raise ValueError(e)

    # Check that PhaseDiff and magnitude1 have the same header
    # Otherwise, FSL in FugueExtrapolationFromMask will crash
    img_phasediff = nib.load(fmap_phasediff_filename)
    img_magnitude = nib.load(fmap_magnitude_filename)
    if img_phasediff.shape != img_magnitude.shape:
        now = datetime.datetime.now().strftime("%H:%M:%S")
        error_msg = (
            f"[{now}] Error: Headers of PhaseDiff and Magnitude1 are not the same "
            f"for {image_id.replace('_', ' | ')} ({img_phasediff.shape} vs {img_magnitude.shape})"
        )
        cprint(error_msg, lvl="error")
        raise NotImplementedError(error_msg)

    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        dwi_json_filename,
        [
            "TotalReadoutTime",
            "PhaseEncodingDirection",
        ],
        handle_missing_keys=handle_missing_keys_dwi,
    )
    phase_encoding_direction = bids_dir_to_fsl_dir(phase_encoding_direction)

    [echo_time1, echo_time2] = extract_metadata_from_json(
        fmap_phasediff_json_filename, ["EchoTime1", "EchoTime2"]
    )
    delta_echo_time = abs(echo_time2 - echo_time1)

    print_begin_image(
        image_id,
        ["TotalReadoutTime", "PhaseEncodingDirection", "DeltaEchoTime"],
        [str(total_readout_time), phase_encoding_direction, str(delta_echo_time)],
    )

    return (
        image_id,
        dwi_filename,
        b_vectors_filename,
        b_values_filename,
        total_readout_time,
        phase_encoding_direction,
        fmap_magnitude_filename,
        fmap_phasediff_filename,
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
