import math
from pathlib import Path
from typing import Optional

import nibabel as nib
import numpy as np


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
            dwi_preproc_filename: "_space-b0_desc-preproc_dwi.nii.gz",
            b_values_preproc_filename: "_space-b0_desc-preproc_dwi.bval",
            b_vectors_preproc_filename: "_space-b0_desc-preproc_dwi.bvec",
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

    from clinica.utils.dwi import (
        DWIDataset,
        check_dwi_volume,
        get_readout_time_and_phase_encoding_direction,
    )
    from clinica.utils.filemanip import extract_metadata_from_json, get_subject_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(dwi_filename)
    try:
        check_dwi_volume(
            DWIDataset(
                dwi=dwi_filename,
                b_values=b_values_filename,
                b_vectors=b_vectors_filename,
            )
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
    (
        total_readout_time,
        phase_encoding_direction,
    ) = get_readout_time_and_phase_encoding_direction(dwi_json_filename)

    # Read metadata from PhaseDiff JSON file:
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


def convert_phase_difference_to_hertz(
    phase_diff_filename: Path,
    delta_echo_time: float,
    working_dir: Optional[Path] = None,
) -> Path:
    """Convert input phase difference map to Hz.

    Parameters
    ----------
    phase_diff_filename : Path
        The path to the phase difference map image.

    delta_echo_time : float
        The DeltaEchoTime from BIDS specifications.

    working_dir : Path, optional
        The path to the working directory. If None, the output image
        will be written in current folder. The file name
        will have the same base name as in_file, but with
        the "_radsec.nii.gz" suffix.

    Returns
    -------
    out_file : Path
        The path to output file.
    """
    im = nib.load(phase_diff_filename)
    working_dir = working_dir or Path.cwd()
    out_file = working_dir / _get_output_file(phase_diff_filename, "radsec")
    data = im.get_fdata().astype(np.float32) * (
        1.0 / (float(delta_echo_time) * 2 * math.pi)
    )
    nib.Nifti1Image(data, im.affine, im.header).to_filename(out_file)

    return out_file


def _get_output_file(input_file: Path, suffix: str) -> str:
    filename = input_file.stem
    if input_file.suffix == ".gz":
        filename = Path(filename).stem

    return f"{filename}_{suffix}.nii.gz"


def convert_phase_difference_to_hertz_task(
    phase_diff_filename: str,
    delta_echo_time: float,
    working_dir: str = None,
) -> str:
    """Wrapper for Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi_preprocessing_using_fmap.utils import (
        convert_phase_difference_to_hertz,  # noqa
    )

    if working_dir:
        working_dir = Path(working_dir)
    return str(
        convert_phase_difference_to_hertz(
            Path(phase_diff_filename), delta_echo_time, working_dir
        )
    )


def demean_image(
    input_image: Path, mask: Optional[Path] = None, working_dir: Optional[Path] = None
) -> Path:
    """Demean image data inside mask.

    This function was taken from: https://github.com/niflows/nipype1-workflows/

    Parameters
    ----------
    input_image : Path
        The image to demean.

    mask : Path, optional
        If provided, will be used to mask the image.

    working_dir : Path, optional
        The path to the working directory. If None, the output image
        will be written in current folder. The file name
        will have the same base name as in_file, but with
        the "_demean.nii.gz" suffix.

    Returns
    -------
    out_file : Path
        The path to output file.
    """
    image = nib.load(input_image)
    working_dir = working_dir or Path.cwd()
    out_file = working_dir / _get_output_file(input_image, "demean")
    data = image.get_fdata().astype(np.float32)
    if mask:
        mask = nib.load(mask).get_fdata().astype(np.float32)
        mask[mask > 0] = 1.0
        mask[mask < 1] = 0.0
    else:
        mask = np.ones_like(data)
    mean_image = np.median(data[mask == 1].reshape(-1))
    data[mask == 1] = data[mask == 1] - mean_image
    nib.Nifti1Image(data, image.affine, image.header).to_filename(out_file)

    return out_file


def demean_image_task(
    input_image: str, mask: str = None, working_dir: str = None
) -> str:
    """Wrapper for Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi_preprocessing_using_fmap.utils import (
        demean_image,  # noqa
    )

    if working_dir:
        working_dir = Path(working_dir)
    if mask:
        mask = Path(mask)
    return str(demean_image(Path(input_image), mask, working_dir))


def convert_phase_difference_to_rads(
    phase_diff_filename: Path, working_dir: Optional[Path] = None
) -> Path:
    """Converts input phase difference map to rads.

    This function was taken from: https://github.com/niflows/nipype1-workflows/

    Parameters
    ----------
    phase_diff_filename : Path
        The path to the phase difference map image.

    working_dir : Path, optional
        The path to the working directory. If None, the output image
        will be written in current folder. The file name
        will have the same base name as in_file, but with
        the "_rads.nii.gz" suffix.

    Returns
    -------
    out_file : Path
        The path to output file.
    """
    working_dir = working_dir or Path.cwd()
    out_file = working_dir / _get_output_file(phase_diff_filename, "rads")
    in_file = np.atleast_1d(phase_diff_filename).tolist()
    image = nib.load(in_file[0])
    data = image.get_fdata().astype(np.float32)
    header = image.header.copy()

    if len(in_file) == 2:
        data = nib.load(in_file[1]).get_fdata().astype(np.float32) - data
    elif (data.ndim == 4) and (data.shape[-1] == 2):
        data = np.squeeze(data[..., 1] - data[..., 0])
        header.set_data_shape(data.shape[:3])

    data_min = data.min()
    data_max = data.max()
    data = (2.0 * math.pi * (data - data_min) / (data_max - data_min)) - math.pi
    header.set_data_dtype(np.float32)
    header.set_xyzt_units("mm")
    header["datatype"] = 16
    nib.Nifti1Image(data, image.affine, header).to_filename(out_file)

    return out_file


def convert_phase_difference_to_rads_task(
    phase_diff_filename: str, working_dir: str = None
) -> str:
    """Wrapper for Nipype."""
    from pathlib import Path

    from clinica.pipelines.dwi_preprocessing_using_fmap.utils import (
        convert_phase_difference_to_rads,  # noqa
    )

    if working_dir:
        working_dir = Path(working_dir)
    return str(convert_phase_difference_to_rads(Path(phase_diff_filename), working_dir))
