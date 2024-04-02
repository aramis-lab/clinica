"""This module contains utilities common to DWI preprocessing pipelines."""

from os import PathLike
from pathlib import Path
from typing import Optional, Tuple, Union

import nibabel as nib
import numpy as np

from ..utils import DWIDataset

__all__ = [
    "generate_index_file",
    "generate_acq_file",
    "compute_average_b0",
    "check_file",
    "add_suffix_to_filename",
    "get_b0_filter",
    "check_b_value_threshold",
    "get_readout_time_and_phase_encoding_direction",
    "check_dwi_volume",
    "check_dwi_dataset",
]


def generate_index_file(
    b_values_filename: Path,
    image_id: Optional[str] = None,
    output_dir: Optional[Path] = None,
) -> Path:
    """Generate [`image_id`]_index.txt file for FSL eddy command.

    At the moment, all volumes are assumed to be acquired with the
    same parameters. The generate_acq_file function writes a single
    line, and this function writes a vector of ones linking each
    DWI volume to this first line.

    Parameters
    ----------
    b_values_filename : Path
        The path to the b-values file.

    image_id : str, optional
        An optional prefix for the output file name.
        Defaults to None.

    output_dir : Path, optional
        The path to the directory in which the index file
        should be written. If not provided, it will be written
        in the same folder as the provided b values filename.

    Returns
    -------
    index_filename: Path
        The path to output index file. [`image_id`]_index.txt or index.txt file.
    """
    if not b_values_filename.is_file():
        raise FileNotFoundError(f"Unable to find b-values file: {b_values_filename}.")

    b_values = np.loadtxt(b_values_filename)
    index_filename = f"{image_id}_index.txt" if image_id else "index.txt"
    output_dir = output_dir or b_values_filename.parent
    index_filename = output_dir / index_filename
    np.savetxt(index_filename, np.ones(len(b_values)).T)

    return index_filename


def generate_acq_file(
    dwi_filename: Path,
    fsl_phase_encoding_direction: str,
    total_readout_time: str,
    image_id: Optional[str] = None,
    output_dir: Optional[Path] = None,
) -> Path:
    """Generate [`image_id`]_acq.txt file for FSL eddy command.

    Parameters
    ----------
    dwi_filename : str
        Path to the DWI file.

    fsl_phase_encoding_direction : str
        Phase encoding direction from the BIDS specifications in FSL format
        (i.e. x/y/z instead of i/j/k).

    total_readout_time : str
        Total readout time from BIDS specifications.

    image_id : str, optional
        The prefix for the output file. Defaults to None.

    output_dir : Path, optional
        The path to the directory in which the acquisition file
        should be written. If not provided, it will be written
        in the same folder as the provided dwi filename.

    Returns
    -------
    acq_filename : Path
        The path to the acquisition file.
    """
    if fsl_phase_encoding_direction not in ("x", "y", "z", "x-", "y-", "z-"):
        raise RuntimeError(
            f"FSL PhaseEncodingDirection (found value: {fsl_phase_encoding_direction}) "
            f"is unknown, it should be a value in (x, y, z, x-, y-, z-)"
        )
    acq_filename = f"{image_id}_acq.txt" if image_id else "acq.txt"
    output_dir = output_dir or dwi_filename.parent
    acq_filename = output_dir / acq_filename
    basis_vector = _get_phase_basis_vector(fsl_phase_encoding_direction)
    basis_vector.append(float(total_readout_time))
    np.savetxt(acq_filename, np.array([basis_vector]), fmt="%d " * 3 + "%f")

    return acq_filename


def _get_phase_basis_vector(phase: str) -> list:
    """Returns the unit vector corresponding to the given phase."""
    mult = -1 if phase.endswith("-") else 1
    idx = ["x", "y", "z"].index(phase[0])
    result = [0] * 3
    result[idx] = mult

    return result


def compute_average_b0(
    dwi_filename: PathLike,
    b_value_filename: Optional[PathLike] = None,
    b_value_threshold: float = 5.0,
    squeeze: bool = False,
    out_file: Optional[str] = None,
) -> PathLike:
    """Computes the average of the b0 volumes from DWI dataset.

    Parameters
    ----------
    dwi_filename : PathLike
        The path to the DWI files containing the volumes of interest.

    b_value_filename : PathLike, optional
        The path to the b-values file. This will be used to filter the volumes.
        Only volumes with a corresponding b-value smaller than b_value_threshold
        will be kept for computations. If None, all volumes are kept.
        Default is None.

    b_value_threshold : float, optional
        The threshold to apply to b-values to get the volumes which should be kept.
        Default=5.0.

    squeeze : bool, optional
        Whether the returned volume should be squeezed or not.
        If True, it will be a 3D array, if False, a 4D array with a dummy 4th dimension.
        Default=False

    out_file : str, optional
        Name of the output file.
        If None, the output file will be built from the input DWI file base with
        the suffix '_avg_b0.nii[.gz]'.
        Defaults to None.

    Returns
    -------
    out_file : str
        The path to the nifti image file containing the mean of the b0 volumes.

    Raises
    ------
    FileNotFoundError:
        If the DWI input file does not exist.

    ValueError:
        If b_value_threshold < 0.
    """
    from clinica.utils.image import compute_aggregated_volume, get_new_image_like

    dwi_filename = check_file(dwi_filename)
    check_b_value_threshold(b_value_threshold)
    out_file = out_file or add_suffix_to_filename(dwi_filename, "avg_b0")
    volume_filter = (
        get_b0_filter(b_value_filename, b_value_threshold=b_value_threshold)
        if b_value_filename
        else None
    )
    b0 = compute_aggregated_volume(dwi_filename, np.average, volume_filter)
    if not squeeze:
        b0 = b0[..., np.newaxis]
    b0_img = get_new_image_like(dwi_filename, b0)
    b0_img.to_filename(out_file)

    return out_file


def check_file(filename: Union[str, PathLike], resolve: bool = False) -> Path:
    """Check that filename exists and return a Path object."""
    filename = Path(filename).resolve() if resolve else Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"File not found : {filename}.")

    return filename


def add_suffix_to_filename(filename: Path, suffix: str) -> Path:
    """Adds the provided suffix to the given file name."""
    ext = filename.suffix
    if ext == ".gz":
        ext = "." + filename.stem.split(".")[-1] + ext
    return filename.parent / f"{filename.name.replace(ext, '')}_{suffix}{ext}"


def get_b0_filter(
    b_value_filename: PathLike, b_value_threshold: float = 5.0
) -> np.ndarray:
    """Return the index for the volumes where b<=low_bval.

    Parameters
    ----------
    b_value_filename : PathLike
        Path to the bval file.

    b_value_threshold : float, optional
        Defines the threshold for the b0 volumes as all volumes which
        satisfy b-value <= b_value_threshold.
        Defaults to 5.0.

    Returns
    -------
    np.ndarray :
        Index of the b0 volumes.

    Raises
    ------
    FileNotFoundError:
        If b_value_filename file cannot be found.
    """
    b_value_filename = check_file(b_value_filename)
    values = np.loadtxt(b_value_filename)

    return np.where(values <= b_value_threshold)[0]


def check_b_value_threshold(b_value_threshold: float) -> None:
    if b_value_threshold < 0:
        raise ValueError(
            f"b_value_threshold should be >=0. You provided {b_value_threshold}."
        )


def get_readout_time_and_phase_encoding_direction(
    dwi_json_filename: str | Path,
) -> tuple[str, str]:
    """Extract the readout time and phase encoding direction from the DWI JSON file."""
    from clinica.utils.filemanip import extract_metadata_from_json

    [total_readout_time, phase_encoding_direction] = extract_metadata_from_json(
        dwi_json_filename,
        [
            "TotalReadoutTime",
            "PhaseEncodingDirection",
        ],
        handle_missing_keys=_handle_missing_keys_dwi,
    )
    phase_encoding_direction = _bids_dir_to_fsl_dir(phase_encoding_direction)

    return total_readout_time, phase_encoding_direction


def _handle_missing_keys_dwi(data: dict, missing_keys: set[str]) -> dict:
    """Find alternative fields from the bids/sub-X/ses-Y/dwi/sub-X_ses-Y_dwi.json
    file to replace those which were not found in this very json.


    Parameters
    ----------
    data: dict
        Dictionary containing the json data.

    missing_keys: set of str
        Set of keys that are required and were not found in the json.

    Returns
    -------
    dict:
        Contains the values for the requested fields.
    """
    handlers = {
        "TotalReadoutTime": _handle_missing_total_readout_time,
        "PhaseEncodingDirection": _handle_missing_phase_encoding_direction,
    }
    try:
        return {k: handlers[k](data, missing_keys) for k in missing_keys}
    except KeyError:
        raise ValueError(
            f"Could not recover the missing keys {missing_keys} from JSON file."
        )


def _handle_missing_total_readout_time(data: dict, missing_keys: set) -> float:
    """Find an alternative field in the json to replace the TotalReadoutTime.

    Parameters
    ----------
    data: dict
        Dictionary containing the json data.

    missing_keys: set
        Set of keys that are required and were not found in the json.

    Returns
    -------
    float:
        Value for TotalReadoutTime.
    """
    from clinica.utils.exceptions import ClinicaException

    if "EstimatedTotalReadoutTime" in data:
        return data["EstimatedTotalReadoutTime"]
    if "PhaseEncodingSteps" in data and "PixelBandwidth" in data:
        if data["PixelBandwidth"] != 0:
            return data["PhaseEncodingSteps"] / data["PixelBandwidth"]
        raise ValueError("Pixel Bandwidth value is not valid.")
    raise ClinicaException("Could not recover the TotalReadoutTime from JSON file.")


def _handle_missing_phase_encoding_direction(data: dict, missing_keys: set) -> float:
    """Find an alternative field in the json to replace the PhaseEncodingDirection.

    Parameters
    ----------
    data: dict
        Dictionary containing the json data.

    missing_keys: set
        Set of keys that are required and were not found in the json.

    Returns
    -------
    float:
        Value for PhaseEncodingDirection.
    """
    from clinica.utils.exceptions import ClinicaException

    if "PhaseEncodingAxis" in data:
        return data["PhaseEncodingAxis"] + "+"
    raise ClinicaException(
        "Could not recover the PhaseEncodingDirection from JSON file."
    )


def _bids_dir_to_fsl_dir(bids_dir):
    """Converts BIDS PhaseEncodingDirection parameters (i,j,k,i-,j-,k-) to FSL direction (x,y,z,x-,y-,z-)."""
    fsl_dir = bids_dir.lower()
    if "i" not in fsl_dir and "j" not in fsl_dir and "k" not in fsl_dir:
        raise ValueError(
            f"Unknown PhaseEncodingDirection {fsl_dir}: it should be a value in (i, j, k, i-, j-, k-)"
        )

    return (
        fsl_dir.replace("i", "x").replace("j", "y").replace("k", "z").replace("+", "")
    )


def check_dwi_volume(dwi_dataset: DWIDataset) -> None:
    """Checks the consistency of a given DWIDataset.

    More precisely, checks that the number of DWI volumes,
    the number of B-values, and the number of B-vectors are equal.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        The DWI dataset to check.

    Raises
    ------
    ValueError:
        if the number of DWI volumes, the number of B-values,
        and the number of B-vectors are not equal.
    """
    num_b_values = len(np.loadtxt(dwi_dataset.b_values))
    num_b_vectors = np.loadtxt(dwi_dataset.b_vectors).shape[-1]
    num_dwi = nib.load(dwi_dataset.dwi).shape[-1]

    if not (num_b_values == num_b_vectors == num_dwi):
        raise IOError(
            f"Number of DWIs, b-vals and b-vecs mismatch "
            f"(# DWI = {num_dwi}, # B-vec = {num_b_vectors}, #B-val = {num_b_values}) "
        )


def check_dwi_dataset(dwi_dataset: DWIDataset) -> DWIDataset:
    """Checks that provided input files of the DWI dataset exist.

    If they exist, a new dataset object is returned with file paths
    converted to Path objects instead of strings.

    The resulting dataset will contain absolute paths which can be
    required by downstream tools.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        Input DWI dataset to be checked.

    Returns
    -------
    DWIDataset :
        Checked DWI dataset.

    Raises
    ------
    FileNotFoundError
        If one of the files doesn't exist.
    """
    return DWIDataset(*(check_file(f, resolve=True) for f in dwi_dataset))
