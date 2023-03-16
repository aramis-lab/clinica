"""This module contains utilities for DWI handling."""
from collections import namedtuple
from os import PathLike
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

DWIDataset = namedtuple("DWIDataset", "dwi b_values b_vectors")


def count_b0s(in_bval: PathLike, low_bval: float = 5.0) -> int:
    """Counts the number of volumes where b<=low_bval.

    Parameters
    ----------
    in_bval : PathLike
        Path to the bval file.

    low_bval : int, optional
        Defines the threshold for the b0 volumes as all volumes will
        satisfy bval <= lowbval. Defaults to 5.0.

    Returns
    -------
    num_b0s : int
        Number of b0 volumes.
    """
    return len(get_b0_filter(in_bval, low_bval=low_bval))


def get_b0_filter(in_bval: PathLike, low_bval: float = 5.0) -> np.ndarray:
    """Return the index for the volumes where b<=low_bval.

    Parameters
    ----------
    in_bval : PathLike
        Path to the bval file.

    low_bval : int, optional
        Defines the threshold for the b0 volumes as all volumes will
        satisfy bval <= lowbval. Defaults to 5.0.

    Returns
    -------
    np.ndarray :
        Index of the b0 volumes.

    Raises
    ------
    FileNotFoundError:
        If in_bval file cannot be found.
    """
    import os

    if not os.path.isfile(in_bval):
        raise FileNotFoundError(f"Cannot find bval file : {in_bval}.")
    values = np.loadtxt(in_bval)
    return np.where(values <= low_bval)[0]


def compute_average_b0(
    in_dwi: PathLike,
    in_bval: Optional[PathLike] = None,
    low_bval: float = 5.0,
    squeeze: bool = False,
    out_file: Optional[str] = None,
) -> PathLike:
    """Computes the average of the b0 volumes from DWI dataset.

    Parameters
    ----------
    in_dwi : PathLike
        The path to the DWI files containing the volumes of interest.

    in_bval : PathLike, optional
        The path to the bval file. This will be used to filter the volumes.
        Only volumes for which bval is less than low_bval will be kept for computations.
        If None, all volumes are kept. Default is None.

    low_bval : float, optional
        The threshold to apply to bvalues to get the volumes which should be kept.
        Default=5.0.

    squeeze : bool, optional
        Whether the returned volume should be squeezed or not.
        If True, it will be a 3D array, if False, a 4D array with a dummy 4th dimension.
        Default=False

    out_file : str, optional
        Name of the output file.
        If None, the output file will be built from the input DWI file base with
        the suffix '_avg_b0.nii[.gz]'. Defaults to None.

    Returns
    -------
    out_file : str
        The path to the nifti image file containing the mean of the b0 volumes.

    Raises
    ------
    FileNotFoundError:
        If the DWI input file does not exist.

    ValueError:
        If low_bval < 0.
    """
    from pathlib import Path

    from clinica.utils.image import compute_aggregated_volume, get_new_image_like

    in_dwi = Path(in_dwi)

    if not in_dwi.exists():
        raise FileNotFoundError(f"DWI file not found : {in_dwi}.")

    if low_bval < 0:
        raise ValueError(f"low_bval should be >=0. You provided : {low_bval}.")

    out_file = out_file or add_suffix_to_filename(in_dwi, "avg_b0")
    volume_filter = get_b0_filter(in_bval, low_bval=low_bval) if in_bval else None
    b0 = compute_aggregated_volume(in_dwi, np.average, volume_filter)
    if not squeeze:
        b0 = b0[..., np.newaxis]
    b0_img = get_new_image_like(in_dwi, b0)
    b0_img.to_filename(out_file)

    return out_file


def check_dwi_dataset(dwi_dataset: DWIDataset) -> DWIDataset:
    """Check that provided input files of the DWI dataset exist.

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
    file_paths = []
    for input_file in dwi_dataset:
        input_file_path = Path(input_file).resolve()
        if not input_file_path.exists():
            raise FileNotFoundError(f"File {input_file_path} could not be found.")
        file_paths.append(input_file_path)

    return DWIDataset(*file_paths)


def remove_entity_from_filename(filename: Path, entity: str) -> Path:
    """Removes the provided entity from the given file name."""
    return Path(filename.parent / filename.name.replace(f"_{entity}", ""))


def add_suffix_to_filename(filename: Path, suffix: str) -> Path:
    """Adds the provided suffix to the given file name."""
    ext = filename.suffix
    if ext == ".gz":
        ext = "." + filename.stem.split(".")[-1] + ext
    return filename.parent / f"{filename.name.replace(ext, '')}_{suffix}{ext}"


def b0_dwi_split(
    dwi_dataset: DWIDataset, low_bval: float = 5.0
) -> Tuple[DWIDataset, DWIDataset]:
    """Splits the DWI dataset.

    Split the DWI volumes into two datasets :
     - the first dataset is relative to volumes having a b-value <= low_bval.
     - the second dataset is relative to volumes having a b-value > low_bval.

    The function writes 6 files (3 files for each dataset), and returns a
    length 2 tuple containing the two DWI datasets.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        DWI image dataset to split.

    low_bval : float, optional
        Defines the b0 volumes as all volumes bval <= lowbval.
        Defaults to 5.0.

    Returns
    -------
    small_b_dataset : DWIDataset
        DWI dataset for volumes for which b<=low_bval.

    large_b_dataset : DWIDataset
        DWI dataset for volumes for which b>low_bval.

    Raises
    ------
    ValueError:
        If low_bval < 0.
    """
    if low_bval < 0:
        raise ValueError(f"low_bval should be >=0. You provided {low_bval}.")

    dwi_dataset = check_dwi_dataset(dwi_dataset)
    b_values, b_vectors = _check_b_values_and_b_vectors(dwi_dataset)

    small_b_filter = get_b0_filter(dwi_dataset.b_values, low_bval=low_bval)
    large_b_filter = np.array(
        [i for i in range(len(b_values)) if i not in small_b_filter]
    )

    return (
        _build_dwi_dataset_from_filter(dwi_dataset, "small_b", small_b_filter),
        _build_dwi_dataset_from_filter(dwi_dataset, "large_b", large_b_filter),
    )


def _check_b_values_and_b_vectors(
    dwi_dataset: DWIDataset,
) -> Tuple[np.ndarray, np.ndarray]:
    """Opens the b-values and b-vectors files and transpose the b-vectors if needed."""
    import warnings

    b_values = np.loadtxt(dwi_dataset.b_values)
    b_vectors = np.loadtxt(dwi_dataset.b_vectors)
    if b_values.shape[0] == b_vectors.shape[0]:
        warnings.warn(
            "Warning: The b-vectors file should be column-wise. The b-vectors will be transposed",
            UserWarning,
        )
        b_vectors = b_vectors.T

    return b_values, b_vectors


def _build_dwi_dataset_from_filter(
    dwi_dataset: DWIDataset, filter_name: str, filter_array: np.ndarray
) -> DWIDataset:
    """Builds a new DWI dataset from a given DWI dataset and a filter.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        The DWI dataset to filter.

    filter_name : str
        The name of the filter. This will be used to build the
        file names associated with the new dataset.

    filter_array : np.ndarray
        1D array of indices to filter the DWI dataset.

    Returns
    -------
    DWIDataset :
        The new filtered DWI dataset.
    """
    return check_dwi_dataset(
        DWIDataset(
            dwi=_filter_dwi(dwi_dataset, filter_name, filter_array),
            b_values=_filter_b_values(dwi_dataset, filter_name, filter_array),
            b_vectors=_filter_b_vectors(dwi_dataset, filter_name, filter_array),
        )
    )


def _filter_dwi(
    dwi_dataset: DWIDataset, filter_name: str, filter_array: np.ndarray
) -> Path:
    """Filter the dwi component of the provided DWI dataset."""
    from clinica.utils.image import compute_aggregated_volume, get_new_image_like

    dwi_filename = add_suffix_to_filename(dwi_dataset.dwi, filter_name)
    data = compute_aggregated_volume(
        dwi_dataset.dwi, aggregator=None, volumes_to_keep=filter_array
    )
    img = get_new_image_like(dwi_dataset.dwi, data)
    img.to_filename(dwi_filename)

    return dwi_filename


def _filter_b_values(
    dwi_dataset: DWIDataset, filter_name: str, filter_array: np.ndarray
) -> Path:
    """Filter the b-values component of the provided DWI dataset."""
    b_values, _ = _check_b_values_and_b_vectors(dwi_dataset)
    b_values_filename = add_suffix_to_filename(dwi_dataset.b_values, filter_name)
    np.savetxt(b_values_filename, b_values[filter_array], fmt="%d", delimiter=" ")

    return b_values_filename


def _filter_b_vectors(
    dwi_dataset: DWIDataset, filter_name: str, filter_array: np.ndarray
) -> Path:
    """Filter the b-vectors component of the provided DWI dataset."""
    b_values, b_vectors = _check_b_values_and_b_vectors(dwi_dataset)
    b_vectors_filename = add_suffix_to_filename(dwi_dataset.b_vectors, filter_name)
    np.savetxt(
        b_vectors_filename,
        np.array([b[filter_array] for b in b_vectors]),
        fmt="%10.5f",
        delimiter=" ",
    )

    return b_vectors_filename


def insert_b0_into_dwi(in_b0: PathLike, dwi_dataset: DWIDataset) -> DWIDataset:
    """Inserts a b0 volume into the DWI dataset as the first volume and update the bvals and bvecs files.

    Parameters
    ----------
    in_b0 : PathLike
        Path to image file containing one b=0 volume (could be the average of a b0 dataset).

    dwi_dataset : DWIDataset
        DWI dataset in which to insert the b0 volume.

    Returns
    -------
    DWIDataset :
        The diffusion dataset : b0 volume + dwi volumes.
    """
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    dwi_dataset = check_dwi_dataset(dwi_dataset)
    out_dwi = merge_nifti_images_in_time_dimension(
        (in_b0, dwi_dataset.dwi),
        out_file=add_suffix_to_filename(
            remove_entity_from_filename(dwi_dataset.dwi, "large_b"),
            "merged",
        ),
    )

    bvals = np.loadtxt(dwi_dataset.b_values)
    bvals = np.insert(bvals, 0, 0)
    out_bvals = add_suffix_to_filename(
        remove_entity_from_filename(dwi_dataset.b_values, "large_b"), "merged"
    )
    np.savetxt(out_bvals, bvals, fmt="%d", delimiter=" ")

    bvecs = np.loadtxt(dwi_dataset.b_vectors)
    bvecs = np.insert(bvecs, 0, 0.0, axis=1)
    out_bvecs = add_suffix_to_filename(
        remove_entity_from_filename(dwi_dataset.b_vectors, "large_b"), "merged"
    )
    np.savetxt(out_bvecs, bvecs, fmt="%10.5f", delimiter=" ")

    return DWIDataset(dwi=out_dwi, b_values=out_bvals, b_vectors=out_bvecs)


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
    import nibabel as nib
    import numpy as np

    num_b_vals = len(np.loadtxt(dwi_dataset.b_values))
    num_b_vecs = np.loadtxt(dwi_dataset.b_vectors).shape[-1]
    num_dwis = nib.load(dwi_dataset.dwi).shape[-1]

    if not (num_b_vals == num_b_vecs == num_dwis):
        raise IOError(
            f"Number of DWIs, b-vals and b-vecs mismatch "
            f"(# DWI = {num_dwis}, # B-vec = {num_b_vecs}, #B-val = {num_b_vals}) "
        )


def generate_index_file(b_values_filename: str, image_id: str = None) -> str:
    """Generate [`image_id`]_index.txt file for FSL eddy command.

    At the moment, all volumes are assumed to be acquired with the
    same parameters. The generate_acq_file function writes a single
    line, and this function writes a vector of ones linking each
    DWI volume to this first line.

    Parameters
    ----------
    b_values_filename : str
        Path to the b-values file.

    image_id : str, optional
        Optional prefix for the output file name.
        Defaults to None.

    Returns
    -------
    index_filename: str
        Path to output index file. [`image_id`]_index.txt or index.txt file.
    """
    from pathlib import Path

    import numpy as np

    b_values_filename = Path(b_values_filename)
    if not b_values_filename.is_file():
        raise FileNotFoundError(f"Unable to find b-values file: {b_values_filename}.")

    b_values = np.loadtxt(b_values_filename)
    index_filename = f"{image_id}_index.txt" if image_id else "index.txt"
    index_filename = b_values_filename.parent / index_filename
    np.savetxt(index_filename, np.ones(len(b_values)).T)

    return str(index_filename)


def generate_acq_file(
    dwi_filename: str,
    fsl_phase_encoding_direction: str,
    total_readout_time: str,
    image_id=None,
) -> str:
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
        Optional prefix for the output file. Defaults to None.

    Returns
    -------
    acq_filename : str
        Path to the acq.txt file.
    """
    from pathlib import Path

    import numpy as np

    from clinica.utils.dwi import _get_phase_basis_vector  # noqa

    if fsl_phase_encoding_direction not in ("x", "y", "z", "x-", "y-", "z-"):
        raise RuntimeError(
            f"FSL PhaseEncodingDirection (found value: {fsl_phase_encoding_direction}) "
            f"is unknown, it should be a value in (x, y, z, x-, y-, z-)"
        )
    dwi_filename = Path(dwi_filename)
    acq_filename = f"{image_id}_acq.txt" if image_id else "acq.txt"
    acq_filename = dwi_filename.parent / acq_filename
    basis_vector = _get_phase_basis_vector(fsl_phase_encoding_direction)
    basis_vector.append(float(total_readout_time))
    np.savetxt(acq_filename, np.array([basis_vector]), fmt="%d " * 3 + "%f")

    return str(acq_filename)


def _get_phase_basis_vector(phase: str) -> list:
    """Returns the unit vector corresponding to the given phase."""
    mult = -1 if phase.endswith("-") else 1
    idx = ["x", "y", "z"].index(phase[0])
    result = [0] * 3
    result[idx] = mult

    return result


def bids_dir_to_fsl_dir(bids_dir):
    """Converts BIDS PhaseEncodingDirection parameters (i,j,k,i-,j-,k-) to FSL direction (x,y,z,x-,y-,z-)."""
    fsl_dir = bids_dir.lower()
    if "i" not in fsl_dir and "j" not in fsl_dir and "k" not in fsl_dir:
        raise ValueError(
            f"Unknown PhaseEncodingDirection {fsl_dir}: it should be a value in (i, j, k, i-, j-, k-)"
        )

    return fsl_dir.replace("i", "x").replace("j", "y").replace("k", "z")


def extract_bids_identifier_from_filename(dwi_filename: str) -> str:
    """Extract BIDS identifier from CAPS filename.

    Parameters
    ----------
    dwi_filename : str
        DWI file name for which to extract the bids identifier.

    Returns
    -------
    str :
        The corresponding BIDS identifier.

    Examples
    --------
    >>> extract_bids_identifier_from_filename("sub-01_ses-M000_dwi_space-b0_preproc.bval")
    'sub-01_ses-M000_dwi'
    >>> extract_bids_identifier_from_filename("sub-01_ses-M000_dwi.bvec")
    'sub-01_ses-M000_dwi'
    >>> extract_bids_identifier_from_filename("foo/bar/sub-01_ses-M000_dwi_baz.foo.bar")
    'sub-01_ses-M000_dwi'
    """
    import re

    m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+).*_dwi", dwi_filename)
    if not m:
        raise ValueError(
            f"Could not extract the BIDS identifier from the DWI input filename {dwi_filename}."
        )

    return m.group(0)


def rename_files(in_caps_dwi: str, mapping: dict) -> tuple:
    """Rename files provided.

    The new files are symbolic links to old files.
    For this reason, the old files still exists after renaming.

    Parameters
    ----------
    in_caps_dwi : str
        A DWI file from the CAPS folder.
        This is used only to extract the BIDS identifier.

    mapping : dict
        Mapping between original file names and suffixes for
        new file names.

    Returns
    -------
    tuple :
        New file names.
    """
    import os

    from nipype.interfaces.utility import Rename
    from nipype.utils.filemanip import split_filename

    bids_id = extract_bids_identifier_from_filename(in_caps_dwi)
    renamed_files = []
    for original_file, suffix in mapping.items():
        base_dir, _, _ = split_filename(original_file)
        rename = Rename()
        rename.inputs.in_file = original_file
        rename.inputs.format_string = os.path.join(base_dir, f"{bids_id}{suffix}")
        renamed_files.append(rename.run().outputs.out_file)

    return tuple(renamed_files)
