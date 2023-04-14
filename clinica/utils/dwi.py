"""This module contains utilities for DWI handling."""
import functools
from collections import namedtuple
from os import PathLike
from pathlib import Path
from typing import Optional, Tuple, Union

import numpy as np

DWIDataset = namedtuple("DWIDataset", "dwi b_values b_vectors")


def count_b0s(b_value_filename: PathLike, b_value_threshold: float = 5.0) -> int:
    """Counts the number of volumes where b<=low_bval.

    Parameters
    ----------
    b_value_filename : PathLike
        Path to the b-value file.

    b_value_threshold : float, optional
        Defines the threshold for the b0 volumes as all volumes which
        satisfy b-value <= b_value_threshold.
        Defaults to 5.0.

    Returns
    -------
    num_b0s : int
        Number of b0 volumes.
    """
    return len(get_b0_filter(b_value_filename, b_value_threshold=b_value_threshold))


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
    b_value_filename = _check_file(b_value_filename)
    values = np.loadtxt(b_value_filename)

    return np.where(values <= b_value_threshold)[0]


def _check_file(filename: Union[str, PathLike], resolve: bool = False) -> Path:
    """Check that filename exists and return a Path object."""
    filename = Path(filename).resolve() if resolve else Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"File not found : {filename}.")

    return filename


def compute_average_b0_task(
    dwi_filename: str,
    b_value_filename: str,
    b_value_threshold: float = 5.0,
    squeeze: bool = False,
    out_file: str = None,
) -> str:
    """Nipype task for compute_average_b0."""
    from clinica.utils.dwi import _check_file, compute_average_b0  # noqa

    return str(
        compute_average_b0(
            _check_file(dwi_filename),
            _check_file(b_value_filename),
            b_value_threshold,
            squeeze,
            out_file,
        )
    )


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

    dwi_filename = _check_file(dwi_filename)
    _check_b_value_threshold(b_value_threshold)
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


def _check_b_value_threshold(b_value_threshold: float) -> None:
    if b_value_threshold < 0:
        raise ValueError(
            f"b_value_threshold should be >=0. You provided {b_value_threshold}."
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
    return DWIDataset(*(_check_file(f, resolve=True) for f in dwi_dataset))


def remove_entity_from_filename(filename: Path, entity: str) -> Path:
    """Removes the provided entity from the given file name.

    Parameters
    ----------
    filename : Path
        The path object on which to operate removal of entity.

    entity : str
        The entity which should be removed.

    Returns
    -------
    Path :
        The new Path object without the entity.
    """
    return Path(filename.parent / filename.name.replace(f"_{entity}", ""))


def add_suffix_to_filename(filename: Path, suffix: str) -> Path:
    """Adds the provided suffix to the given file name."""
    ext = filename.suffix
    if ext == ".gz":
        ext = "." + filename.stem.split(".")[-1] + ext
    return filename.parent / f"{filename.name.replace(ext, '')}_{suffix}{ext}"


def split_dwi_dataset_with_b_values(
    dwi_dataset: DWIDataset,
    b_value_threshold: float = 5.0,
) -> Tuple[DWIDataset, DWIDataset]:
    """Splits the DWI dataset in two through the B-values.

    Splits the DWI volumes into two datasets :
     - the first dataset is relative to volumes having a b-value <= b_value_threshold.
     - the second dataset is relative to volumes having a b-value > b_value_threshold.

    The function writes 6 files (3 files for each dataset), and returns a
    length 2 tuple containing the two DWI datasets.

    Parameters
    ----------
    dwi_dataset : DWIDataset
        DWI image dataset to split.

    b_value_threshold : float, optional
        Defines the b0 volumes as all volumes b-value <= b_value_threshold.
        Defaults to 5.0.

    Returns
    -------
    small_b_dataset : DWIDataset
        DWI dataset for volumes for which b-value <= b_value_threshold.

    large_b_dataset : DWIDataset
        DWI dataset for volumes for which b-value > b_value_threshold.

    Raises
    ------
    ValueError:
        If b_value_threshold < 0.
    """
    _check_b_value_threshold(b_value_threshold)
    dwi_dataset = check_dwi_dataset(dwi_dataset)
    b_values, _ = _check_b_values_and_b_vectors(dwi_dataset)
    small_b_filter = get_b0_filter(
        dwi_dataset.b_values, b_value_threshold=b_value_threshold
    )
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
    """Filters the dwi component of the provided DWI dataset."""
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
    """Filters the b-values component of the provided DWI dataset."""
    b_values, _ = _check_b_values_and_b_vectors(dwi_dataset)
    b_values_filename = add_suffix_to_filename(dwi_dataset.b_values, filter_name)
    _write_b_values(b_values_filename, b_values[filter_array])

    return b_values_filename


def _filter_b_vectors(
    dwi_dataset: DWIDataset, filter_name: str, filter_array: np.ndarray
) -> Path:
    """Filters the b-vectors component of the provided DWI dataset."""
    _, b_vectors = _check_b_values_and_b_vectors(dwi_dataset)
    b_vectors_filename = add_suffix_to_filename(dwi_dataset.b_vectors, filter_name)
    _write_b_vectors(b_vectors_filename, np.array([b[filter_array] for b in b_vectors]))

    return b_vectors_filename


def _write_numpy(filename: Path, data: np.ndarray, fmt: str, delimiter: str) -> None:
    """Writes the provided array to the provided filename using the provided formatting."""
    import numpy as np

    np.savetxt(filename, data, fmt=fmt, delimiter=delimiter)


_write_b_vectors = functools.partial(_write_numpy, fmt="%10.5f", delimiter=" ")
_write_b_values = functools.partial(_write_numpy, fmt="%d", delimiter=" ")


def insert_b0_into_dwi(b0_filename: PathLike, dwi_dataset: DWIDataset) -> DWIDataset:
    """Inserts a b0 volume into the DWI dataset as the first volume.

    Also updates the b-values and b-vectors files accordingly.

    Parameters
    ----------
    b0_filename : PathLike
        Path to image file containing one b=0 volume (could be the average of a b0 dataset).

    dwi_dataset : DWIDataset
        DWI dataset in which to insert the b0 volume.

    Returns
    -------
    DWIDataset :
        The diffusion dataset : b0 volume + dwi volumes.
    """
    dwi_dataset = check_dwi_dataset(dwi_dataset)
    b0_filename = _check_file(b0_filename)

    return DWIDataset(
        dwi=_insert_b0_into_dwi_image(dwi_dataset.dwi, b0_filename),
        b_values=_insert_b0_into_b_values(dwi_dataset.b_values),
        b_vectors=_insert_b0_into_b_vectors(dwi_dataset.b_vectors),
    )


def _insert_b0_into_dwi_image(dwi_filename: Path, b0_filename: Path) -> Path:
    """Insert the provided B0 volume into the DWI image.

    This insertion is done at index 0 along the 4th / time dimension.
    """
    from clinica.utils.image import merge_nifti_images_in_time_dimension

    output_dwi_filename = merge_nifti_images_in_time_dimension(
        (b0_filename, dwi_filename),
        out_file=add_suffix_to_filename(
            remove_entity_from_filename(dwi_filename, "large_b"),
            "merged",
        ),
    )

    return output_dwi_filename


def _insert_b0_into_b_values(b_values_filename: Path) -> Path:
    """Insert a 0 value at index 0 into the b-values file."""
    b_values = np.loadtxt(b_values_filename)
    b_values = np.insert(b_values, 0, 0)
    output_b_values_filename = add_suffix_to_filename(
        remove_entity_from_filename(b_values_filename, "large_b"), "merged"
    )
    _write_b_values(output_b_values_filename, b_values)

    return output_b_values_filename


def _insert_b0_into_b_vectors(b_vectors_filename: Path) -> Path:
    """Insert a 0 vector into the b-vectors file."""
    b_vectors = np.loadtxt(b_vectors_filename)
    b_vectors = np.insert(b_vectors, 0, 0.0, axis=1)
    output_b_vectors_filename = add_suffix_to_filename(
        remove_entity_from_filename(b_vectors_filename, "large_b"), "merged"
    )
    _write_b_vectors(output_b_vectors_filename, b_vectors)

    return output_b_vectors_filename


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

    num_b_values = len(np.loadtxt(dwi_dataset.b_values))
    num_b_vectors = np.loadtxt(dwi_dataset.b_vectors).shape[-1]
    num_dwi = nib.load(dwi_dataset.dwi).shape[-1]

    if not (num_b_values == num_b_vectors == num_dwi):
        raise IOError(
            f"Number of DWIs, b-vals and b-vecs mismatch "
            f"(# DWI = {num_dwi}, # B-vec = {num_b_vectors}, #B-val = {num_b_values}) "
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
