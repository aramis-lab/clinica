"""This module contains utilities for DWI handling."""


def count_b0s(in_bval: str, low_bval: float = 5.0) -> int:
    """Count the number of volumes where b<=low_bval.

    Parameters
    ----------
    in_bval : str
        Path to the bval file.

    low_bval : int, optional
        Define the threshold for the b0 volumes as all volumes will
        satisfy bval <= lowbval. Defaults to 5.0.

    Returns
    -------
    num_b0s : int
        Number of b0 volumes.
    """
    return len(get_b0_filter(in_bval, low_bval=low_bval))


def get_b0_filter(in_bval: str, low_bval: float = 5.0):
    """Return the index for the volumes where b<=low_bval.

    Parameters
    ----------
    in_bval : str
        Path to the bval file.

    low_bval : int, optional
        Define the threshold for the b0 volumes as all volumes will
        satisfy bval <= lowbval. Defaults to 5.0.

    Returns
    -------
    np.array :
        Index of the b0 volumes.

    Raises
    ------
    FileNotFoundError:
        If in_bval file cannot be found.
    """
    import os

    import numpy as np

    if not os.path.isfile(in_bval):
        raise FileNotFoundError(f"Cannot find bval file : {in_bval}.")
    values = np.loadtxt(in_bval)
    return np.where(values <= low_bval)[0]


def compute_average_b0(
    in_dwi: str, in_bval: str = None, low_bval: float = 5.0, out_file: str = None
) -> str:
    """Compute the average of the b0 volumes from DWI dataset.

    Parameters
    ----------
    in_dwi : str
        The path to the DWI files containing the volumes of interest.

    in_bval : str, optional
        The path to the bval file. This will be used to filter the volumes.
        Only volumes for which bval is less than low_bval will be kept for computations.
        If None, all volumes are kept. Default is None.

    low_bval : float, optional
        The threshold to apply to bvalues to get the volumes which should be kept.
        Default=5.0.

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

    import numpy as np

    from clinica.utils.image import compute_aggregated_volume, get_new_image_like

    in_dwi = Path(in_dwi)

    if not in_dwi.exists():
        raise FileNotFoundError(f"DWI file not found : {in_dwi}.")

    if low_bval < 0:
        raise ValueError(f"low_bval should be >=0. You provided : {low_bval}.")

    if not out_file:
        ext = in_dwi.suffix
        if ext == ".gz":
            ext = "." + in_dwi.stem.split(".")[-1] + ext
        out_file = in_dwi.parent / f"{in_dwi.name.rstrip(ext)}_avg_b0{ext}"

    volume_filter = get_b0_filter(in_bval, low_bval=low_bval) if in_bval else None
    b0 = compute_aggregated_volume(in_dwi, np.average, volume_filter)
    b0_img = get_new_image_like(in_dwi, b0)
    b0_img.to_filename(out_file)

    return out_file


def b0_dwi_split(in_dwi: str, in_bval: str, in_bvec: str, low_bval: float = 5.0):
    """Split DWI dataset.

    Split the DWI volumes into two datasets :
     - the first dataset contains the set of b<=low_bval volumes.
     - the second dataset contains the set of DWI volumes.

    Args:
        in_dwi (str): DWI dataset.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.
        low_bval (float, optional): Define the b0 volumes as all volume bval <= lowbval. Defaults to 5.0.

    Returns:
        out_b0 (str): The set of b<=low_bval volumes.
        out_dwi (str): Output. The set of b>low_bval volumes.
        out_bvals (str): The b-values corresponding to the out_dwi.
        out_bvecs (str): The b-vecs corresponding to the out_dwi.
    """
    import os
    import warnings

    import nibabel as nib
    import numpy as np

    assert os.path.isfile(in_dwi)
    assert os.path.isfile(in_bval)
    assert os.path.isfile(in_bvec)
    assert low_bval >= 0

    im = nib.load(in_dwi)
    data = im.get_fdata(dtype="float32")
    hdr = im.get_header().copy()
    bvals = np.loadtxt(in_bval)
    bvecs = np.loadtxt(in_bvec)

    if bvals.shape[0] == bvecs.shape[0]:
        warnings.warn(
            "Warning: The b-vectors file should be column-wise. The b-vectors will be transposed",
            UserWarning,
        )
        bvecs = bvecs.T

    lowbs = np.where(bvals <= low_bval)[0]

    fname_b0, ext_b0 = os.path.splitext(os.path.basename(in_dwi))
    if ext_b0 == ".gz":
        fname_b0, ext2 = os.path.splitext(fname_b0)
        ext_b0 = ext2 + ext_b0
    out_b0 = os.path.abspath(f"{fname_b0}_b0{ext_b0}")
    # out_b0 = op.abspath('b0.nii.gz')
    b0data = data[..., lowbs]
    hdr.set_data_shape(b0data.shape)
    nib.Nifti1Image(b0data, im.get_affine(), hdr).to_filename(out_b0)

    dwi_bvals = np.where(bvals > low_bval)[0]
    out_dwi = os.path.abspath("dwi.nii.gz")
    dwi_data = data[..., dwi_bvals]
    hdr.set_data_shape(dwi_data.shape)
    nib.Nifti1Image(dwi_data, im.get_affine(), hdr).to_filename(out_dwi)

    bvals_dwi = bvals[dwi_bvals]
    out_bvals = os.path.abspath("bvals")
    np.savetxt(out_bvals, bvals_dwi, fmt="%d", delimiter=" ")

    bvecs_dwi = np.array(
        [
            bvecs[0][dwi_bvals].tolist(),
            bvecs[1][dwi_bvals].tolist(),
            bvecs[2][dwi_bvals].tolist(),
        ]
    )
    out_bvecs = os.path.abspath("bvecs")
    np.savetxt(out_bvecs, bvecs_dwi, fmt="%10.5f", delimiter=" ")

    return out_b0, out_dwi, out_bvals, out_bvecs


def insert_b0_into_dwi(in_b0, in_dwi, in_bval, in_bvec):
    """Insert a b0 volume into the dwi dataset as the first volume and update the bvals and bvecs files.

    Args:
        in_b0 (str): One b=0 volume (could be the average of a b0 dataset).
        in_dwi (str): The set of DWI volumes.
        in_bval (str): File describing the b-values of the DWI dataset.
        in_bvec (str): File describing the directions of the DWI dataset.

    Returns:
        out_dwi (str): Diffusion dataset : b0 volume + dwi volumes.
        out_bval (str): B-values update.
        out_bvec (str): Directions of diffusion update.
    """
    import os

    import numpy as np

    from clinica.utils.image import merge_volumes_time_dimension

    assert os.path.isfile(in_b0)
    assert os.path.isfile(in_dwi)
    assert os.path.isfile(in_bval)
    assert os.path.isfile(in_bvec)

    out_dwi = merge_volumes_time_dimension(in_b0, in_dwi)

    lst = np.loadtxt(in_bval).tolist()
    lst.insert(0, 0)
    out_bvals = os.path.abspath("bvals")
    np.savetxt(out_bvals, np.matrix(lst), fmt="%d", delimiter=" ")

    bvecs = np.loadtxt(in_bvec)
    bvecs_0 = bvecs[0].tolist()
    bvecs_0.insert(0, 0.0)
    bvecs_1 = bvecs[1].tolist()
    bvecs_1.insert(0, 0.0)
    bvecs_2 = bvecs[2].tolist()
    bvecs_2.insert(0, 0.0)
    bvecs_dwi = np.array([bvecs_0, bvecs_1, bvecs_2])
    out_bvecs = os.path.abspath("bvecs")
    np.savetxt(out_bvecs, bvecs_dwi, fmt="%10.5f", delimiter=" ")

    return out_dwi, out_bvals, out_bvecs


def check_dwi_volume(in_dwi, in_bvec, in_bval):
    """Check that # DWI = # B-val = # B-vec.

    Raises
        ValueError if # DWI, # B-val and # B-vec mismatch.
    """
    import nibabel as nib
    import numpy as np

    bvals = np.loadtxt(in_bval)
    num_b_vals = len(bvals)

    bvecs = np.loadtxt(in_bvec)
    _, num_b_vecs = bvecs.shape

    img = nib.load(in_dwi)
    _, _, _, num_dwis = img.shape

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
