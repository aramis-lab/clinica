from os import PathLike
from pathlib import Path
from typing import List, Tuple, Union

import nibabel as nib
import numpy as np
from nibabel import Nifti1Header

__all__ = [
    "apply_binary_mask",
    "compute_atlas_statistics",
    "create_binary_mask",
    "create_pvc_mask",
    "get_from_list",
    "init_input_node",
    "normalize_to_reference",
    "build_pet_pvc_name",
]


def apply_binary_mask(image: Path, binary_mask: Path) -> Path:
    """Apply the provided `binary_mask` to the provided `image`.

    Parameters
    ----------
    image : Path
        The path to the Nifti1Image to apply the mask on.

    binary_mask : Path
        The path to the Nifti1Image containing the mask.

    Returns
    -------
    masked_image_path : Path
        The path to the masked Nifti1Image.
    """
    import nibabel as nib

    original_image = nib.load(image)
    mask = nib.load(binary_mask)
    data = original_image.get_fdata(dtype="float32") * mask.get_fdata(dtype="float32")
    masked_image_path = Path.cwd() / f"masked_{image.name}"
    masked_image = nib.Nifti1Image(
        data, original_image.affine, header=original_image.header
    )
    nib.save(masked_image, masked_image_path)

    return masked_image_path


def compute_atlas_statistics(image: Path, atlas_names: List[str]) -> List[Path]:
    """Generate regional measure from atlas_list in TSV files.

    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Parameters
    ----------
    image : Path
        The path to the Nifti image.

    atlas_names : List of str
        List of names of atlas to be applied.

    Returns
    -------
    atlas_statistics : List of path
        List of paths to TSV files.
    """
    from clinica.utils.filemanip import get_filename_no_ext
    from clinica.utils.statistics import statistics_on_atlas

    atlas_statistics_list = []
    for atlas in atlas_names:
        out_atlas_statistics = (
            Path.cwd() / f"{get_filename_no_ext(image)}_space-{atlas}_statistics.tsv"
        )
        statistics_on_atlas(image, atlas, out_atlas_statistics)
        atlas_statistics_list.append(out_atlas_statistics)
        break  # Why is there a break here ????

    return atlas_statistics_list


def create_binary_mask(
    tissues: List[PathLike],
    threshold: float = 0.3,
) -> Path:
    """Create a binary mask Nifti1Image from the list of tissues.

    Tissue images are summed and the result is thresholded with the
    provided `threshold` input.

    Parameters
    ----------
    tissues : list of PathLike
        List of paths to tissue Nifti1Images. Must be non-empty.

    threshold : float, optional
        Threshold to apply when binarizing the Nifti1Image.
        Default=0.3.

    Returns
    -------
    out_mask : Path
        The path to the binary mask Nifti1Image as a string.
    """
    from clinica.utils.filemanip import get_filename_no_ext

    data, affine, header = _aggregate_tissue_images(tissues)
    data = (data > threshold) * 1.0
    out_mask = Path.cwd() / f"{get_filename_no_ext(tissues[0])}_brainmask.nii"
    mask = nib.Nifti1Image(data, affine, header=header)
    nib.save(mask, out_mask)

    return out_mask


def _aggregate_tissue_images(
    tissues: List[PathLike],
) -> Tuple[np.ndarray, np.ndarray, Nifti1Header]:
    """Aggregates the image data contained in the tissue images provided.

    Parameters
    ----------
    tissues : list of PathLike
        List of tissue images to aggregate.

    Returns
    -------
    data : np.ndarray
        Aggregated data.

    affine : np.ndarray
        Affine of the first image, acting as the affine
        of the aggregated image under the assumption that
        all images have the same affine.

    header : Nifti1Header
        Header of the aggregated image.
    """
    _check_non_empty_tissue_list(tissues)
    first_image = nib.load(tissues[0])
    data = np.zeros(shape=list(first_image.get_fdata(dtype="float32").shape))
    for image in tissues:
        data += nib.load(image).get_fdata(dtype="float32")
    return data, first_image.affine, first_image.header


def _check_non_empty_tissue_list(tissues: List[PathLike]) -> None:
    """Check that provided list is non-empty."""
    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )


def create_pvc_mask(tissues: List[PathLike]) -> Path:
    """Create a pvc mask from tissue list.

    Parameters
    ----------
    tissues : list of PathLike
        List of paths to tissue Nifti1Images. Must be non-empty.

    Returns
    -------
    out_mask : Path
        The path to the resulting mask Nifti1Image.
    """
    background, affine, header = _aggregate_tissue_images(tissues)
    shape = background.shape
    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)
    for i, tissue in enumerate(tissues):
        image = nib.load(tissue)
        data[..., i] = np.array(image.get_fdata(dtype="float32"))
    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)
    out_mask = Path.cwd() / "pvc_mask.nii"
    mask = nib.Nifti1Image(data, affine, header=header)
    nib.save(mask, out_mask)

    return out_mask


def get_from_list(in_list: list, index: int):
    return in_list[index]


def init_input_node(pet_nii: str) -> str:
    from clinica.utils.filemanip import get_subject_id, load_volume
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    image_id = get_subject_id(pet_nii)
    try:
        load_volume(pet_nii)
    except ValueError as e:
        error_msg = f"Clinica could not load volumes for {image_id.replace('_', ' | ')}. {str(e)}"
        cprint(error_msg, lvl="error")
        raise ValueError(error_msg)
    print_begin_image(image_id)

    return pet_nii


def normalize_to_reference(pet_image: Path, region_mask: Path) -> Path:
    """Normalize the provided `pet_image` by dividing by the mean
    value of the region defined by the provided `region_mask`.

    Parameters
    ----------
    pet_image : Path
        The path to the Nifti1Image which should be normalized.

    region_mask : Path
        The path to the mask to be used to define the region.

    Returns
    -------
    suvr_pet_path : Path
        The path to the normalized Nifti1Image.
    """
    pet = nib.load(pet_image)
    ref = nib.load(region_mask)
    region = pet.get_fdata(dtype="float32") * ref.get_fdata(dtype="float32")
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))
    data = pet.get_fdata(dtype="float32") / region_mean
    suvr_pet_path = Path.cwd() / f"suvr_{pet_image.name}"
    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def build_pet_pvc_name(pet_image: Union[str, PathLike], pvc_method: str) -> str:
    """Build the name for the PET PVC interface.

    Parameters
    ----------
    pet_image : str
        Path to the PET scan file.

    pvc_method : str
        Name of the PVC method. This will be concatenated
        with the `pet_image` filename.

    Returns
    -------
    pet_pvc_path : str
        Name for the PET PVC interface.

    Examples
    --------
    >>> build_pet_pvc_name(
    ...     "/home/bids/sub-01/ses-M00/pet/sub-01_ses-M00_task-rest_trc-av45_pet.nii.gz",
    ...     "RBV"
    ...)
    'pvc-rbv_sub-01_ses-M00_task-rest_trc-av45_pet.nii.gz'

    """
    pet_image = Path(pet_image)

    return f"pvc-{pvc_method.lower()}_{pet_image.name}"
