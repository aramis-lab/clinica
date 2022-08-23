import typing as ty
from os import PathLike, getcwd
from os.path import basename, join
from pathlib import PurePath

import nibabel as nib
import numpy as np
from pydra.mark import annotate, task


@task
@annotate({"return": {"suvr_pet_path": PurePath}})
def normalize_to_reference(pet_image: str, region_mask: str) -> PurePath:
    """Normalize the provided `pet_image` by dividing by the mean
    value of the region defined by the provided `region_mask`.

    Parameters
    ----------
    pet_image : Path to the Nifti1Image which should be normalized.
    region_mask : Path to the mask to be used to define the region.

    Returns
    -------
    suvr_pet_path : Path to the normalized Nifti1Image.
    """
    pet = nib.load(pet_image)
    ref = nib.load(region_mask)

    region = pet.get_fdata(dtype="float32") * ref.get_fdata(dtype="float32")
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    data = pet.get_fdata(dtype="float32") / region_mean

    suvr_pet_path = join(getcwd(), "suvr_" + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def _check_non_empty_tissue_list(tissues: ty.List):
    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )


def _load_tissues(tissues: ty.List) -> np.ndarray:
    _check_non_empty_tissue_list(tissues)
    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_fdata(dtype="float32").shape)
    data = np.zeros(shape=shape)
    for image in tissues:
        data += nib.load(image).get_fdata(dtype="float32")
    return data


@task
@annotate({"return": {"out_mask": PurePath}})
def create_binary_mask(tissues: ty.List, threshold: float = 0.3) -> PurePath:
    """Create a binary mask Nifti1Image from the list of tissues.

    Tissue images are summed and the result is thresholded with the
    provided `threshold` input.

    Parameters
    ----------
    tissues : List of paths to tissue Nifti1Images. Must be non-empty.
    threshold : Threshold to apply when binarizing the Nifti1Image.

    Returns
    -------
    out_mask : Path to the binary mask Nifti1Image.
    """
    data = _load_tissues(tissues)
    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + "_brainmask.nii")
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


@task
@annotate({"return": {"masked_image_path": PurePath}})
def apply_binary_mask(image: str, binary_mask: str) -> PurePath:
    """Apply the provided `binary_mask` to the provided `image`.

    Parameters
    ----------
    image : Path to the Nifti1Image to apply the mask on.
    binary_mask : Path to the Nifti1Image containing the mask.

    Returns
    -------
    masked_image_path : Path to the masked Nifti1Image.
    """
    original_image = nib.load(image)
    mask = nib.load(binary_mask)

    data = original_image.get_fdata(dtype="float32") * mask.get_fdata(dtype="float32")

    masked_image_path = join(getcwd(), "masked_" + basename(image))
    masked_image = nib.Nifti1Image(
        data, original_image.affine, header=original_image.header
    )
    nib.save(masked_image, masked_image_path)
    return masked_image_path


@task
@annotate({"return": {"atlas_statistics_list": ty.List}})
def atlas_statistics(in_image: str, in_atlas_list: ty.List) -> ty.List:
    """Generate regional measure from atlas_list in TSV files.

    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Parameters
    ----------
    in_image : A Nifti image
    in_atlas_list : List of names of atlas to be applied

    Returns
    -------
    atlas_statistics_list : List of paths to TSV files
    """
    from os.path import abspath, join

    from nipype.utils.filemanip import split_filename

    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas

    orig_dir, base, ext = split_filename(in_image)
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics_list = []
    for atlas in in_atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(
                    join(getcwd(), base + "_space-" + atlas + "_statistics.tsv")
                )
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break
    return atlas_statistics_list


@task
@annotate({"return": {"out_mask": PurePath}})
def create_pvc_mask(tissues: ty.List) -> PurePath:
    """Create a pvc mask from tissue list.

    Parameters
    ----------
    tissues : List of paths to tissue Nifti1Images. Must be non-empty.

    Returns
    -------
    out_mask : Path to the resulting mask Nifti1Image.
    """
    background = _load_tissues(tissues)
    shape = background.shape
    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)
    for i, tissue in enumerate(tissues):
        image = nib.load(tissue)
        data[..., i] = np.array(image.get_fdata(dtype="float32"))
    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)
    out_mask = join(getcwd(), "pvc_mask.nii")
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask
