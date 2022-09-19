import typing as ty
from os import PathLike, getcwd
from os.path import basename, join
from pathlib import PurePath

import nibabel as nib
import numpy as np
from nibabel.nifti1 import Nifti1Header
from pydra.mark import annotate, task


@task
@annotate({"return": {"suvr_pet_path": PurePath}})
def normalize_to_reference(pet_image: PurePath, region_mask: PurePath) -> PurePath:
    """Normalize the provided `pet_image` by dividing by the mean
    value of the region defined by the provided `region_mask`.

    Parameters
    ----------
    pet_image : PurePath
        Path to the Nifti1Image which should be normalized.

    region_mask : PurePath
        Path to the mask to be used to define the region.

    Returns
    -------
    suvr_pet_path : PurePath
        Path to the normalized Nifti1Image.
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


def _check_non_empty_tissue_list(tissues: ty.List[PurePath]) -> None:
    """Check that provided list is non-empty."""
    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )


def _load_tissues(
    tissues: ty.List[PurePath],
) -> ty.Tuple[np.ndarray, np.ndarray, Nifti1Header]:
    """Aggregates the image data contained in the tissue images provided.

    Parameters
    ----------
    tissues : List[PurePath]
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
    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_fdata(dtype="float32").shape)
    data = np.zeros(shape=shape)
    for image in tissues:
        data += nib.load(image).get_fdata(dtype="float32")
    return data, img_0.affine, img_0.header


@task
@annotate({"return": {"out_mask": PurePath}})
def create_binary_mask(
    tissues: ty.List[PurePath],
    threshold: float = 0.3,
) -> PurePath:
    """Create a binary mask Nifti1Image from the list of tissues.

    Tissue images are summed and the result is thresholded with the
    provided `threshold` input.

    Parameters
    ----------
    tissues : List[PurePath]
        List of paths to tissue Nifti1Images. Must be non-empty.

    threshold : float, optional
        Threshold to apply when binarizing the Nifti1Image.
        Default=0.3.

    Returns
    -------
    out_mask : PurePath
        Path to the binary mask Nifti1Image.
    """
    data, affine, header = _load_tissues(tissues)
    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + "_brainmask.nii")
    mask = nib.Nifti1Image(data, affine, header=header)
    nib.save(mask, out_mask)
    return out_mask


@task
@annotate({"return": {"masked_image_path": PurePath}})
def apply_binary_mask(image: PurePath, binary_mask: PurePath) -> PurePath:
    """Apply the provided `binary_mask` to the provided `image`.

    Parameters
    ----------
    image : PurePath
        Path to the Nifti1Image to apply the mask on.

    binary_mask : PurePath
        Path to the Nifti1Image containing the mask.

    Returns
    -------
    masked_image_path : PurePath
        Path to the masked Nifti1Image.
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
@annotate({"return": {"atlas_statistics": ty.List}})
def atlas_statistics(in_image: str, in_atlas_list: ty.List) -> ty.List:
    """Generate regional measure from atlas_list in TSV files.

    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Parameters
    ----------
    in_image : str
        Path to the Nifti image.

    in_atlas_list : List
        List of names of atlas to be applied.

    Returns
    -------
    atlas_statistics : List
        List of paths to TSV files.
    """
    from os.path import abspath, join

    from nipype.utils.filemanip import split_filename

    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas

    orig_dir, base, ext = split_filename(in_image)
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics = []
    for atlas in in_atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(
                    join(getcwd(), base + "_space-" + atlas + "_statistics.tsv")
                )
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics.append(out_atlas_statistics)
                break
    return atlas_statistics


@task
@annotate({"return": {"out_mask": PurePath}})
def create_pvc_mask(tissues: ty.List) -> PurePath:
    """Create a pvc mask from tissue list.

    Parameters
    ----------
    tissues : List
        List of paths to tissue Nifti1Images. Must be non-empty.

    Returns
    -------
    out_mask : PurePath
        Path to the resulting mask Nifti1Image.
    """
    background, affine, header = _load_tissues(tissues)
    shape = background.shape
    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)
    for i, tissue in enumerate(tissues):
        image = nib.load(tissue)
        data[..., i] = np.array(image.get_fdata(dtype="float32"))
    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)
    out_mask = join(getcwd(), "pvc_mask.nii")
    mask = nib.Nifti1Image(data, affine, header=header)
    nib.save(mask, out_mask)
    return out_mask


@task
@annotate({"return": {"pet_pvc_path": str}})
def pet_pvc_name(pet_image: PurePath, pvc_method: str) -> str:
    """Build the name for the PET PVC interface.

    Parameters
    ----------
    pet_image : PurePath
        Path to the PET scan file.

    pvc_method : str
        Name of the PVC method. This will be concatenated
        with the `pet_image` filename.

    Returns
    -------
    pet_pvc_path : str
        Name for the PET PVC interface.
    """
    from os.path import basename

    return "pvc-" + pvc_method.lower() + "_" + basename(pet_image)


@task
@annotate({"return": {"psf_x": int, "psf_y": int, "psf_z": int}})
def get_psf_task(
    pvc_psf_tsv: PurePath,
    pet_filename: PurePath,
    pet_tracer: str,
) -> ty.Tuple[int]:
    """Returns the tuple (psf_x, psf_y, psf_z) for the subject at the
    session currently treated by the workflow.

    .. note::
        This function uses the name of the pet scan file to retrieve
        the subject and session IDs which are needed by the function
        `read_psf_information`.

    Parameters
    ----------
    pvc_psf_tsv : PurePath
        Path to the custom pvc_psf TSV file.
        See docstring of function `read_psf_information` to get more
        details on this file format.

    pet_filename : PurePath
        Path to the PET image for the subject and session treated by
        the workflow. See note above for details.

    pet_tracer : str
        Label of the PET tracer for which to retrieve psf.

    Returns
    -------
    (psf_x, psf_y, psf_z) : Tuple[int]
        The length 3 tuple for the subject session.
    """
    import re

    from clinica.utils.pet import read_psf_information

    m = re.search(r"(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)", str(filename.name))
    if not m:
        raise ValueError(
            f"Could not extract subject and session IDs from filename {filename}"
        )
    sub = m.group(1)
    ses = m.group(2)
    psf_x, psf_y, psf_z = read_psf_information(pvc_psf_tsv, [sub], [ses], pet_tracer)[0]
    return psf_x, psf_y, psf_z
