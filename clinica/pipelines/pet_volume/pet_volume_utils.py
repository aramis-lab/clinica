def init_input_node(pet_nii):
    import nibabel as nib

    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(pet_nii)

    # Check that the PET file is a 3D volume
    img = nib.load(pet_nii)
    if len(img.shape) == 4:
        error_msg = (
            f"Clinica does not handle 4D volumes for {image_id.replace('_', ' | ')}"
        )
        cprint(error_msg, lvl="error")
        raise NotImplementedError(error_msg)

    # Print begin message
    print_begin_image(image_id)

    return pet_nii


def _check_non_empty_tissue_list(tissues: list) -> None:
    """Check that provided list is non-empty."""
    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )


def _load_tissues(tissues: list):
    """Aggregates the image data contained in the tissue images provided.

    Parameters
    ----------
    tissues : list
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
    import nibabel as nib
    import numpy as np

    from clinica.pipelines.pet_volume.pet_volume_utils import (  # noqa
        _check_non_empty_tissue_list,
    )

    _check_non_empty_tissue_list(tissues)
    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_fdata(dtype="float32").shape)
    data = np.zeros(shape=shape)
    for image in tissues:
        data += nib.load(image).get_fdata(dtype="float32")
    return data, img_0.affine, img_0.header


def create_binary_mask(
    tissues: list,
    threshold: float = 0.3,
) -> str:
    """Create a binary mask Nifti1Image from the list of tissues.

    Tissue images are summed and the result is thresholded with the
    provided `threshold` input.

    Parameters
    ----------
    tissues : list
        List of paths to tissue Nifti1Images. Must be non-empty.

    threshold : float, optional
        Threshold to apply when binarizing the Nifti1Image.
        Default=0.3.

    Returns
    -------
    out_mask : str
        Path to the binary mask Nifti1Image as a string.
    """
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib

    from clinica.pipelines.pet_volume.pet_volume_utils import _load_tissues  # noqa

    data, affine, header = _load_tissues(tissues)
    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + "_brainmask.nii")
    mask = nib.Nifti1Image(data, affine, header=header)
    nib.save(mask, out_mask)
    return out_mask


def apply_binary_mask(image: str, binary_mask: str) -> str:
    """Apply the provided `binary_mask` to the provided `image`.

    Parameters
    ----------
    image : str
        Path to the Nifti1Image to apply the mask on.

    binary_mask : str
        Path to the Nifti1Image containing the mask.

    Returns
    -------
    masked_image_path : str
        Path to the masked Nifti1Image.
    """
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib

    original_image = nib.load(image)
    mask = nib.load(binary_mask)

    data = original_image.get_fdata(dtype="float32") * mask.get_fdata(dtype="float32")

    masked_image_path = join(getcwd(), "masked_" + basename(image))
    masked_image = nib.Nifti1Image(
        data, original_image.affine, header=original_image.header
    )
    nib.save(masked_image, masked_image_path)
    return masked_image_path


def create_pvc_mask(tissues: list) -> str:
    """Create a pvc mask from tissue list.

    Parameters
    ----------
    tissues : list
        List of paths to tissue Nifti1Images. Must be non-empty.

    Returns
    -------
    out_mask : str
        Path to the resulting mask Nifti1Image.
    """
    from os import getcwd
    from os.path import join

    import nibabel as nib
    import numpy as np

    from clinica.pipelines.pet_volume.pet_volume_utils import _load_tissues  # noqa

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


def pet_pvc_name(pet_image: str, pvc_method: str) -> str:
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
    >>> pet_pvc_name(
    ...     "/home/bids/sub-01/ses-M00/pet/sub-01_ses-M00_task-rest_trc-av45_pet.nii.gz",
    ...     "RBV"
    ...)
    'pvc-rbv_sub-01_ses-M00_task-rest_trc-av45_pet.nii.gz'

    """
    from os.path import basename

    return "pvc-" + pvc_method.lower() + "_" + basename(pet_image)


def normalize_to_reference(pet_image: str, region_mask: str) -> str:
    """Normalize the provided `pet_image` by dividing by the mean
    value of the region defined by the provided `region_mask`.

    Parameters
    ----------
    pet_image : str
        Path to the Nifti1Image which should be normalized.

    region_mask : str
        Path to the mask to be used to define the region.

    Returns
    -------
    suvr_pet_path : str
        Path to the normalized Nifti1Image.
    """
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib
    import numpy as np

    pet = nib.load(pet_image)
    ref = nib.load(region_mask)

    region = pet.get_fdata(dtype="float32") * ref.get_fdata(dtype="float32")
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    data = pet.get_fdata(dtype="float32") / region_mean

    suvr_pet_path = join(getcwd(), "suvr_" + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def atlas_statistics(in_image: str, in_atlas_list: list) -> list:
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
    from os import getcwd
    from os.path import abspath, join

    from nipype.utils.filemanip import split_filename

    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas

    orig_dir, base, ext = split_filename(str(in_image))
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics_list = []
    for atlas in in_atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(
                    join(getcwd(), base + "_space-" + atlas + "_statistics.tsv")
                )
                statistics_on_atlas(str(in_image), atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break
    return atlas_statistics_list


def get_from_list(in_list, index):
    return in_list[index]
