# Functions used by nipype interface.
import os

from nibabel.nifti1 import Nifti1Image


def init_input_node(pet):
    """Initiate the pipeline."""
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(pet)
    print_begin_image(image_id)
    return pet


def concatenate_transforms(
    pet_to_t1w_transform: str, t1w_to_mni_transform: str
) -> list:
    """Concatenate two input transformation files into a list.

    Parameters
    ----------
    pet_to_t1w_transform : str
        First transformation to apply.

    t1w_to_mni_transform : str
        Second transformation to apply.

    Returns
    -------
    list :
        Both transform files path in a list.
    """
    return [t1w_to_mni_transform, pet_to_t1w_transform]


# Normalize the images based on the reference mask region
def suvr_normalization(
    input_img: str,
    norm_img: str,
    ref_mask: str,
) -> str:
    """Normalize the input image according to the reference region.

    It uses nilearn `resample_to_img` and scipy `trim_mean` functions.
    This function is different than the one in other PET pipelines
    because there is a downsampling step.

    Parameters
    ----------
    input_img : str
        Path to the image to be processed.

    norm_img : str
        Path to the image used to compute the mean of the reference region.

    ref_mask : str
        Path to the mask of the reference region.

    Returns
    -------
    output_img : str
        Path to the normalized nifti image.
    """
    import os

    import nibabel as nib
    import numpy as np
    from nilearn.image import resample_to_img
    from scipy.stats import trim_mean

    pet = nib.load(input_img)
    norm = nib.load(norm_img)
    mask = nib.load(ref_mask)

    # Downsample the pet image used for normalization so we can multiply it with the mask
    ds_img = resample_to_img(norm, mask, interpolation="nearest")

    # Compute the mean of the region
    region = np.multiply(ds_img.get_fdata(), mask.get_fdata(dtype="float32"))
    array_region = np.where(region != 0, region, np.nan).flatten()
    region_mean = trim_mean(array_region[~np.isnan(array_region)], 0.1)

    from clinica.utils.stream import cprint

    cprint(region_mean)

    # Divide the value of the image voxels by the computed mean
    data = pet.get_fdata(dtype="float32") / region_mean

    # Create and save the normalized image
    output_img = os.path.join(
        os.getcwd(),
        os.path.basename(input_img).split(".nii")[0] + "_suvr_normalized.nii.gz",
    )

    normalized_img = nib.Nifti1Image(data, pet.affine, header=pet.header)
    normalized_img.to_filename(output_img)

    return output_img


def crop_nifti(input_img: str, ref_img: str) -> str:
    """Crop input image based on the reference.

    It uses nilearn `resample_to_img` function.

    Parameters
    ----------
    input_img : str
        Path to the input image.

    ref_img : str
        Path to the reference image used for cropping.

    Returns
    -------
    output_img : str
        Path to the cropped image.
    """
    from pathlib import Path

    from nilearn.image import resample_to_img

    basedir = Path.cwd()
    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_img, force_resample=True)
    crop_filename = Path(input_img.split(".nii")[0] + "_cropped.nii.gz")
    output_img = basedir / crop_filename
    crop_img.to_filename(str(output_img))
    return str(output_img)


def rename_into_caps(
    pet_filename_bids: str,
    pet_filename_raw: str,
    transformation_filename_raw: str,
    suvr_reference_region: str,
    uncropped_image: bool,
    pet_filename_in_t1w_raw: str = None,
    output_dir: str = None,
):
    """
    Rename the outputs of the pipelines into CAPS format.

    More precisely, the following files will be renamed:

        - 'pet_filename_raw' will be renamed to 'pet_filename_caps'
        - 'transformation_filename_raw' will be renamed to 'transformation_filename_caps'
        - 'pet_filename_in_t1w_raw' will be renamed (if provided) to 'pet_filename_in_t1w_caps'

    Parameters
    ----------
    pet_filename_bids : str
        Input PET image file from the BIDS dataset.
        This file is used to extract the entities from the source
        file like subject, session, run, tracer...

    pet_filename_raw : str
        Preprocessed PET file as outputted by Nipype.

    transformation_filename_raw : str
        Transformation file from PET to MRI space as outputted by Nipype.

    suvr_reference_region : str
        SUVR mask name for file name output.
        This will be used to derive the prefix of the PET image.

    uncropped_image : bool
        Pipeline argument for image cropping.
        This will be used to derive the prefix of the PET image.

    pet_filename_in_t1w_raw : str, optional
        Intermediate PET in T1w MRI space.
        If not provided, no renaming will be done.

    output_dir : str, optional
        Specify the output folder where the renamed files should
        be written. This is mostly used for testing purposes in order
        to enable this function to write to Pytest's temporary folders.
        When used in the pipeline, this is left as None as other nodes
        are responsible for adding the missing pieces to the output paths
        and writing the files to disk.

    Returns
    -------
    pet_filename_caps : str
        The renamed preprocessed PET file to match CAPS conventions.

    transformation_filename_caps : str
        The transformation file from PET to MRI space renamed to match CAPS conventions.

    pet_filename_in_t1w_caps : str or None
        Intermediate PET in T1w MRI space renamed to match CAPS conventions.
        If 'pet_filename_in_t1w_raw' is None, this will be None.
    """
    import os

    from clinica.pipelines.pet_linear.pet_linear_utils import (  # noqa
        _get_bids_entities_without_suffix,
        _rename_intermediate_pet_in_t1w_space_into_caps,
        _rename_pet_into_caps,
        _rename_transformation_into_caps,
    )

    bids_entities = _get_bids_entities_without_suffix(pet_filename_bids, suffix="pet")
    if output_dir:
        bids_entities = os.path.join(output_dir, bids_entities)
    pet_filename_caps = _rename_pet_into_caps(
        bids_entities, pet_filename_raw, not uncropped_image, suvr_reference_region
    )
    transformation_filename_caps = _rename_transformation_into_caps(
        bids_entities, transformation_filename_raw
    )
    pet_filename_in_t1w_caps = None
    if pet_filename_in_t1w_raw is not None:
        pet_filename_in_t1w_caps = _rename_intermediate_pet_in_t1w_space_into_caps(
            bids_entities, pet_filename_in_t1w_raw
        )

    return pet_filename_caps, transformation_filename_caps, pet_filename_in_t1w_caps


def _get_bids_entities_without_suffix(filename: str, suffix: str) -> str:
    """Return the BIDS entities without the suffix from a BIDS path."""
    from nipype.utils.filemanip import split_filename

    _, stem, _ = split_filename(filename)
    return stem.rstrip(f"_{suffix}")


def _rename_pet_into_caps(
    entities: str, filename: str, cropped: bool, suvr_reference_region: str
) -> str:
    """Rename into CAPS PET."""
    return _rename(
        filename, entities, _get_pet_bids_components(cropped, suvr_reference_region)
    )


def _rename_transformation_into_caps(entities: str, filename: str) -> str:
    """Rename into CAPS transformation file."""
    return _rename(filename, entities, "_space-T1w_rigid.mat")


def _rename_intermediate_pet_in_t1w_space_into_caps(
    entities: str, filename: str
) -> str:
    """Rename intermediate PET in T1w MRI space."""
    return _rename(filename, entities, "_space-T1w_pet.nii.gz")


def _rename(filename: str, entities: str, suffix: str):
    """Rename 'filename' into '{entities}{suffix}'."""
    from nipype.interfaces.utility import Rename

    rename = Rename()
    rename.inputs.in_file = filename
    rename.inputs.format_string = entities + suffix

    return rename.run().outputs.out_file


def _get_pet_bids_components(cropped: bool, suvr_reference_region: str) -> str:
    """Return a string composed of the PET-specific entities (space, resolution,
    desc, and suvr), suffix, and extension.
    """
    space = "_space-MNI152NLin2009cSym"
    resolution = "_res-1x1x1"
    desc = "_desc-Crop" if cropped else ""
    suvr = f"_suvr-{suvr_reference_region}"

    return f"{space}{desc}{resolution}{suvr}_pet.nii.gz"


def print_end_pipeline(pet, final_file):
    """
    Display end message for <subject_id> when <final_file> is connected.
    """
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(pet))
