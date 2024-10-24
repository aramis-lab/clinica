from pathlib import Path
from typing import Optional, Tuple, Union

from clinica.utils.pet import SUVRReferenceRegion

_all__ = [
    "init_input_node",
    "concatenate_transforms",
    "perform_suvr_normalization",
    "rename_into_caps",
    "print_end_pipeline",
]


def init_input_node(pet: str) -> str:
    """Initiate the pipeline."""
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

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


def perform_suvr_normalization(
    pet_image_path: Path,
    normalizing_image_path: Path,
    reference_mask_path: Path,
) -> Path:
    """Normalize the input image according to the reference region.

    It uses nilearn `resample_to_img` and scipy `trim_mean` functions.
    This function is different than the one in other PET pipelines
    because there is a down-sampling step.

    Parameters
    ----------
    pet_image_path : Path
        The path to the image to be processed.

    normalizing_image_path : Path
        The path to the image used to compute the mean of the reference region.

    reference_mask_path : Path
        The path to the mask of the reference region.

    Returns
    -------
    output_img : Path
        The path to the normalized nifti image.
    """
    import nibabel as nib
    import numpy as np
    from nilearn.image import resample_to_img
    from scipy.stats import trim_mean

    from clinica.utils.filemanip import get_filename_no_ext
    from clinica.utils.stream import cprint

    pet_image = nib.load(pet_image_path)
    normalizing_image = nib.load(normalizing_image_path)
    reference_mask = nib.load(reference_mask_path)

    # Down-sample the pet image used for normalization so we can multiply it with the mask
    down_sampled_image = resample_to_img(
        normalizing_image, reference_mask, interpolation="nearest"
    )

    # Compute the mean of the region
    region = np.multiply(
        down_sampled_image.get_fdata(), reference_mask.get_fdata(dtype="float32")
    )
    array_region = np.where(region != 0, region, np.nan).flatten()
    region_mean = trim_mean(array_region[~np.isnan(array_region)], 0.1)
    cprint(region_mean, lvl="info")

    # Divide the value of the image voxels by the computed mean
    data = pet_image.get_fdata(dtype="float32") / region_mean

    # Create and save the normalized image
    output_image = (
        Path.cwd() / f"{get_filename_no_ext(pet_image_path)}_suvr_normalized.nii.gz"
    )
    normalized_img = nib.Nifti1Image(data, pet_image.affine, header=pet_image.header)
    normalized_img.to_filename(output_image)

    return output_image


def remove_mni_background(
    pet_image_path: Path,
    mni_mask_path: Path,
) -> Path:
    """ """
    import nibabel as nib

    from clinica.utils.filemanip import get_filename_no_ext

    pet_image = nib.load(pet_image_path)
    mni_mask = nib.load(mni_mask_path)

    # TODO: check if it is the best way to do the operation
    data = pet_image.get_fdata(dtype="float32") * mni_mask.get_fdata(dtype="float32")

    output_image = (
        Path.cwd() / f"{get_filename_no_ext(pet_image_path)}_remove_background.nii.gz"
    )
    normalized_img = nib.Nifti1Image(data, pet_image.affine, header=pet_image.header)
    normalized_img.to_filename(output_image)

    return output_image


def rename_into_caps(
    pet_bids_image_filename: Path,
    pet_preprocessed_image_filename: Path,
    pet_to_mri_transformation_filename: Path,
    suvr_reference_region: Union[str, SUVRReferenceRegion],
    uncropped_image: bool,
    pet_filename_in_t1w_raw: Optional[Path] = None,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path, Optional[Path]]:
    """
    Rename the outputs of the pipelines into CAPS format.

    More precisely, the following files will be renamed:

        - 'pet_filename_raw' will be renamed to 'pet_filename_caps'
        - 'transformation_filename_raw' will be renamed to 'transformation_filename_caps'
        - 'pet_filename_in_t1w_raw' will be renamed (if provided) to 'pet_filename_in_t1w_caps'

    Parameters
    ----------
    pet_bids_image_filename : Path
        Input PET image file from the BIDS dataset.
        This file is used to extract the entities from the source
        file like subject, session, run, tracer...

    pet_preprocessed_image_filename : Path
        Preprocessed PET file as outputted by Nipype.

    pet_to_mri_transformation_filename : Path
        Transformation file from PET to MRI space as outputted by Nipype.

    suvr_reference_region : str or SUVRReferenceRegion
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
    pet_filename_caps : Path
        The renamed preprocessed PET file to match CAPS conventions.

    transformation_filename_caps : Path
        The transformation file from PET to MRI space renamed to match CAPS conventions.

    pet_filename_in_t1w_caps : Path or None
        Intermediate PET in T1w MRI space renamed to match CAPS conventions.
        If 'pet_filename_in_t1w_raw' is None, this will be None.
    """
    bids_entities = _get_bids_entities_without_suffix(
        pet_bids_image_filename, suffix="pet"
    )
    if output_dir:
        bids_entities = str(output_dir / bids_entities)
    pet_filename_caps = _rename_pet_into_caps(
        bids_entities,
        pet_preprocessed_image_filename,
        not uncropped_image,
        SUVRReferenceRegion(suvr_reference_region),
    )
    transformation_filename_caps = _rename_transformation_into_caps(
        bids_entities, pet_to_mri_transformation_filename
    )
    pet_filename_in_t1w_caps = None
    if pet_filename_in_t1w_raw is not None:
        pet_filename_in_t1w_caps = _rename_intermediate_pet_in_t1w_space_into_caps(
            bids_entities, pet_filename_in_t1w_raw
        )

    return pet_filename_caps, transformation_filename_caps, pet_filename_in_t1w_caps


def _get_bids_entities_without_suffix(filename: Path, suffix: str) -> str:
    """Return the BIDS entities without the suffix from a BIDS path."""
    from clinica.utils.filemanip import get_filename_no_ext

    return get_filename_no_ext(filename).rstrip(f"_{suffix}")


def _rename_pet_into_caps(
    entities: str, filename: Path, cropped: bool, region: SUVRReferenceRegion
) -> Path:
    """Rename into CAPS PET."""
    return _rename(filename, entities, _get_pet_bids_components(cropped, region))


def _rename_transformation_into_caps(entities: str, filename: Path) -> Path:
    """Rename into CAPS transformation file."""
    return _rename(filename, entities, "_space-T1w_rigid.mat")


def _rename_intermediate_pet_in_t1w_space_into_caps(
    entities: str, filename: Path
) -> Path:
    """Rename intermediate PET in T1w MRI space."""
    return _rename(filename, entities, "_space-T1w_pet.nii.gz")


def _rename(filename: Path, entities: str, suffix: str) -> Path:
    """Rename 'filename' into '{entities}{suffix}'."""
    from nipype.interfaces.utility import Rename

    rename = Rename()
    rename.inputs.in_file = str(filename)
    rename.inputs.format_string = f"{entities}{suffix}"

    return Path(rename.run().outputs.out_file)


def _get_pet_bids_components(cropped: bool, region: SUVRReferenceRegion) -> str:
    """Return a string composed of the PET-specific entities (space, resolution,
    desc, and suvr), suffix, and extension.
    """
    space = "_space-MNI152NLin2009cSym"
    resolution = "_res-1x1x1"
    desc = "_desc-Crop" if cropped else ""
    suvr = f"_suvr-{region.value}"

    return f"{space}{desc}{resolution}{suvr}_pet.nii.gz"


def print_end_pipeline(pet: str, final_file):
    """
    Display end message for <subject_id> when <final_file> is connected.
    """
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(pet))
