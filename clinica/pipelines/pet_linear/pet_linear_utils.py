# Functions used by nipype interface.

# Initiate the pipeline
def init_input_node(pet):
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(pet)
    print_begin_image(image_id)
    return pet


# Concatenate two transformation in one transformation list
def concatenate_transforms(pet_to_t1w_tranform, t1w_to_mni_tranform):
    """Concatenate two input transformation files into a list.
    Args:
       transform1 (str): first transformation to apply
       transform2 (str): second transformation to apply
    Returns:
       transform_list (list of string): both transform files path in a list
    """
    return [t1w_to_mni_tranform, pet_to_t1w_tranform]


# Normalize the images based on the reference mask region
def suvr_normalization(input_img, norm_img, ref_mask):
    """Normalize the input image according to the reference region.
    It uses nilearn `resample_to_img` and scipy `trim_mean` functions.
    This function is different than the one in other PET pipelines
    because there is a downsampling step.
    Args:
       input_img (str): image to be processed
       norm_img (str): image used to compute the mean of the reference region
       ref_mask (str): mask of the reference region
    Returns:
       output_img (nifty image): normalized nifty image
       mask_template (nifty image): output mask on disk
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


# It crops an image based on the reference.
def crop_nifti(input_img, ref_crop):
    """Crop input image based on the reference. It uses nilearn
    `resample_to_img` function.
    Args:
       input_img (str): image to be processed
       ref_img (str): template used to crop the image
    Returns:
       output_img (nifty image): crop image on disk.
       crop_template (nifty image): output template on disk.
    """

    import os

    import nibabel as nib
    import numpy as np
    from nilearn.image import resample_to_img

    basedir = os.getcwd()

    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_crop, force_resample=True)

    output_img = os.path.join(
        basedir, os.path.basename(input_img).split(".nii")[0] + "_cropped.nii.gz"
    )

    crop_img.to_filename(output_img)

    return output_img


def rename_into_caps(
    in_bids_pet,
    fname_pet,
    fname_trans,
    suvr_reference_region,
    uncropped_image,
    fname_pet_in_t1w=None,
):
    """
    Rename the outputs of the pipelines into CAPS format.
    Args:
        in_bids_pet (str): Input BIDS PET to extract the <source_file>
        fname_pet (str): Preprocessed PET file.
        fname_trans (str): Transformation file from PET to MRI space
        suvr_reference_region (str): SUVR mask name for file name output
        uncropped_image (bool): Pipeline argument for image cropping
        fname_pet_in_t1w (bool): Pipeline argument for saving intermediate file
    Returns:
        The different outputs in CAPS format
    """
    import os

    from nipype.interfaces.utility import Rename
    from nipype.utils.filemanip import split_filename

    _, source_file_pet, _ = split_filename(in_bids_pet)

    # Rename into CAPS PET:
    rename_pet = Rename()
    rename_pet.inputs.in_file = fname_pet
    if not uncropped_image:
        suffix = f"_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-{suvr_reference_region}_pet.nii.gz"
        rename_pet.inputs.format_string = source_file_pet + suffix
    else:
        suffix = f"_space-MNI152NLin2009cSym_res-1x1x1_suvr-{suvr_reference_region}_pet.nii.gz"
        rename_pet.inputs.format_string = source_file_pet + suffix
    out_caps_pet = rename_pet.run().outputs.out_file

    # Rename into CAPS transformation file:
    rename_trans = Rename()
    rename_trans.inputs.in_file = fname_trans
    rename_trans.inputs.format_string = source_file_pet + "_space-T1w_rigid.mat"
    out_caps_trans = rename_trans.run().outputs.out_file

    # Rename intermediate PET in T1w MRI space
    if fname_pet_in_t1w is not None:
        rename_pet_in_t1w = Rename()
        rename_pet_in_t1w.inputs.in_file = fname_pet_in_t1w
        rename_pet_in_t1w.inputs.format_string = (
            source_file_pet + "_space-T1w_pet.nii.gz"
        )
        out_caps_pet_in_t1w = rename_pet_in_t1w.run().outputs.out_file
    else:
        out_caps_pet_in_t1w = None

    return out_caps_pet, out_caps_trans, out_caps_pet_in_t1w


def print_end_pipeline(pet, final_file):
    """
    Display end message for <subject_id> when <final_file> is connected.
    """
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(pet))
