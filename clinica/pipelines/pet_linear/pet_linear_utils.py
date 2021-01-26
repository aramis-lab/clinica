# coding: utf8

"""pet_linear - Clinica Utilities.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


# Function used by nipype interface.
# Normalize the images based on the reference mask region
def suvr_normalization(input_img, ref_mask):
    """Normalize the input image according to the reference region.
    It uses nilearn `resample_to_img` function.
    Args:
       input_img (str): image to be processed
       ref_mask (str): mask of the reference region
    Returns:
       output_img (nifty image): normalized nifty image
    """

    import nibabel as nib
    import os
    import numpy as np
    from nilearn.image import resample_to_img

    basedir = os.getcwd()

    # Downsample the input image so we can multiply it with the mask
    ds_img = resample_to_img(img, mask, interpolation="nearest")

    # Compute the mean of the region
    region = np.multiply(ds_img.get_data(), mask.get_data())
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    # Divide the value of the image voxels by the computed mean
    data = img.get_data() / region_mean
    normalized_img = nib.Nifti1Image(data, img.affine, header=img.header)

    output_img = os.path.join(
       basedir,
       os.path.basename(input_img).split('.nii')[0] + '_suvr_normalized.nii.gz')

    normalized_img.to_filename(output_img)
    mask_template = ref_mask

    return output_img, mask_template


# It crops an image based on the reference.
def crop_nifti(input_img, ref_crop):
    """Crop input image based on the reference. It uses nilearn
    `resample_to_img` function.
    Args:
       input_img (str): image to be processed
       ref_img (str): template used to crop the image
    Returns:
       output_img (nifty image): crop image on disk.
       crop_template: (nifty image): output template on disk.
    """

    import nibabel as nib
    import os
    import numpy as np
    from nilearn.image import resample_to_img

    basedir = os.getcwd()

    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_crop)

    output_img = os.path.join(
       basedir,
       os.path.basename(input_img).split('.nii')[0] + '_cropped.nii.gz')

    crop_img.to_filename(output_img)
    crop_template = ref_crop

    return output_img, crop_template


def print_end_pipeline(pet, final_file):
    """
    Display end message for <subject_id> when <final_file> is connected.
    """
    from clinica.utils.ux import print_end_image
    from clinica.utils.filemanip import get_subject_id
    print_end_image(get_subject_id(pet))
