# coding: utf8


def get_substitutions_datasink(bids_file):

    substitutions_ls = [  # registration
            (bids_file + '_T1w_corrected.nii.gz',
                bids_file + '_desc-BiasCorrected_T1w.nii.gz'),
            (bids_file + 'Warped_cropped_intensity_norm.nii.gz',
                bids_file + '_space-MNI152NLin2009cSym_res-1x1x1_intensity_norm_T1w.nii.gz'),
            (bids_file + 'Warped_cropped.nii.gz',
                bids_file + '_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_T1w.nii.gz'),
            (bids_file + '0GenericAffine.mat',
                bids_file + '_space-MNI152NLin2009cSym_res-1x1x1_affine.mat'),
            (bids_file + 'Warped.nii.gz',
                bids_file + '_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz')
            ]
    return bids_file, substitutions_ls


# Function used by the nipype interface.
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
    from nilearn.image import resample_img, crop_img, resample_to_img
    from nibabel.spatialimages import SpatialImage

    basedir = os.getcwd()
    # crop_ref = crop_img(ref_img, rtol=0.5)
    # crop_ref.to_filename(os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped_template.nii.gz'))
    # crop_template = os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped_template.nii.gz')

    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_crop, force_resample=True)
    crop_img.to_filename(os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped.nii.gz'))

    output_img = os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped.nii.gz')
    crop_template = ref_crop

    return output_img, crop_template


def print_end_pipeline(t1w, final_file):
    """
    Display end message for <subject_id> when <final_file> is connected.
    """
    from clinica.utils.ux import print_end_image
    from clinica.utils.filemanip import get_subject_id
    print_end_image(get_subject_id(t1w))
