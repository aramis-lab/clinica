def get_substitutions_datasink(bids_file, pipeline_name):
    suffix = "T1w" if pipeline_name == "t1-linear" else "flair"
    substitutions = [
        (
            f"{bids_file}Warped_cropped.nii.gz",
            f"{bids_file}_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_{suffix}.nii.gz",
        ),
        (
            f"{bids_file}0GenericAffine.mat",
            f"{bids_file}_space-MNI152NLin2009cSym_res-1x1x1_affine.mat",
        ),
        (
            f"{bids_file}Warped.nii.gz",
            f"{bids_file}_space-MNI152NLin2009cSym_res-1x1x1_{suffix}.nii.gz",
        ),
    ]
    return bids_file, substitutions


# Function used by the nipype interface.
# It crops an image based on the reference.
def crop_nifti(input_img, ref_crop):
    """Crop input image based on the reference.

    It uses nilearn `resample_to_img` function.

    Args:
        input_img (str): image to be processed
        ref_crop (str): template used to crop the image

    Returns:
       output_img (NIfTI image): crop image on disk.
       crop_template: (NIfTI image): output template on disk.
    """
    import os

    from nilearn.image import crop_img, resample_to_img

    basedir = os.getcwd()
    # crop_ref = crop_img(ref_img, rtol=0.5)
    # crop_ref.to_filename(os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped_template.nii.gz'))
    # crop_template = os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped_template.nii.gz')

    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_crop, force_resample=True)
    crop_img.to_filename(
        os.path.join(
            basedir, os.path.basename(input_img).split(".nii")[0] + "_cropped.nii.gz"
        )
    )

    output_img = os.path.join(
        basedir, os.path.basename(input_img).split(".nii")[0] + "_cropped.nii.gz"
    )
    crop_template = ref_crop

    return output_img, crop_template


def print_end_pipeline(anat, final_file):
    """Display end message for <subject_id> when <final_file> is connected."""
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_end_image

    print_end_image(get_subject_id(anat))
