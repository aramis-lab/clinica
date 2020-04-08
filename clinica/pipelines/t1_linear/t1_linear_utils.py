# coding: utf8

"""T1 Linear - Clinica Utilities.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


# Get containers to produce the CAPS structure
# Warning!!! This file should be in the future in the utils package of Clinica
def container_from_filename(bids_or_caps_filename):
    """Extract container from BIDS or CAPS file.
    Args:
       bids_or_caps_filename (str): full path to BIDS or CAPS filename.
    Returns:
       Container path of the form "subjects/<participant_id>/<session_id>"
    Examples:
       >>> from clinica.utils.nipype import container_from_filename
       >>> container_from_filename('/path/to/bids/sub-CLNC01/ses-M00/anat/sub-CLNC01_ses-M00_T1w.nii.gz')
               'subjects/sub-CLNC01/ses-M00'
       >>> container_from_filename('caps/subjects/sub-CLNC01/ses-M00/dwi/preprocessing/sub-CLNC01_ses-M00_preproc.nii')
               'subjects/sub-CLNC01/ses-M00'
    """

    import os
    import re
    m = re.search(r'(sub-[a-zA-Z0-9]+)/(ses-[a-zA-Z0-9]+)', bids_or_caps_filename)
    if m is None:
        raise ValueError(
                'Input filename is not in a BIDS or CAPS compliant format.'
                'It does not contain the participant and session ID.'
                )
    subject = m.group(1)
    session = m.group(2)
    return os.path.join('subjects', subject, session)


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
