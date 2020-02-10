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


def get_data_datasink(image_id):
    substitutions_ls = [  # registration
            (image_id + '_T1w_corrected.nii.gz',
                image_id + '_corrected_T1w.nii.gz'),
            (image_id + 'Warped_cropped_intensity_norm.nii.gz',
                image_id + '_space-MNI152NLin2009cSym_res-1x1x1_intensity_norm_T1w.nii.gz'),
            (image_id + 'Warped_cropped.nii.gz',
                image_id + '_space-MNI152NLin2009cSym_res-1x1x1_T1w.nii.gz'),
            (image_id + 'Warped_cropped.pt',
                image_id + '_space-MNI152NLin2009cSym_res-1x1x1_T1w.pt'),
            (image_id + 'Warped.nii.gz',
                image_id + '_space-MNI152NLin2009cSym_res-1x1x1_linear_registration_T1w.nii.gz')
            ]
    return image_id, substitutions_ls


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
    crop_img = resample_to_img(input_img, ref_crop)
    crop_img.to_filename(os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped.nii.gz'))

    output_img = os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_cropped.nii.gz')

    return output_img, crop_template


def ants_histogram_intensity_normalization(crop_template, input_img, image_dimension):
    """Histogram-based intensity normalization. It uses the `ImageMath` binary
    with the parameter `HistogramMatch` provieded by the ANTS package.
    This is a function to do histogram-based intensity normalization.
    Normalize the grayscale values for a source image by matching the shape of
    the source image histogram to a reference histogram.
    Args:
       crop_template (str): reference histogram obtained from this image.
       input_img (str): source image (histogram to be normalized).
       image_dimension (int): 2 or 3, for the image channels.
    Returns:
       output_img (Nifty image): image with histogram normalized.
    """

    import os

    basedir = os.getcwd()
    output_img = os.path.join(basedir, os.path.basename(input_img).split('.nii')[0] + '_intensity_norm.nii.gz')
    commands = [
            'ImageMath',
            str(image_dimension),
            output_img,
            'HistogramMatch',
            input_img,
            crop_template
            ]
    cmd = ' '.join(commands)
    os.system(cmd)

    return output_img
