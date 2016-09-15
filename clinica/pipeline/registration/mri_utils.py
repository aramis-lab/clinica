#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This module contains functions used for the registration pipeline or other pipelines."""

def convert_flirt_transformation_to_mrtrix_transformation(in_souce_image, in_reference_image, in_flirt_matrix, name_output_matrix=None):
    """
    Convert flirt matrix to mrtrix matrix.

    This function converts converts a transformation matrix produced by FSL's
    flirt command into a format usable by MRtrix. The output of this function
    is usually for the mrtransform command.

    Args:
        in_souce_image (str): File containing the source image used in FSL
            flirt with the -in flag.
        in_reference_image (str): File containing the reference image used in
            FSL flirt with the -ref flag.
        in_flirt_matrix (str): File containing the transformation matrix
            obtained by FSL flirt.
        name_output_matrix (Optional[str]): Name of the output matrix
            (default=deformed_image.nii.gz).

    Returns:
        out_mrtrix_matrix (str): Transformation matrix in MRtrix format.
    """
    import os.path as op
    import os

    assert(op.isfile(in_souce_image))
    assert(op.isfile(in_reference_image))
    assert(op.isfile(in_flirt_matrix))

    if name_output_matrix is None:
        out_mrtrix_matrix = op.abspath('mrtrix_matrix.mat')
    else:
        out_mrtrix_matrix = op.abspath(name_output_matrix);

    cmd = 'transformconvert ' + in_flirt_matrix + ' ' + in_souce_image + ' ' + in_reference_image + ' flirt_import ' + out_mrtrix_matrix
    os.system(cmd)

    return out_mrtrix_matrix



def apply_mrtrix_transform_without_resampling(in_image, in_mrtrix_matrix, name_output_image=None):
    """
    Apply a transformation without resampling.

    This function applies a linear transform on the input image without
    reslicing: it only modifies the transform matrix in the image header.

    Args:
        in_image (str): File containing the input image to be transformed.
        in_mrtrix_matrix (str): File containing the transformation matrix
            obtained by the MRtrix transformconvert command.
        name_output_image (Optional[str]): Name of the output image
            (default=deformed_image.nii.gz).

    Returns:
        out_deformed_image (str): File containing the deformed image according
            to in_mrtrix_matrix transformation.

    Example:
        >>> from clinica.pipeline.registration.mri_utils import apply_mrtrix_transform_without_resampling
        >>> apply_mrtrix_transform_without_resampling('t1_image.nii', 't1_to_diffusion.txt', 't1_in_diffusion_space.nii')
    """
    import os.path as op
    import os

    assert(op.isfile(in_image))
    assert(op.isfile(in_mrtrix_matrix))

    if name_output_image is None:
        out_deformed_image = op.abspath('deformed_image.nii.gz')
    else:
        out_deformed_image = op.abspath(name_output_image);

    cmd = 'mrtransform -linear ' + in_mrtrix_matrix + ' ' + in_image + ' ' + out_deformed_image
    os.system(cmd)

    return out_deformed_image
